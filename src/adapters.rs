//! Finding library adapters and primers anywhere in a read, by shared 12-mers.

use std::collections::HashMap;
use std::fs;
use std::path::Path;

use anyhow::{ensure, Context, Result};

use crate::read::{ASCII_TO_CODE, BASE_N};

/// k-mer length: 4^12 k-mers fit a 2 MiB bit table, so the scan is one bit test per base.
const K: usize = 12;

/// Hits on alignment diagonals this close belong to one occurrence (sequencing indels
/// shift the diagonal).
const DIAGONAL_SLACK: i64 = 8;

/// An occurrence needs at least this many k-mer hits on one diagonal...
const MIN_HITS: usize = 2;

/// ...and at least this share of the adapter's indexed k-mers. Repeats such as (GA)n share a
/// few k-mers with adapters; a real adapter, even at Q10, keeps a quarter of its 12-mers.
const MIN_SHARE: f64 = 0.25;

/// What an adapter or primer is, for checks that care which one was found.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Role {
    Library,
    /// ONT cDNA strand-switching primer, at the 5' end of full-length transcripts.
    CdnaSsp,
    /// ONT cDNA VN primer, at the poly(A) end of full-length transcripts.
    CdnaVnp,
    /// From `--adapters`, `--kinnex-adapters` or a barcode list.
    Custom,
}

/// Built-in adapters and primers. ONT's are dorado's (`adapter_primer_kits.cpp`,
/// `barcode_kits.cpp`); PacBio's are the SMRTbell hairpin and the Illumina-style ends of the
/// Twist UDI adapters PacBio lists for Amplifi (`Amplifi_TwistUDIadapters_noP7P5.fasta`).
const BUILT_IN: [(&str, &str, Role); 11] = [
    (
        "ONT ligation adapter (LSK110/114)",
        "CCTGTACTTCGTTCAGTTACGTATTGC",
        Role::Library,
    ),
    (
        "ONT rapid adapter (RAD/RBK)",
        "GTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        Role::Library,
    ),
    ("ONT RNA004 adapter", "GGTTGTTTCTGTTGGTGCTG", Role::Library),
    ("ONT MAB114 primer", "CCATATCCGTGTCGCCCTT", Role::Library),
    (
        "ONT DRB004 barcode adapter",
        "CCAATTTGTGGGTTCGT",
        Role::Library,
    ),
    (
        "ONT cDNA SSP primer",
        "TTTCTGTTGGTGCTGATATTGCTGGG",
        Role::CdnaSsp,
    ),
    (
        "ONT cDNA VNP primer",
        "ACTTGCCTGTCGCTCTATCTTCTTT",
        Role::CdnaVnp,
    ),
    (
        "ONT PCS110 VNP primer",
        "ACTTGCCTGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAGTTTT",
        Role::CdnaVnp,
    ),
    (
        "Illumina read 1 primer (10x, Twist UDI)",
        "CTACACGACGCTCTTCCGATCT",
        Role::Library,
    ),
    (
        "Illumina read 2 primer (Twist UDI)",
        "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT",
        Role::Library,
    ),
    (
        "PacBio SMRTbell adapter",
        "ATCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGAT",
        Role::Library,
    ),
];

/// One adapter occurrence in a read.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Hit {
    pub adapter: usize,
    /// First and one-past-last read position covered by matching k-mers.
    pub start: usize,
    pub end: usize,
}

/// Where an indexed k-mer lies in an adapter.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Owner {
    adapter: u16,
    strand: u8,
    offset: u16,
}

/// A k-mer of a read matching an adapter k-mer.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct Seed {
    adapter: u16,
    strand: u8,
    /// Read position minus adapter offset: constant along one occurrence.
    diagonal: i64,
    position: usize,
}

/// Buffers reused between scans.
#[derive(Debug, Default)]
pub struct ScanBuffer {
    seeds: Vec<Seed>,
    candidates: Vec<(usize, Hit)>,
}

/// Index of every adapter's k-mers, on both strands.
#[derive(Debug, Clone)]
pub struct Adapters {
    names: Vec<String>,
    present: Vec<u64>,
    owners: HashMap<u32, Vec<Owner>>,
    /// Hits an occurrence of each adapter needs.
    required: Vec<usize>,
    roles: Vec<Role>,
}

impl Adapters {
    /// The built-in adapters plus those in `extra`, a FASTA file.
    ///
    /// # Errors
    ///
    /// Fails when `extra` cannot be read or holds a sequence shorter than 12 bases.
    pub fn load(extra: Option<&Path>) -> Result<Self> {
        let mut list: Vec<(String, Vec<u8>, Role)> = BUILT_IN
            .iter()
            .map(|&(name, seq, role)| (name.to_owned(), seq.as_bytes().to_vec(), role))
            .collect();
        if let Some(path) = extra {
            for (name, seq) in read_fasta(path)? {
                ensure!(seq.len() >= K, "adapter {name} is shorter than {K} bases");
                list.push((name, seq, Role::Custom));
            }
        }
        Ok(Self::new(
            list.into_iter()
                .map(|(name, seq, role)| (name, seq, role, MIN_SHARE))
                .collect(),
        ))
    }

    /// An index of sequences given by name, such as barcodes, where an occurrence needs
    /// `min_share` of a sequence's informative k-mers; sequences shorter than 12 bases are
    /// left out.
    pub fn from_sequences(list: Vec<(String, Vec<u8>)>, min_share: f64) -> Self {
        Self::with_shares(
            list.into_iter()
                .map(|(name, seq)| (name, seq, min_share))
                .collect(),
        )
    }

    /// Like [`Adapters::from_sequences`], with a share for each sequence.
    pub fn with_shares(list: Vec<(String, Vec<u8>, f64)>) -> Self {
        Self::new(
            list.into_iter()
                .filter(|(_, seq, _)| seq.len() >= K)
                .map(|(name, seq, share)| (name, seq, Role::Custom, share))
                .collect(),
        )
    }

    fn new(list: Vec<(String, Vec<u8>, Role, f64)>) -> Self {
        let mut index = Self {
            names: Vec::with_capacity(list.len()),
            present: vec![0; (1 << (2 * K)) / 64],
            owners: HashMap::new(),
            required: Vec::new(),
            roles: Vec::new(),
        };
        for (id, (name, seq, role, min_share)) in list.into_iter().enumerate() {
            index.roles.push(role);
            let id = u16::try_from(id).expect("fewer than 65536 adapters");
            let codes: Vec<u8> = seq.iter().map(|&b| ASCII_TO_CODE[usize::from(b)]).collect();
            let reverse: Vec<u8> = codes
                .iter()
                .rev()
                .map(|&c| if c == BASE_N { c } else { 3 - c })
                .collect();
            let mut indexed = 0;
            for (strand, sequence) in [&codes, &reverse].into_iter().enumerate() {
                for (offset, window) in sequence.windows(K).enumerate() {
                    let Some(kmer) = encode(window) else { continue };
                    if is_low_complexity(window) {
                        continue;
                    }
                    indexed += 1;
                    index.present[kmer as usize / 64] |= 1 << (kmer % 64);
                    index.owners.entry(kmer).or_default().push(Owner {
                        adapter: id,
                        strand: u8::try_from(strand).expect("two strands"),
                        offset: u16::try_from(offset).expect("adapters are short"),
                    });
                }
            }
            // `indexed` counts both strands; an occurrence lies on one.
            let per_strand = indexed as f64 / 2.0;
            index
                .required
                .push(MIN_HITS.max((per_strand * min_share).ceil() as usize));
            index.names.push(name);
        }
        index
    }

    /// What adapter `id` is.
    pub fn role(&self, id: usize) -> Role {
        self.roles[id]
    }

    /// Adapter names, in index order.
    pub fn names(&self) -> &[String] {
        &self.names
    }

    /// Finds adapter occurrences in `bases`, replacing the contents of `hits`.
    ///
    /// Matching k-mers are grouped by adapter, strand and alignment diagonal; a group is an
    /// occurrence when it has enough distinct hits. Where occurrences of different adapters
    /// overlap (adapters share sequence), the best supported is kept.
    pub fn scan(&self, bases: &[u8], buffer: &mut ScanBuffer, hits: &mut Vec<Hit>) {
        hits.clear();
        let seeds = &mut buffer.seeds;
        seeds.clear();
        let mask = (1u32 << (2 * K)) - 1;
        let mut kmer = 0u32;
        let mut valid = 0;
        for (i, &base) in bases.iter().enumerate() {
            if base == BASE_N {
                valid = 0;
                continue;
            }
            kmer = ((kmer << 2) | u32::from(base)) & mask;
            valid += 1;
            if valid >= K && self.present[kmer as usize / 64] & (1 << (kmer % 64)) != 0 {
                if let Some(owners) = self.owners.get(&kmer) {
                    let position = i + 1 - K;
                    seeds.extend(owners.iter().map(|o| Seed {
                        adapter: o.adapter,
                        strand: o.strand,
                        diagonal: position as i64 - i64::from(o.offset),
                        position,
                    }));
                }
            }
        }
        if seeds.is_empty() {
            return;
        }
        seeds.sort_unstable();
        let candidates = &mut buffer.candidates;
        candidates.clear();
        for group in seeds.chunk_by(|a, b| {
            a.adapter == b.adapter
                && a.strand == b.strand
                && b.diagonal - a.diagonal <= DIAGONAL_SLACK
        }) {
            let adapter = usize::from(group[0].adapter);
            let mut positions: Vec<usize> = group.iter().map(|s| s.position).collect();
            positions.sort_unstable();
            positions.dedup();
            if positions.len() >= self.required[adapter] {
                let hit = Hit {
                    adapter,
                    start: positions[0],
                    end: positions[positions.len() - 1] + K,
                };
                candidates.push((positions.len(), hit));
            }
        }
        candidates.sort_by(|a, b| b.0.cmp(&a.0).then(a.1.adapter.cmp(&b.1.adapter)));
        for &(_, hit) in candidates.iter() {
            if !hits.iter().any(|h| h.start < hit.end && hit.start < h.end) {
                hits.push(hit);
            }
        }
        hits.sort_by_key(|h| h.start);
    }
}

/// 2-bit packs a k-mer of base codes, or `None` if it holds an N.
fn encode(window: &[u8]) -> Option<u32> {
    window.iter().try_fold(0u32, |acc, &c| {
        (c != BASE_N).then_some((acc << 2) | u32::from(c))
    })
}

/// Homopolymer- or repeat-rich k-mers (`TTTTTTTTCCTG`, `GAGGAGGAGGAG`) match genomic
/// repeats by the thousand, so they are left out of the index.
fn is_low_complexity(window: &[u8]) -> bool {
    let mut counts = [0usize; 4];
    for &c in window {
        counts[usize::from(c.min(3))] += 1;
    }
    let mut dimers = [false; 16];
    for pair in window.windows(2) {
        dimers[usize::from(pair[0].min(3) * 4 + pair[1].min(3))] = true;
    }
    counts.iter().any(|&n| n >= 7) || dimers.iter().filter(|&&d| d).count() < 5
}

/// Reads a small FASTA file (adapters, barcodes).
///
/// # Errors
///
/// Fails when the file cannot be read.
pub fn read_fasta(path: &Path) -> Result<Vec<(String, Vec<u8>)>> {
    let text = fs::read_to_string(path).with_context(|| format!("reading {}", path.display()))?;
    Ok(parse_fasta(&text))
}

/// Minimal FASTA parsing for small adapter files.
pub fn parse_fasta(text: &str) -> Vec<(String, Vec<u8>)> {
    let mut records = Vec::new();
    for line in text.lines().map(str::trim) {
        if let Some(name) = line.strip_prefix('>') {
            let name = name
                .split_whitespace()
                .next()
                .unwrap_or_default()
                .to_owned();
            records.push((name, Vec::new()));
        } else if let Some((_, seq)) = records.last_mut() {
            seq.extend(line.bytes().map(|b| b.to_ascii_uppercase()));
        }
    }
    records
}

#[cfg(test)]
mod tests {
    use super::*;

    fn codes(seq: &str) -> Vec<u8> {
        seq.bytes().map(|b| ASCII_TO_CODE[usize::from(b)]).collect()
    }

    fn scan(seq: &str) -> Vec<Hit> {
        let adapters = Adapters::load(None).expect("built-ins load");
        let mut hits = Vec::new();
        adapters.scan(&codes(seq), &mut ScanBuffer::default(), &mut hits);
        hits
    }

    #[test]
    fn finds_adapters_on_both_strands() {
        let genomic = "GATTACAGGCTTAGCCATGCAAGTCCGATTGCAAGT";
        let read = format!("{genomic}CCTGTACTTCGTTCAGTTACGTATTGC{genomic}");
        let hits = scan(&read);
        assert_eq!(hits.len(), 1, "{hits:?}");
        assert_eq!(hits[0].adapter, 0);
        assert_eq!(hits[0].start, genomic.len());
        assert_eq!(hits[0].end, genomic.len() + 27);

        let reverse = "GCAATACGTAACTGAACGAAGTACAGG";
        let hits = scan(&format!("{genomic}{reverse}"));
        assert!(hits.iter().any(|h| h.adapter == 0), "{hits:?}");
    }

    #[test]
    fn ignores_simple_repeats() {
        let repeat = "GA".repeat(200) + &"AATGG".repeat(80) + &"CT".repeat(200);
        assert!(scan(&repeat).is_empty(), "{:?}", scan(&repeat));
    }

    #[test]
    fn finds_adapters_with_errors() {
        // The SMRTbell adapter with a substitution and a deletion.
        let adapter = "ATCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGAT";
        let mut damaged = adapter.as_bytes().to_vec();
        damaged[28] = b'T';
        damaged.remove(38);
        let genomic = "GATTACAGGCTTAGCCATGCAAGTCCGATTGCAAGTCATGCATTGACC";
        let read = format!(
            "{genomic}{}{genomic}",
            String::from_utf8(damaged).expect("ASCII")
        );
        let hits = scan(&read);
        assert_eq!(hits.len(), 1, "{hits:?}");
        let adapters = Adapters::load(None).expect("built-ins load");
        assert_eq!(adapters.names()[hits[0].adapter], "PacBio SMRTbell adapter");
    }

    #[test]
    fn ignores_plain_sequence_and_homopolymers() {
        assert!(scan("GATTACAGGCTTAGCCATGCAAGTCCGATTGCAAGTTTTTTTTTTTTTTTTTTTTTTTTT").is_empty());
    }

    #[test]
    fn parses_extra_adapters() {
        let records = parse_fasta(">my_adapter desc\nacgtacgtacgt\nACGT\n>two\nAAAA\n");
        assert_eq!(
            records[0],
            ("my_adapter".to_owned(), b"ACGTACGTACGTACGT".to_vec())
        );
        assert_eq!(records.len(), 2);
    }
}
