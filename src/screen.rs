//! Screening reads against built-in references: spike-in controls, common contaminants,
//! the human mitochondrion and the human nuclear genome.
//!
//! References and reads are sampled the same way: the canonical 19-mers whose hash falls in
//! the lowest 1/24 of the range (FracMinHash). The human genome is too large for that and is
//! sampled at 1 in 1,000, a subset of the same hashes. A read is assigned to the reference
//! that holds the largest share of its sampled k-mers (at that reference's density), when
//! that is at least 3 k-mers and at least half the share a read of its Q keeps,
//! `(1 - error)^19`, and never under 8%.

use std::collections::HashMap;
use std::hash::{BuildHasherDefault, Hasher};
use std::io::{Read as _, Write as _};
use std::path::Path;

use anyhow::{bail, ensure, Context, Result};
use flate2::read::{GzDecoder, MultiGzDecoder};
use flate2::write::GzEncoder;
use serde::{Deserialize, Serialize};

use crate::read::{ASCII_TO_CODE, BASE_N};

/// The small references, as gzipped FASTA whose headers are `key\tkind\taccession\tname`;
/// sequences of several records are joined by `N`.
const REFERENCES: &[u8] = include_bytes!("../assets/screen.fa.gz");

/// The human nuclear genome's sampled k-mers (T2T-CHM13v2.0, without chrM), built by
/// `lrqc build-human-sketch` (see [`write_human_sketch`]).
const HUMAN_SKETCH: &[u8] = include_bytes!("../assets/human.sketch.gz");

/// Header of the human sketch file.
const HUMAN_MAGIC: &[u8; 8] = b"LRQCHS1\n";

/// Human k-mers are sampled at or below this hash: one in 1,000, 1.8 million distinct
/// informative k-mers. A Q20 read keeps about 0.7 of them per kb, so most reads of 4 kb and
/// longer can be called.
const HUMAN_MAX: u64 = u64::MAX / 1000;

/// A read with at least this many informative k-mers at human density can be tested against
/// a sparse reference: a Q20 human read keeps about 5 hits of 6, past [`MIN_HITS`] and the
/// share bar. That takes about 7 kb, so shorter reads say nothing either way.
const SPARSE_TESTABLE_KMERS: usize = 6;

/// Sampled human hashes are stored by their bits from 23 up: 32 bits of a 55-bit range. A
/// read's hash matches a different k-mer's stored bits with probability under 0.1%.
const HUMAN_SHIFT: u32 = 23;

/// k-mer length: 4^19 k-mers make chance matches against 11 Mb of references rare, while a
/// Q10 read still keeps an eighth of its 19-mers.
const K: usize = 19;

/// k-mers whose hash is at or below this are sampled: one in 24. Sampling by hash needs one
/// comparison per base, and the index of sampled reference k-mers (about 450,000) fits in the
/// CPU cache, so a read makes a lookup every 24 bases into memory that is mostly cached.
const SAMPLED_MAX: u64 = u64::MAX / 24;

/// A read needs this many sampled k-mers in one reference...
const MIN_HITS: u32 = 3;

/// ...and at least this share of all its sampled k-mers. A human read holding a 1 kb nuclear
/// copy of mitochondrial sequence (a NUMT) in 20 kb shares about 5%, below this.
const MIN_FRACTION: f64 = 0.08;

/// ...and at least this fraction of the share expected for its read Q. A read from the
/// reference keeps `(1 - error)^19` of its 19-mers: about 98% at Q30, 83% at Q20, 13% at Q10.
/// Related species share conserved genes: on HiFi reads from the PacBio 96-plex microbial
/// run, Salmonella and Klebsiella reads shared under 30% of their 19-mers with E. coli
/// K-12, where reads of K-12 itself shared 99%. Half the expected share separates the two
/// for accurate reads and falls back to [`MIN_FRACTION`] for noisy ones.
const EXPECTED_SHARE_FACTOR: f64 = 0.5;

/// The index records which references hold each k-mer as bits of a `u16`.
const MAX_REFERENCES: usize = 16;

/// What a reference is for.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum Kind {
    /// Spiked in on purpose: ONT's DNA control strand (lambda) and RNA calibration strand,
    /// ERCC and SIRV spike-ins.
    Control,
    Contaminant,
    /// Organelle DNA of the sample itself.
    Organelle,
    /// Cloning vectors, which also hold natural plasmid backbones (IncQ RSF1010, IncP RK2)
    /// that bacterial isolates carry, so reported but not graded.
    Vector,
}

/// One screened reference.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Reference {
    pub key: String,
    pub name: String,
    pub kind: Kind,
    pub accession: String,
    pub length: usize,
    /// Sampled sparsely (1 in 1,000, the human genome): only reads long enough to test (see
    /// [`Screen::sparse_testable`]) say anything about it, so its share is of those.
    #[serde(default)]
    pub sparse: bool,
}

/// k-mer hashes are already mixed, so the map uses them as they are.
#[derive(Debug, Default)]
struct IdentityHasher(u64);

impl Hasher for IdentityHasher {
    fn finish(&self) -> u64 {
        self.0
    }

    fn write(&mut self, bytes: &[u8]) {
        for &b in bytes {
            self.0 = (self.0 << 8) | u64::from(b);
        }
    }

    fn write_u64(&mut self, value: u64) {
        self.0 = value;
    }
}

/// Sampled k-mer hash to the references holding it, one bit each.
type KmerMap = HashMap<u64, u16, BuildHasherDefault<IdentityHasher>>;

/// The sampled k-mer index of all references.
#[derive(Debug)]
pub struct Screen {
    references: Vec<Reference>,
    index: KmerMap,
    /// The human sketch's stored bits, sorted, and the human reference's id.
    human: Vec<u32>,
    human_id: usize,
}

/// Buffers reused between reads.
#[derive(Debug, Default)]
pub struct ScreenBuffer {
    sampled: Vec<u64>,
    /// The sampled hashes at human density whose k-mers can tell human from other genomes
    /// (see [`HUMAN_MIN_TRINUCLEOTIDES`]).
    human: Vec<u64>,
}

impl Screen {
    /// Builds the index of the built-in references.
    ///
    /// # Errors
    ///
    /// Fails only if the embedded references are corrupt.
    pub fn load() -> Result<Self> {
        let mut text = String::new();
        GzDecoder::new(REFERENCES)
            .read_to_string(&mut text)
            .context("reading the built-in screen references")?;
        let mut references = Vec::new();
        let mut index = KmerMap::default();
        let mut buffer = ScreenBuffer::default();
        let mut lines = text.lines();
        while let Some(header) = lines.next() {
            let fields: Vec<&str> = header.trim_start_matches('>').split('\t').collect();
            let sequence = lines.next().context("reference without sequence")?;
            let kind = match fields.get(1).copied() {
                Some("control") => Kind::Control,
                Some("organelle") => Kind::Organelle,
                Some("vector") => Kind::Vector,
                _ => Kind::Contaminant,
            };
            let id = references.len();
            assert!(
                id < MAX_REFERENCES,
                "at most {MAX_REFERENCES} screen references"
            );
            references.push(Reference {
                key: fields.first().copied().unwrap_or_default().to_owned(),
                name: fields.get(3).copied().unwrap_or_default().to_owned(),
                kind,
                accession: fields.get(2).copied().unwrap_or_default().to_owned(),
                length: sequence.len(),
                sparse: false,
            });
            let codes: Vec<u8> = sequence
                .bytes()
                .map(|b| ASCII_TO_CODE[usize::from(b)])
                .collect();
            sample(&codes, &mut buffer);
            for &m in &buffer.sampled {
                *index.entry(m).or_insert(0) |= 1 << id;
            }
        }
        let (genome_length, human) = read_human_sketch(HUMAN_SKETCH)?;
        let human_id = references.len();
        references.push(Reference {
            key: "human".to_owned(),
            name: "Human nuclear genome (T2T-CHM13, sampled 1 in 1,000)".to_owned(),
            kind: Kind::Contaminant,
            accession: "GCF_009914755.1".to_owned(),
            length: genome_length,
            sparse: true,
        });
        assert!(
            references.len() <= MAX_REFERENCES,
            "at most {MAX_REFERENCES} screen references"
        );
        Ok(Self {
            references,
            index,
            human,
            human_id,
        })
    }

    pub fn references(&self) -> &[Reference] {
        &self.references
    }

    /// The reference a read comes from, if any; `read_q` raises the bar for accurate reads
    /// (see [`EXPECTED_SHARE_FACTOR`]).
    pub fn classify(
        &self,
        bases: &[u8],
        read_q: Option<f32>,
        buffer: &mut ScreenBuffer,
    ) -> Option<usize> {
        sample(bases, buffer);
        let total = buffer.sampled.len() as u64;
        if total == 0 {
            return None;
        }
        // Each reference's share is of the read's k-mers sampled at its density.
        let human_total = buffer.human.len() as u64;
        let sampled_at = |r: usize| {
            if r == self.human_id {
                human_total
            } else {
                total
            }
        };
        // Hits per reference, and hits on k-mers no other reference holds. E. coli K-12
        // carries lambda-like prophage remnants, so some k-mers belong to both.
        let mut hits = [0u64; MAX_REFERENCES];
        let mut unique = [0u64; MAX_REFERENCES];
        for &m in &buffer.sampled {
            let mut owners = self.index.get(&m).copied().unwrap_or(0);
            if m <= HUMAN_MAX && buffer.human.binary_search(&m).is_ok() && self.is_human(m) {
                owners |= 1 << self.human_id;
            }
            if owners == 0 {
                continue;
            }
            for (r, count) in hits.iter_mut().enumerate() {
                *count += u64::from(owners >> r & 1);
            }
            if owners.is_power_of_two() {
                unique[owners.trailing_zeros() as usize] += 1;
            }
        }
        // Shares compare as fractions, cross-multiplied: (hits, unique) over sampled.
        let share = |r: usize| (hits[r], unique[r], sampled_at(r).max(1));
        let greater = |a: usize, b: usize| {
            let ((ha, ua, sa), (hb, ub, sb)) = (share(a), share(b));
            (ha * sb, ua * sb).cmp(&(hb * sa, ub * sa))
        };
        // Only references with enough hits compete: a 1 kb read samples about one k-mer at
        // human density, and one chance hit there would be a 100% share.
        let candidates = || (0..self.references.len()).filter(|&r| hits[r] >= u64::from(MIN_HITS));
        let best = candidates().max_by(|&a, &b| greater(a, b).then(b.cmp(&a)))?;
        // A read whose hits all fall on k-mers another reference holds just as well is
        // ambiguous.
        let tied = candidates().any(|r| r != best && greater(r, best) == std::cmp::Ordering::Equal);
        let count = hits[best];
        let min_fraction = read_q.map_or(MIN_FRACTION, |q| {
            let error = 10f64.powf(-f64::from(q) / 10.0).min(1.0);
            let expected = (1.0 - error).powi(K as i32);
            (EXPECTED_SHARE_FACTOR * expected).max(MIN_FRACTION)
        });
        (!tied
            && count >= u64::from(MIN_HITS)
            && count as f64 >= min_fraction * sampled_at(best) as f64)
            .then_some(best)
    }

    /// Whether the read [`Screen::classify`] last saw was long enough to test against the
    /// sparse references.
    pub fn sparse_testable(buffer: &ScreenBuffer) -> bool {
        buffer.human.len() >= SPARSE_TESTABLE_KMERS
    }

    /// Whether a sampled hash at human density is in the human sketch.
    fn is_human(&self, hash: u64) -> bool {
        self.human.binary_search(&human_key(hash)).is_ok()
    }
}

/// Whether a FASTA header names the mitochondrion: NCBI's "mitochondrion", UCSC's `chrM`,
/// Ensembl's `MT`.
fn is_mitochondrial(header: &str) -> bool {
    let lower = header.to_ascii_lowercase();
    let name = lower.split_whitespace().next().unwrap_or_default();
    lower.contains("mitochondri") || matches!(name, "chrm" | "chrmt" | "mt" | "m")
}

/// The stored bits of a human-density hash.
fn human_key(hash: u64) -> u32 {
    u32::try_from(hash >> HUMAN_SHIFT).expect("hashes at human density fit 32 bits")
}

/// Reads a human sketch: the genome length and the sorted stored bits.
fn read_human_sketch(gz: &[u8]) -> Result<(usize, Vec<u32>)> {
    let mut bytes = Vec::new();
    GzDecoder::new(gz)
        .read_to_end(&mut bytes)
        .context("reading the built-in human sketch")?;
    ensure!(
        bytes.len() >= 20 && &bytes[..8] == HUMAN_MAGIC,
        "the human sketch is corrupt"
    );
    let genome_length = u64::from_le_bytes(bytes[8..16].try_into().expect("8 bytes"));
    let count = u32::from_le_bytes(bytes[16..20].try_into().expect("4 bytes")) as usize;
    let mut keys = Vec::with_capacity(count);
    let mut previous = 0u32;
    let mut rest = &bytes[20..];
    while !rest.is_empty() {
        let mut delta = 0u32;
        let mut shift = 0;
        loop {
            let (&byte, tail) = rest
                .split_first()
                .context("the human sketch is truncated")?;
            rest = tail;
            // The fifth byte of a u32 varint carries its top 4 bits.
            ensure!(shift < 28 || byte <= 0x0f, "the human sketch is corrupt");
            delta |= u32::from(byte & 0x7f) << shift;
            if byte & 0x80 == 0 {
                break;
            }
            shift += 7;
        }
        previous = previous
            .checked_add(delta)
            .context("the human sketch is corrupt")?;
        keys.push(previous);
    }
    ensure!(keys.len() == count, "the human sketch is truncated");
    Ok((usize::try_from(genome_length)?, keys))
}

/// Builds the human sketch from a FASTA of the genome (gzipped or not), leaving out records
/// whose header names the mitochondrion, which is screened on its own.
///
/// # Errors
///
/// Fails when the FASTA cannot be read or the sketch cannot be written.
pub fn write_human_sketch(fasta: &Path, out: &Path) -> Result<()> {
    let file =
        std::fs::File::open(fasta).with_context(|| format!("opening {}", fasta.display()))?;
    let mut text = String::new();
    let reader: Box<dyn std::io::Read> = if fasta.extension().is_some_and(|e| e == "gz") {
        Box::new(MultiGzDecoder::new(std::io::BufReader::new(file)))
    } else {
        Box::new(std::io::BufReader::new(file))
    };
    let mut keys: Vec<u32> = Vec::new();
    let mut genome_length = 0u64;
    let mut buffer = ScreenBuffer::default();
    let mut codes: Vec<u8> = Vec::new();
    let mut skip = false;
    let mut flush = |codes: &mut Vec<u8>, keys: &mut Vec<u32>, length: &mut u64| {
        if codes.is_empty() {
            return;
        }
        *length += codes.len() as u64;
        sample(codes, &mut buffer);
        keys.extend(buffer.human.iter().map(|&m| human_key(m)));
        codes.clear();
    };
    let mut lines = std::io::BufReader::new(reader);
    loop {
        text.clear();
        if std::io::BufRead::read_line(&mut lines, &mut text)? == 0 {
            break;
        }
        let line = text.trim_end();
        if let Some(header) = line.strip_prefix('>') {
            flush(&mut codes, &mut keys, &mut genome_length);
            skip = is_mitochondrial(header);
            eprintln!(
                "lrqc: {} {header}",
                if skip { "skipping" } else { "sketching" }
            );
        } else if !skip {
            codes.extend(line.bytes().map(|b| ASCII_TO_CODE[usize::from(b)]));
        }
    }
    flush(&mut codes, &mut keys, &mut genome_length);
    if keys.is_empty() {
        bail!("no sequence in {}", fasta.display());
    }
    keys.sort_unstable();
    keys.dedup();
    let mut bytes = Vec::with_capacity(keys.len() * 2 + 20);
    bytes.extend_from_slice(HUMAN_MAGIC);
    bytes.extend_from_slice(&genome_length.to_le_bytes());
    bytes.extend_from_slice(&u32::try_from(keys.len())?.to_le_bytes());
    let mut previous = 0u32;
    for &key in &keys {
        let mut delta = key - previous;
        previous = key;
        loop {
            let byte = (delta & 0x7f) as u8;
            delta >>= 7;
            if delta == 0 {
                bytes.push(byte);
                break;
            }
            bytes.push(byte | 0x80);
        }
    }
    let mut encoder = GzEncoder::new(
        std::fs::File::create(out).with_context(|| format!("creating {}", out.display()))?,
        flate2::Compression::best(),
    );
    encoder.write_all(&bytes)?;
    encoder.finish()?;
    eprintln!(
        "lrqc: {} sampled human k-mers from {genome_length} bases",
        keys.len()
    );
    Ok(())
}

/// Mixes a canonical k-mer so sampling is not biased toward low-complexity k-mers.
fn mix(mut x: u64) -> u64 {
    x ^= x >> 33;
    x = x.wrapping_mul(0xff51_afd7_ed55_8ccd);
    x ^= x >> 33;
    x = x.wrapping_mul(0xc4ce_b9fe_1a85_ec53);
    x ^ (x >> 33)
}

/// Fills `buffer.sampled` with the distinct sampled canonical k-mer hashes of `bases`.
fn sample(bases: &[u8], buffer: &mut ScreenBuffer) {
    buffer.human.clear();
    let out = &mut buffer.sampled;
    out.clear();
    let mask = (1u64 << (2 * K)) - 1;
    let shift = 2 * (K - 1);
    let (mut forward, mut reverse) = (0u64, 0u64);
    let mut valid = 0;
    for &base in bases {
        if base == BASE_N {
            valid = 0;
            continue;
        }
        let b = u64::from(base);
        forward = ((forward << 2) | b) & mask;
        reverse = (reverse >> 2) | ((3 - b) << shift);
        valid += 1;
        if valid >= K {
            let h = mix(forward.min(reverse));
            if h <= SAMPLED_MAX {
                let distinct = distinct_trinucleotides(forward);
                if distinct >= MIN_TRINUCLEOTIDES {
                    out.push(h);
                    if h <= HUMAN_MAX
                        && distinct >= HUMAN_MIN_TRINUCLEOTIDES
                        && HUMAN_GC_BASES.contains(&gc_bases(forward))
                    {
                        buffer.human.push(h);
                    }
                }
            }
        }
    }
    out.sort_unstable();
    out.dedup();
    buffer.human.sort_unstable();
    buffer.human.dedup();
}

/// A sampled k-mer needs this many distinct trinucleotides of its 17. Poly-A, (AT)n,
/// (CTT)n and their kin occur in every genome and in junk reads: on an RNA004 run of human
/// cells, reads of nothing but such repeats matched AT-rich Mycoplasma by chance.
const MIN_TRINUCLEOTIDES: u32 = 8;

/// Human k-mers must be more complex still...
const HUMAN_MIN_TRINUCLEOTIDES: u32 = 10;

/// ...and of middling G+C. The AT-rich 19-mers that fill human repeats turn up by chance in
/// AT-rich genomes: with no bound, 1% of reads from a 27% GC Borrelia isolate were called
/// human; with 5 to 14 G+C of 19, 14 reads still were, two thirds of their hits on k-mers
/// with 5 or 6. The sketch holds 5 to 14 (built so); reads are tested on 7 to 14.
const HUMAN_GC_BASES: std::ops::RangeInclusive<u32> = 7..=14;

/// Distinct trinucleotides in a k-mer (the same on both strands).
fn distinct_trinucleotides(kmer: u64) -> u32 {
    let mut seen = 0u64;
    for i in 0..K - 2 {
        seen |= 1 << ((kmer >> (2 * i)) & 0x3f);
    }
    seen.count_ones()
}

/// G and C bases in a k-mer: a base is G or C when its two bits differ (C = 01, G = 10).
fn gc_bases(kmer: u64) -> u32 {
    ((kmer ^ (kmer >> 1)) & 0x15_5555_5555).count_ones()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn reference_codes(screen_key: &str) -> Vec<u8> {
        let mut text = String::new();
        GzDecoder::new(REFERENCES)
            .read_to_string(&mut text)
            .expect("valid gzip");
        let mut lines = text.lines();
        while let Some(header) = lines.next() {
            let sequence = lines.next().expect("sequence");
            if header[1..].starts_with(screen_key) {
                return sequence
                    .bytes()
                    .map(|b| ASCII_TO_CODE[usize::from(b)])
                    .collect();
            }
        }
        panic!("no reference {screen_key}");
    }

    /// Substitutes about one base in `every`, at pseudo-random positions (fixed seed).
    fn damage(codes: &[u8], every: u64) -> Vec<u8> {
        let mut state = 0x9e37_79b9_7f4a_7c15u64;
        codes
            .iter()
            .map(|&c| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                if state.is_multiple_of(every) {
                    (c + 1) % 4
                } else {
                    c
                }
            })
            .collect()
    }

    #[test]
    fn assigns_noisy_reads_to_their_reference() {
        let screen = Screen::load().expect("built-in references");
        let keys: Vec<&str> = screen.references().iter().map(|r| r.key.as_str()).collect();
        assert_eq!(
            keys,
            [
                "lambda",
                "eno2",
                "phix",
                "ecoli",
                "chrM",
                "ercc",
                "sirv",
                "m_hyorhinis",
                "m_orale",
                "m_fermentans",
                "m_arginini",
                "a_laidlawii",
                "univec",
                "human"
            ]
        );
        let mut buffer = ScreenBuffer::default();
        let ecoli = reference_codes("ecoli");
        // An E. coli read with 5% substitutions (about Q13), reverse-complemented.
        let mut read = damage(&ecoli[1_000_000..1_010_000], 20);
        read.reverse();
        for c in &mut read {
            *c = 3 - *c;
        }
        assert_eq!(screen.classify(&read, None, &mut buffer), Some(3));
        let lambda = reference_codes("lambda");
        assert_eq!(
            screen.classify(&damage(&lambda[100..3700], 25), None, &mut buffer),
            Some(0)
        );
    }

    #[test]
    fn raises_the_bar_for_accurate_reads() {
        let screen = Screen::load().expect("built-in references");
        let mut buffer = ScreenBuffer::default();
        let ecoli = reference_codes("ecoli");
        // A relative: 20% of a read is E. coli sequence, the rest is not.
        let mut state = 11u64;
        let mut relative: Vec<u8> = (0..8_000)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                (state % 4) as u8
            })
            .collect();
        relative.extend_from_slice(&ecoli[2_000_000..2_002_000]);
        assert_eq!(
            screen.classify(&relative, None, &mut buffer),
            Some(3),
            "8% bar without Q"
        );
        assert_eq!(
            screen.classify(&relative, Some(30.0), &mut buffer),
            None,
            "HiFi bar"
        );
        // A real E. coli HiFi read (0.1% errors) still matches at Q30.
        let hifi = damage(&ecoli[3_000_000..3_010_000], 1000);
        assert_eq!(screen.classify(&hifi, Some(30.0), &mut buffer), Some(3));
        // A Q10 read is held only to the 8% bar.
        let noisy = damage(&ecoli[3_000_000..3_010_000], 10);
        assert_eq!(screen.classify(&noisy, Some(10.0), &mut buffer), Some(3));
    }

    #[test]
    fn assigns_spike_in_and_mycoplasma_reads() {
        let screen = Screen::load().expect("built-in references");
        let mut buffer = ScreenBuffer::default();
        for (key, start, len) in [
            ("ercc", 0, 900),
            ("sirv", 0, 3000),
            ("m_hyorhinis", 200_000, 8000),
            ("m_orale", 200_000, 8000),
            ("m_fermentans", 200_000, 8000),
            ("m_arginini", 200_000, 8000),
            ("a_laidlawii", 200_000, 8000),
        ] {
            let codes = reference_codes(key);
            let id = screen
                .references()
                .iter()
                .position(|r| r.key == key)
                .expect("reference");
            // About Q14: 4% substitutions.
            let read = damage(&codes[start..start + len], 25);
            assert_eq!(screen.classify(&read, None, &mut buffer), Some(id), "{key}");
        }
    }

    /// Pseudo-random bases at 40% G+C (fixed seed).
    fn genome(length: usize, seed: u64) -> Vec<u8> {
        let mut state = seed;
        (0..length)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                match state % 10 {
                    0..=2 => 0,
                    3..=5 => 3,
                    6 | 7 => 1,
                    _ => 2,
                }
            })
            .collect()
    }

    /// The built-in screen with its human sketch replaced by one of `human`.
    fn screen_with_human(human: &[u8]) -> Screen {
        let mut screen = Screen::load().expect("built-in references");
        let mut buffer = ScreenBuffer::default();
        sample(human, &mut buffer);
        let mut keys: Vec<u32> = buffer.human.iter().map(|&m| human_key(m)).collect();
        keys.sort_unstable();
        keys.dedup();
        screen.human = keys;
        screen
    }

    #[test]
    fn calls_long_human_reads_and_leaves_short_ones_untested() {
        let human = genome(3_000_000, 0x1234_5678_9abc_def1);
        let screen = screen_with_human(&human);
        let mut buffer = ScreenBuffer::default();
        // A 20 kb read at about Q20 (1% substitutions).
        let read = damage(&human[1_000_000..1_020_000], 100);
        assert_eq!(
            screen.classify(&read, Some(20.0), &mut buffer),
            Some(screen.human_id)
        );
        assert!(Screen::sparse_testable(&buffer));
        // A 1 kb read samples about one k-mer at human density: too few to test.
        let short = damage(&human[2_000_000..2_001_000], 100);
        assert_eq!(screen.classify(&short, Some(20.0), &mut buffer), None);
        assert!(!Screen::sparse_testable(&buffer));
        // An unrelated genome of the same composition is not human.
        let other = genome(20_000, 0x0fed_cba9_8765_4321);
        assert_eq!(screen.classify(&other, Some(20.0), &mut buffer), None);
    }

    #[test]
    fn one_chance_human_hit_does_not_outrank_lambda() {
        let lambda = reference_codes("lambda");
        let mut buffer = ScreenBuffer::default();
        // A noisy 1 kb lambda read that samples one k-mer at human density, and a "human"
        // genome that happens to hold that k-mer: one human hit of one is a 100% share,
        // above lambda's, but too few hits to count.
        let read = (0..lambda.len() - 1000)
            .step_by(500)
            .map(|s| damage(&lambda[s..s + 1000], 25))
            .find(|read| {
                sample(read, &mut buffer);
                buffer.human.len() == 1
            })
            .expect("a lambda read with one human-density k-mer");
        let mut human = genome(500_000, 99);
        human.push(BASE_N);
        human.extend_from_slice(&read);
        let screen = screen_with_human(&human);
        assert_eq!(screen.classify(&read, None, &mut buffer), Some(0));
    }

    #[test]
    fn recognizes_mitochondrial_headers() {
        assert!(is_mitochondrial(
            "NC_012920.1 Homo sapiens mitochondrion, complete genome"
        ));
        assert!(is_mitochondrial("chrM"));
        assert!(is_mitochondrial(
            "MT dna:chromosome chromosome:GRCh38:MT:1:16569:1 REF"
        ));
        assert!(!is_mitochondrial(
            "NC_060925.1 Homo sapiens isolate CHM13 chromosome 1"
        ));
        assert!(!is_mitochondrial("chr1"));
    }

    #[test]
    fn rejects_overlong_sketch_varints() {
        let mut bytes = HUMAN_MAGIC.to_vec();
        bytes.extend_from_slice(&0u64.to_le_bytes());
        bytes.extend_from_slice(&1u32.to_le_bytes());
        bytes.extend_from_slice(&[0xff, 0xff, 0xff, 0xff, 0x7f]);
        let mut gz = GzEncoder::new(Vec::new(), flate2::Compression::fast());
        gz.write_all(&bytes).expect("in memory");
        assert!(read_human_sketch(&gz.finish().expect("in memory")).is_err());
    }

    #[test]
    fn the_human_sketch_holds_about_two_million_kmers() {
        let screen = Screen::load().expect("built-in references");
        assert!((1_500_000..3_000_000).contains(&screen.human.len()));
        assert!(screen.human.is_sorted());
        let human = &screen.references()[screen.human_id];
        assert!(human.length > 3_000_000_000);
    }

    #[test]
    fn samples_about_one_kmer_in_24() {
        let mut state = 7u64;
        let bases: Vec<u8> = (0..240_000)
            .map(|_| {
                state = state
                    .wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1);
                (state >> 62) as u8
            })
            .collect();
        let mut buffer = ScreenBuffer::default();
        sample(&bases, &mut buffer);
        let n = buffer.sampled.len();
        assert!(
            (9_000..11_000).contains(&n),
            "{n} sampled from 240,000 k-mers"
        );
    }

    #[test]
    fn assigns_reads_from_lambda_regions_shared_with_e_coli() {
        let screen = Screen::load().expect("built-in references");
        let mut buffer = ScreenBuffer::default();
        let lambda = reference_codes("lambda");
        // Every 1 kb window of lambda, with 4% substitutions, goes to lambda, including
        // those overlapping the prophage-like stretches E. coli K-12 shares.
        let missed: Vec<usize> = (0..lambda.len() - 1000)
            .step_by(1000)
            .filter(|&start| {
                screen.classify(&damage(&lambda[start..start + 1000], 25), None, &mut buffer)
                    != Some(0)
            })
            .collect();
        assert!(
            missed.is_empty(),
            "lambda windows not assigned at {missed:?}"
        );
    }

    #[test]
    fn leaves_unrelated_reads_alone() {
        let screen = Screen::load().expect("built-in references");
        let mut buffer = ScreenBuffer::default();
        // A pseudo-random sequence stands in for human DNA.
        let mut state = 0x2545_f491_4f6c_dd1du64;
        let read: Vec<u8> = (0..20_000)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                (state % 4) as u8
            })
            .collect();
        assert_eq!(screen.classify(&read, None, &mut buffer), None);
    }
}
