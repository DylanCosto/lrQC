//! Barcode sequences, for finding barcodes left in reads: untrimmed barcodes of the read's
//! own sample, and barcodes of other samples (cross-talk, or chimeras across samples).

use crate::adapters::{self, Adapters};
use crate::read::Platform;

/// Every ONT barcode in its kit's flanks, from dorado's `barcode_kits.cpp`: label, dorado
/// name, kits, front flank, barcode, rear flank.
///
/// Each barcode is matched inside its flanks, as dorado matches it. A bare 24-base barcode
/// is too short: a 15-base chance match among hundreds of sequences turns up in a few
/// percent of long reads, and one sequencing error in its middle breaks every 12-mer but
/// one. Flanks are cut to the 12 bases next to the barcode, so that a kit's flanks around
/// another barcode stay far below the share an occurrence needs.
const ONT_BARCODES: &str = include_str!("../assets/ont_barcodes.tsv");

/// PacBio's SMRTbell adapter indexes (16 of bc1001-bc1022, and bc2001-bc2384).
const PACBIO_BARCODES: &str = include_str!("../assets/pacbio_barcodes.fa");

/// Share of a flanked 24-base ONT barcode's k-mers an occurrence needs: about 22 matching
/// bases, which survives one or two sequencing errors and almost never happens by chance.
const ONT_MIN_SHARE: f64 = 0.3;

/// ONT barcodes shorter than this (Twist's 10-base UDIs) are matched exactly: at 0.3 of a
/// 30-base target, chance matches among 800 of them would turn up in a few percent of reads.
const ONT_SHORT_BARCODE: usize = 20;

/// PacBio barcodes (16 bases) are matched exactly: HiFi reads are accurate, and anything
/// looser matches human sequence by chance.
const PACBIO_MIN_SHARE: f64 = 1.0;

/// The flanked ONT barcode targets: label and sequence to match, and the share it needs.
fn ont_targets() -> Vec<(String, Vec<u8>, f64)> {
    ONT_BARCODES
        .lines()
        .filter(|line| !line.starts_with('#'))
        .filter_map(|line| {
            let fields: Vec<&str> = line.split('\t').collect();
            let &[label, _, _, front, barcode, rear] = fields.as_slice() else {
                return None;
            };
            let share = if barcode.len() < ONT_SHORT_BARCODE {
                1.0
            } else {
                ONT_MIN_SHARE
            };
            Some((
                label.to_owned(),
                format!("{front}{barcode}{rear}").into_bytes(),
                share,
            ))
        })
        .collect()
}

/// PacBio's barcodes named in `names` (such as `bc2001`), with their built-in sequences.
///
/// Only the run's own barcodes are searched for: chance 16-base matches to the hundreds of
/// unused ones would read as cross-talk.
pub fn pacbio_built_in<'a>(names: impl IntoIterator<Item = &'a str>) -> Vec<(String, Vec<u8>)> {
    let wanted: std::collections::HashSet<&str> = names.into_iter().collect();
    adapters::parse_fasta(PACBIO_BARCODES)
        .into_iter()
        .filter(|(name, _)| wanted.contains(name.as_str()))
        .collect()
}

/// Barcode indexes by platform.
#[derive(Debug)]
pub struct Barcodes {
    ont: Adapters,
    pacbio: Option<Adapters>,
}

impl Barcodes {
    /// ONT's barcodes, and PacBio's from `pacbio` (lima's FASTA, or [`pacbio_built_in`]).
    pub fn new(pacbio: &[(String, Vec<u8>)]) -> Self {
        let pacbio = Adapters::from_sequences(pacbio.to_vec(), PACBIO_MIN_SHARE);
        Self {
            ont: Adapters::with_shares(ont_targets()),
            pacbio: (!pacbio.names().is_empty()).then_some(pacbio),
        }
    }

    /// The index to scan reads of `platform` with.
    pub fn for_platform(&self, platform: Platform) -> Option<&Adapters> {
        match platform {
            Platform::Ont => Some(&self.ont),
            Platform::PacBio => self.pacbio.as_ref(),
            Platform::Unknown => None,
        }
    }
}

/// Whether barcode `found` belongs to group `label`: ONT `barcode07` is its own group,
/// PacBio `bc2001` belongs to `bc2001--bc2001` and to `bc2001--bc2002`.
pub fn is_own(label: &str, found: &str) -> bool {
    label == found || label.split("--").any(|part| part == found)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::adapters::ScanBuffer;
    use crate::read::ASCII_TO_CODE;

    fn scan(index: &Adapters, seq: &str) -> Vec<String> {
        let read: Vec<u8> = seq.bytes().map(|b| ASCII_TO_CODE[usize::from(b)]).collect();
        let mut hits = Vec::new();
        index.scan(&read, &mut ScanBuffer::default(), &mut hits);
        hits.iter()
            .map(|h| index.names()[h.adapter].clone())
            .collect()
    }

    #[test]
    fn finds_flanked_ont_barcodes_on_either_strand() {
        let barcodes = Barcodes::new(&[]);
        let index = barcodes.for_platform(Platform::Ont).expect("ONT index");
        // Every flanked target in dorado's kits, with Twist's 800 UDIs.
        assert_eq!(index.names().len(), 1371);
        let genomic = "GATTACAGGCTTAGCCATGCAAGTCCGATTGCAAGTCATGCATTGACCATTGACCAGT";
        // Native barcode 5 with its flanks, and the same reverse-complemented.
        let native = "ATTGCTAAGGTTAAAAGGTTACACAAACCCTGGACAAGCAGCACCT";
        let reverse: String = native
            .chars()
            .rev()
            .map(|c| match c {
                'A' => 'T',
                'C' => 'G',
                'G' => 'C',
                _ => 'A',
            })
            .collect();
        for target in [native, reverse.as_str()] {
            assert_eq!(
                scan(index, &format!("{genomic}{target}{genomic}")),
                ["barcode05"]
            );
        }
        // One error in the middle of the barcode still matches.
        let damaged = native.replacen("CAAACCC", "CAATCCC", 1);
        assert_eq!(
            scan(index, &format!("{genomic}{damaged}{genomic}")),
            ["barcode05"]
        );
        // A bare barcode with a different flank's worth of genome around it does not.
        assert!(scan(index, &format!("{genomic}AAGGTTAC{genomic}")).is_empty());
    }

    #[test]
    fn a_kits_flanks_around_another_barcode_never_match() {
        let barcodes = Barcodes::new(&[]);
        let index = barcodes.for_platform(Platform::Ont).expect("ONT index");
        let genomic = "GATTACAGGCTTAGCCATGCAAGTCCGATTGCAAGTCATGCATTGACCATTGACCAGT";
        for line in ONT_BARCODES.lines().filter(|l| !l.starts_with('#')) {
            let fields: Vec<&str> = line.split('\t').collect();
            let (front, barcode, rear) = (fields[3], fields[4], fields[5]);
            // The flanks around a barcode-length run of a sequence no barcode has.
            let filler = "T".repeat(barcode.len());
            let read = format!("{genomic}{front}{filler}{rear}{genomic}");
            assert!(scan(index, &read).is_empty(), "{line}");
        }
    }

    #[test]
    fn finds_every_kits_barcodes() {
        let barcodes = Barcodes::new(&[]);
        let index = barcodes.for_platform(Platform::Ont).expect("ONT index");
        let genomic = "GATTACAGGCTTAGCCATGCAAGTCCGATTGCAAGTCATGCATTGACCATTGACCAGT";
        for line in ONT_BARCODES.lines().filter(|l| !l.starts_with('#')) {
            let fields: Vec<&str> = line.split('\t').collect();
            let target = format!("{}{}{}", fields[3], fields[4], fields[5]);
            let found = scan(index, &format!("{genomic}{target}{genomic}"));
            assert!(
                found.iter().any(|name| name == fields[0]),
                "{line}: {found:?}"
            );
        }
    }

    #[test]
    fn built_in_pacbio_barcodes_are_the_runs_own() {
        let list = pacbio_built_in(["bc1016", "bc2001", "bc9999"]);
        let mut names: Vec<&str> = list.iter().map(|(name, _)| name.as_str()).collect();
        names.sort_unstable();
        assert_eq!(names, ["bc1016", "bc2001"]);
        assert_eq!(
            list[0].1, b"ATCGTGCGACGAGTAT",
            "bc2001 comes first in PacBio's file"
        );
        assert_eq!(adapters::parse_fasta(PACBIO_BARCODES).len(), 400);
    }

    #[test]
    fn matches_pacbio_pairs() {
        assert!(is_own("barcode05", "barcode05"));
        assert!(!is_own("barcode05", "barcode06"));
        assert!(is_own("bc2001--bc2002", "bc2002"));
        assert!(!is_own("bc2001--bc2001", "bc2003"));
    }
}
