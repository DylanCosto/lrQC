//! The per-read view that the parsers fill and the statistics consume.

use std::mem;
use std::sync::Arc;

use serde::{Deserialize, Serialize};

/// Base codes used throughout: A, C, G, T and anything else.
pub const BASE_A: u8 = 0;
pub const BASE_C: u8 = 1;
pub const BASE_G: u8 = 2;
pub const BASE_T: u8 = 3;
pub const BASE_N: u8 = 4;

/// Maps an ASCII base to its code; IUPAC ambiguity codes become [`BASE_N`].
pub const ASCII_TO_CODE: [u8; 256] = {
    let mut table = [BASE_N; 256];
    table[b'A' as usize] = BASE_A;
    table[b'a' as usize] = BASE_A;
    table[b'C' as usize] = BASE_C;
    table[b'c' as usize] = BASE_C;
    table[b'G' as usize] = BASE_G;
    table[b'g' as usize] = BASE_G;
    table[b'T' as usize] = BASE_T;
    table[b't' as usize] = BASE_T;
    table[b'U' as usize] = BASE_T;
    table[b'u' as usize] = BASE_T;
    table
};

/// Sequencing platform a read came from.
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Default, Serialize, Deserialize,
)]
pub enum Platform {
    Ont,
    PacBio,
    #[default]
    Unknown,
}

impl Platform {
    /// Human-readable platform name.
    pub fn label(self) -> &'static str {
        match self {
            Platform::Ont => "Oxford Nanopore",
            Platform::PacBio => "PacBio",
            Platform::Unknown => "unknown",
        }
    }
}

/// The basecaller's or instrument's pass/fail verdict on a read.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Filter {
    Pass,
    Fail,
}

/// How a BAM record relates to the read it holds.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum RecordKind {
    #[default]
    Primary,
    Secondary,
    Supplementary,
}

/// Alignment state of a primary record.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub enum Alignment {
    /// The input carries no alignments (FASTQ, unaligned BAM).
    #[default]
    NotAligned,
    /// An aligned input in which this read did not map.
    Unmapped,
    Mapped(Mapped),
}

/// What the CIGAR and tags of a mapped primary record say.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Mapped {
    pub mapq: u8,
    /// `NM`: mismatches plus inserted and deleted bases.
    pub edit_distance: Option<u32>,
    /// Bases in M, = and X operations.
    pub aligned_columns: u64,
    pub inserted: u64,
    pub deleted: u64,
    pub insertion_runs: u64,
    pub deletion_runs: u64,
    /// The record has an `SA` tag: the read has supplementary alignments.
    pub split: bool,
}

/// One read, with buffers reused from record to record.
#[derive(Debug, Default)]
pub struct Read {
    /// The read name, as in the input.
    pub name: Vec<u8>,
    pub name_hash: u64,
    pub kind: RecordKind,
    /// Length in bases; set even when the input has no sequence (sequencing summaries).
    pub len: u64,
    /// Base codes in sequencing orientation; empty when the input has no sequence.
    pub bases: Vec<u8>,
    /// Phred qualities in sequencing orientation; empty when the input has none.
    pub quals: Vec<u8>,
    pub barcode: Option<Arc<str>>,
    pub platform: Platform,
    pub filter: Option<Filter>,
    /// Mean read Q the basecaller reported (ONT `qs`, summary `mean_qscore_template`).
    pub reported_q: Option<f32>,
    /// PacBio predicted accuracy (`rq`).
    pub accuracy: Option<f32>,
    /// PacBio CCS passes (`np`).
    pub passes: Option<u32>,
    /// PacBio effective coverage (`ec`).
    pub effective_coverage: Option<f32>,
    /// PacBio lima barcode score (`bq`).
    pub barcode_score: Option<u32>,
    pub kinetics: Kinetics,
    pub channel: Option<u32>,
    /// Read start in seconds: Unix time for BAM and FASTQ, run time for summaries.
    pub start_secs: Option<f64>,
    pub duration_secs: Option<f32>,
    pub end_reason: Option<Arc<str>>,
    /// Raw `MM` tag text; empty when the read has none.
    pub mod_positions: Vec<u8>,
    /// Raw `ML` probabilities.
    pub mod_probs: Vec<u8>,
    pub alignment: Alignment,
    /// dorado duplex status (`dx`).
    pub duplex: Option<DuplexStatus>,
    /// dorado split this read from a longer one (`pi` names the parent).
    pub origin: Origin,
    /// ONT mux, the well of the channel (`mx`).
    pub mux: Option<u8>,
    /// ONT signal shift and scale (`sm`, `sd`) and raw samples (`ns`).
    pub signal_shift: Option<f32>,
    pub signal_scale: Option<f32>,
    pub signal_samples: Option<u64>,
    /// dorado's estimated poly(A) tail length (`pt`).
    pub poly_a: Option<u32>,
    /// PacBio ZMW hole number (`zm`, or the read name's middle field).
    pub zmw: Option<u64>,
    /// Sequencing run: flow cell, run ID or PacBio movie.
    pub run: Option<Arc<str>>,
    pub molecule: Molecule,
}

impl Read {
    /// Drops tag values no basecaller writes (NaN, infinities, negative Qs, an `rq`
    /// outside (0, 1]), so a damaged or hand-edited record falls back to what its bases
    /// say instead of skewing every summary.
    pub fn sanitize(&mut self) {
        fn keep(value: &mut Option<f32>, ok: impl Fn(f32) -> bool) {
            if value.is_some_and(|v| !(v.is_finite() && ok(v))) {
                *value = None;
            }
        }
        keep(&mut self.reported_q, |q| (0.0..=100.0).contains(&q));
        keep(&mut self.accuracy, |rq| rq > 0.0 && rq <= 1.0);
        keep(&mut self.effective_coverage, |ec| ec >= 0.0);
        keep(&mut self.duration_secs, |du| du >= 0.0);
        keep(&mut self.signal_shift, |_| true);
        keep(&mut self.signal_scale, |_| true);
        if self.start_secs.is_some_and(|st| !st.is_finite()) {
            self.start_secs = None;
        }
    }
}

/// Whether a `BC` tag holds a barcode's sequence (as the SAM specification has it)
/// rather than a demultiplexer's label (`barcode01`, as dorado writes it).
pub fn is_barcode_sequence(value: &[u8]) -> bool {
    value.len() >= 4
        && value.iter().all(|b| {
            matches!(
                b.to_ascii_uppercase(),
                b'A' | b'C' | b'G' | b'T' | b'N' | b'-' | b'+'
            )
        })
}

/// What was sequenced.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
pub enum Molecule {
    #[default]
    Dna,
    /// ONT direct RNA: single-stranded transcripts, read 3' to 5'.
    Rna,
}

impl Molecule {
    /// Direct RNA when a dorado model or read group name says so (`rna004_130bps_sup`,
    /// `..._rna004_sup@v6.0.0`).
    pub fn from_model_name(name: &str) -> Self {
        let name = name.to_ascii_lowercase();
        if name.contains("rna00") || name.contains("_rna_") {
            Molecule::Rna
        } else {
            Molecule::Dna
        }
    }
}

/// dorado's duplex classification.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DuplexStatus {
    /// A simplex read with no duplex offspring (`dx:i:0`).
    Simplex,
    /// A duplex read (`dx:i:1`).
    Duplex,
    /// A simplex read whose strand pair produced a duplex read (`dx:i:-1`).
    SimplexWithDuplexOffspring,
}

impl DuplexStatus {
    /// The duplex status a `dx` value encodes.
    pub fn from_tag(value: i64) -> Option<Self> {
        match value {
            0 => Some(DuplexStatus::Simplex),
            1 => Some(DuplexStatus::Duplex),
            -1 => Some(DuplexStatus::SimplexWithDuplexOffspring),
            _ => None,
        }
    }
}

/// Whether a read is whole or a piece of a read the basecaller split.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Origin {
    #[default]
    Whole,
    SplitFromParent,
}

/// PacBio ZMW from a read name such as `m84011_220902_175841_s1/182061619/ccs`.
pub fn zmw_from_name(name: &[u8]) -> Option<u64> {
    let mut fields = name.split(|&b| b == b'/');
    let movie = fields.next()?;
    let zmw = fields.next()?;
    if !movie.starts_with(b"m") {
        return None;
    }
    std::str::from_utf8(zmw).ok()?.parse().ok()
}

/// PacBio movie from a read name, the run of the read.
pub fn movie_from_name(name: &[u8]) -> Option<&[u8]> {
    let movie = name.split(|&b| b == b'/').next()?;
    (movie.starts_with(b"m") && name.contains(&b'/')).then_some(movie)
}

/// Whether a PacBio read carries polymerase kinetics tags.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Kinetics {
    #[default]
    Absent,
    Present,
}

impl Read {
    /// Resets every field, keeping the buffers' allocations.
    pub fn clear(&mut self) {
        let mut bases = mem::take(&mut self.bases);
        let mut quals = mem::take(&mut self.quals);
        let mut mod_positions = mem::take(&mut self.mod_positions);
        let mut mod_probs = mem::take(&mut self.mod_probs);
        let mut name = mem::take(&mut self.name);
        name.clear();
        bases.clear();
        quals.clear();
        mod_positions.clear();
        mod_probs.clear();
        *self = Read {
            name,
            bases,
            quals,
            mod_positions,
            mod_probs,
            ..Read::default()
        };
    }

    /// Reverse-complements the bases and reverses the qualities, restoring sequencing
    /// orientation for reads aligned to the reverse strand.
    pub fn reverse_complement(&mut self) {
        self.bases.reverse();
        for base in &mut self.bases {
            if *base != BASE_N {
                *base = BASE_T - *base;
            }
        }
        self.quals.reverse();
    }
}

/// A 64-bit FNV-1a hash: stable across runs and platforms, so sampling is reproducible.
pub fn hash_name(name: &[u8]) -> u64 {
    const OFFSET: u64 = 0xcbf2_9ce4_8422_2325;
    const PRIME: u64 = 0x0000_0100_0000_01b3;
    let mut hash = OFFSET;
    for &byte in name {
        hash ^= u64::from(byte);
        hash = hash.wrapping_mul(PRIME);
    }
    // Finalize (splitmix64) so the low and high bits are both well mixed.
    hash ^= hash >> 30;
    hash = hash.wrapping_mul(0xbf58_476d_1ce4_e5b9);
    hash ^= hash >> 27;
    hash = hash.wrapping_mul(0x94d0_49bb_1331_11eb);
    hash ^ (hash >> 31)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reverse_complement_restores_orientation() {
        let mut read = Read {
            bases: vec![BASE_A, BASE_C, BASE_N, BASE_G],
            quals: vec![1, 2, 3, 4],
            ..Read::default()
        };
        read.reverse_complement();
        assert_eq!(read.bases, vec![BASE_C, BASE_N, BASE_G, BASE_T]);
        assert_eq!(read.quals, vec![4, 3, 2, 1]);
    }

    #[test]
    fn recognizes_direct_rna_models() {
        assert_eq!(
            Molecule::from_model_name("30374973-806b_rna004_sup@v6.0.0"),
            Molecule::Rna
        );
        assert_eq!(
            Molecule::from_model_name("basecall_model=rna004_130bps_sup@v5.1.0"),
            Molecule::Rna
        );
        assert_eq!(
            Molecule::from_model_name("70eb6e24_dna_r10.4.1_e8.2_400bps_sup@v5.0.0"),
            Molecule::Dna
        );
    }

    #[test]
    fn clear_keeps_capacity() {
        let mut read = Read {
            bases: Vec::with_capacity(100),
            len: 5,
            ..Read::default()
        };
        read.clear();
        assert_eq!(read.len, 0);
        assert!(read.bases.capacity() >= 100);
    }
}
