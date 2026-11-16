//! Per-group accumulators: everything the report needs, gathered in one pass over the reads.
//!
//! Every field is a count or a sum, so accumulators from different threads merge exactly.

use std::collections::{HashMap, HashSet};
use std::sync::{Arc, LazyLock};

use crate::adapters::{Adapters, Hit, Role, ScanBuffer};
use crate::barcodes::{self, Barcodes};
use crate::bins::{self, POSITIONS};
use crate::read::{
    Alignment, DuplexStatus, Filter, Kinetics, Molecule, Origin, Platform, Read, RecordKind,
    BASE_A, BASE_C, BASE_G, BASE_T,
};
use crate::sampler::EndSample;
use crate::screen::{Screen, ScreenBuffer};

/// Phred values 0 to 93, the range SAM and FASTQ can store.
pub const QUAL_LEVELS: usize = 94;

/// Read Q in tenths, 0.0 to 60.0.
pub const READ_Q_BINS: usize = 601;

/// Whole read Q values in the length/quality map, 0 to 60.
pub const READ_Q_LEVELS: usize = 61;

/// Translocation speed bins of 10 bases per second, up to 2,000.
pub const SPEED_BINS: usize = 200;
pub const SPEED_BIN_WIDTH: f64 = 10.0;

/// Identity bins of 0.1 percent.
pub const IDENTITY_BINS: usize = 1001;

/// Bases from a read's end searched for poly(A) and poly(T) tails.
const POLY_WINDOW: usize = 100;

/// Shortest homopolymer counted as a poly(A) or poly(T) tail.
const POLY_MIN: usize = 15;

/// Adapters starting within this many bases of a read end count as at that end.
pub const END_WINDOW: usize = 150;

/// ML probabilities at or above this (0.5) count a call as modified.
const MOD_THRESHOLD: u8 = 128;

/// Error probability of each Phred value.
static ERROR_PROB: LazyLock<[f64; QUAL_LEVELS]> = LazyLock::new(|| {
    let mut table = [0.0; QUAL_LEVELS];
    for (q, p) in table.iter_mut().enumerate() {
        *p = 10f64.powf(-(q as f64) / 10.0);
    }
    table
});

/// Accumulated statistics of one group of reads.
#[derive(Debug, Clone)]
pub struct Stats {
    pub platform_reads: [u64; 3],
    pub reads: u64,
    pub bases: u64,
    pub lengths: Vec<u32>,
    pub read_q: Vec<u64>,
    pub read_q_bases: Vec<u64>,
    /// Reads whose exact read Q reaches Q10, Q20 and Q30, and the bases of those reaching
    /// Q20: the histograms round to tenths, which would count a Q19.96 read as Q20.
    pub reads_at_q: [u64; 3],
    pub bases_at_q20: u64,
    pub reads_with_q: u64,
    /// Reads by log length bin (outer) and whole read Q (inner).
    pub length_q: Vec<u64>,
    pub reads_with_seq: u64,
    pub reads_with_quals: u64,
    /// Quality counts by position bin (outer) and Phred value (inner), from each end.
    pub quals_from_start: Vec<u64>,
    pub quals_from_end: Vec<u64>,
    /// Base counts by position bin (outer) and base code (inner), from each end.
    pub bases_from_start: Vec<u64>,
    pub bases_from_end: Vec<u64>,
    pub base_quals: [u64; QUAL_LEVELS],
    pub gc: [u64; 101],
    /// Reads by GC percent, only reads of [`GC_PEAK_MIN_LENGTH`] bases or more.
    pub gc_long: [u64; 101],
    pub pass: Yield,
    pub fail: Yield,
    pub adapters: Vec<AdapterCounts>,
    /// Reads with any adapter or primer within the end window of either end.
    pub reads_adapter_at_end: u64,
    /// Reads with any adapter or primer inside, away from both ends.
    pub reads_adapter_internal: u64,
    pub poly_a_tail: u64,
    pub poly_t_head: u64,
    pub sample: EndSample,
    pub name_hashes: HashSet<u64>,
    pub duplicate_names: u64,
    /// PacBio reads whose read Q comes from their base qualities, the input lacking `rq`.
    pub pacbio_q_estimated: u64,
    pub time: HashMap<i64, TimeBin>,
    pub channels: HashMap<u32, Yield>,
    pub speed: Vec<u64>,
    pub end_reasons: HashMap<Arc<str>, u64>,
    pub passes: [u64; 101],
    pub effective_coverage: [u64; 101],
    pub barcode_score: [u64; 101],
    pub with_kinetics: u64,
    pub reads_with_mods: u64,
    pub mods: HashMap<String, ModCounts>,
    pub mapping: MappingCounts,
    /// Read Q in tenths, per read, in the order of `lengths`; [`NO_Q`] when unknown.
    pub read_q_tenths: Vec<u16>,
    pub ont: OntCounts,
    pub pacbio: PacBioCounts,
    pub cdna: CdnaCounts,
    pub barcodes: BarcodeCounts,
    /// Reads and bases assigned to each screen reference.
    pub screen: Vec<Yield>,
    pub reads_screened: u64,
    /// Reads long enough to test against sparse references, and of them, the reads assigned
    /// to each reference.
    pub sparse_testable: u64,
    pub screen_testable: Vec<u64>,
    /// Reads by Kinnex adapters found, up to 40.
    pub kinnex: Vec<u64>,
    /// Reads of direct RNA.
    pub rna_reads: u64,
}

/// GC peaks are judged on reads at least this long. Shorter reads in long-read runs are
/// largely adapter and barcode fragments (about 45% GC), which form a GC peak of their own;
/// the adapter and overrepresented-sequence checks grade those.
pub const GC_PEAK_MIN_LENGTH: u64 = 500;

/// Marks a read without a read Q in [`Stats::read_q_tenths`].
pub const NO_Q: u16 = u16::MAX;

/// Bins of the ONT signal and poly(A) distributions.
pub const SHIFT_BIN: f64 = 2.0;
pub const SCALE_BIN: f64 = 1.0;
pub const SAMPLES_PER_BASE_BIN: f64 = 0.25;
pub const POLY_A_BIN: u32 = 5;
const POLY_A_BINS: usize = 121;
const KINNEX_MAX: usize = 40;

/// ONT-specific counts.
#[derive(Debug, Clone)]
pub struct OntCounts {
    /// Simplex, duplex, and simplex reads with duplex offspring.
    pub duplex: [Yield; 3],
    /// Read Q (tenths) by duplex class: `READ_Q_BINS` per class.
    pub duplex_q: Vec<u64>,
    pub split_reads: Yield,
    pub mux: HashMap<u8, Yield>,
    /// Signal shift, scale and samples per base, by bin.
    pub shift: HashMap<i64, u64>,
    pub scale: HashMap<i64, u64>,
    pub samples_per_base: HashMap<i64, u64>,
    /// Reads by poly(A) length bin of 5 bases; the last bin holds 600 and more.
    pub poly_a: Vec<u64>,
}

impl Default for OntCounts {
    fn default() -> Self {
        Self {
            duplex: [Yield::default(); 3],
            duplex_q: vec![0; 3 * READ_Q_BINS],
            split_reads: Yield::default(),
            mux: HashMap::new(),
            shift: HashMap::new(),
            scale: HashMap::new(),
            samples_per_base: HashMap::new(),
            poly_a: vec![0; POLY_A_BINS],
        }
    }
}

/// PacBio-specific counts.
#[derive(Debug, Clone, Default)]
pub struct PacBioCounts {
    /// Reads by passes (outer, up to 60) and whole read Q (inner).
    pub passes_q: Vec<u64>,
    /// Reads per ZMW.
    pub zmw_reads: HashMap<u64, u32>,
}

/// cDNA primer structure of reads.
#[derive(Debug, Clone, Copy, Default)]
pub struct CdnaCounts {
    /// Reads with the strand-switching primer at the start and the VN primer at the end.
    pub full_length_forward: u64,
    /// The same, read from the other strand: VN primer first.
    pub full_length_reverse: u64,
    pub ssp_only: u64,
    pub vnp_only: u64,
}

/// Barcodes found inside reads.
#[derive(Debug, Clone, Default)]
pub struct BarcodeCounts {
    pub reads_scanned: u64,
    /// The read's own barcode within 150 bases of an end: left untrimmed.
    pub own_at_end: u64,
    /// The read's own barcode further inside: a chimera within the sample.
    pub own_internal: u64,
    /// Reads holding another sample's barcode.
    pub with_foreign: u64,
    /// Reads holding each foreign barcode.
    pub foreign: HashMap<String, u64>,
}

/// Reads and bases.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Yield {
    pub reads: u64,
    pub bases: u64,
}

impl Yield {
    fn add(&mut self, bases: u64) {
        self.reads += 1;
        self.bases += bases;
    }

    fn merge(&mut self, other: Yield) {
        self.reads += other.reads;
        self.bases += other.bases;
    }
}

/// Where in reads one adapter was found.
#[derive(Debug, Clone)]
pub struct AdapterCounts {
    pub reads_at_start: u64,
    pub reads_at_end: u64,
    pub reads_internal: u64,
    /// Reads by the position bin of the occurrence nearest each end.
    pub nearest_start: Vec<u64>,
    pub nearest_end: Vec<u64>,
}

impl AdapterCounts {
    fn new() -> Self {
        Self {
            reads_at_start: 0,
            reads_at_end: 0,
            reads_internal: 0,
            nearest_start: vec![0; POSITIONS.len()],
            nearest_end: vec![0; POSITIONS.len()],
        }
    }

    fn merge(&mut self, other: &AdapterCounts) {
        self.reads_at_start += other.reads_at_start;
        self.reads_at_end += other.reads_at_end;
        self.reads_internal += other.reads_internal;
        add_into(&mut self.nearest_start, &other.nearest_start);
        add_into(&mut self.nearest_end, &other.nearest_end);
    }
}

/// One minute of a run.
#[derive(Debug, Clone, Copy, Default)]
pub struct TimeBin {
    pub reads: u64,
    pub bases: u64,
    pub pass_reads: u64,
    /// Sum of read Q in tenths: integers, so merges in any order give the same total.
    pub q_tenths: u64,
    pub q_reads: u64,
}

/// Calls of one modification code (`C+m`, `A+a`, ...).
#[derive(Debug, Clone, Default)]
pub struct ModCounts {
    pub calls: u64,
    pub modified: u64,
    /// Calls on a C followed by G in the read (CpG context).
    pub cpg_calls: u64,
    pub cpg_modified: u64,
    /// Calls by ML probability, in 26 bins of about 4%.
    pub probs: [u64; 26],
}

/// Alignment counts, for aligned BAMs.
#[derive(Debug, Clone)]
pub struct MappingCounts {
    pub mapped: u64,
    pub unmapped: u64,
    pub secondary: u64,
    pub supplementary: u64,
    pub split: u64,
    pub mapped_bases: u64,
    pub mapq: [u64; 61],
    pub identity: Vec<u64>,
    pub gap_compressed_identity: Vec<u64>,
    pub aligned_fraction: [u64; 101],
}

impl Default for MappingCounts {
    fn default() -> Self {
        Self {
            mapped: 0,
            unmapped: 0,
            secondary: 0,
            supplementary: 0,
            split: 0,
            mapped_bases: 0,
            mapq: [0; 61],
            identity: vec![0; IDENTITY_BINS],
            gap_compressed_identity: vec![0; IDENTITY_BINS],
            aligned_fraction: [0; 101],
        }
    }
}

/// Settings that decide how reads are counted.
#[derive(Debug, Clone, Copy)]
pub struct Settings<'cfg> {
    /// Read Q at which reads pass when their input does not say; `None` uses the platform's.
    pub min_q: Option<f32>,
    pub adapters: &'cfg Adapters,
    pub barcodes: &'cfg Barcodes,
    /// The contaminant screen; `None` when turned off.
    pub screen: Option<&'cfg Screen>,
    /// Kinnex/MAS segmentation adapters, from `--kinnex-adapters`.
    pub kinnex: Option<&'cfg Adapters>,
}

impl Settings<'_> {
    /// Read Q threshold for a platform: MinKNOW's for super-accuracy calls, HiFi's Q20.
    pub fn min_q_for(&self, platform: Platform) -> f32 {
        self.min_q.unwrap_or(match platform {
            Platform::PacBio => 20.0,
            Platform::Ont | Platform::Unknown => 10.0,
        })
    }
}

/// Buffers reused from read to read.
#[derive(Debug, Default)]
pub struct Scratch {
    adapter_scan: ScanBuffer,
    hits: Vec<Hit>,
    barcode_hits: Vec<Hit>,
    kinnex_hits: Vec<Hit>,
    screen: ScreenBuffer,
    outcome: ReadOutcome,
}

impl Scratch {
    /// What the last [`Stats::add`] measured about its read.
    pub fn outcome(&self) -> &ReadOutcome {
        &self.outcome
    }
}

/// Per-read results of [`Stats::add`], for the per-read table.
#[derive(Debug, Clone, Default)]
pub struct ReadOutcome {
    /// The record was a primary record, so it counted as a read.
    pub counted: bool,
    pub read_q: Option<f32>,
    pub filter: Option<Filter>,
    pub gc_pct: Option<u8>,
    /// Index of the screen reference the read matched.
    pub screen: Option<usize>,
    /// Whether the read was screened at all (sequence present, screen on).
    pub screened: bool,
    /// Indexes of the adapters found, in read order.
    pub adapters: Vec<usize>,
    /// Barcodes found inside the read.
    pub barcodes: Vec<String>,
    pub kinnex_adapters: Option<usize>,
    /// BLAST-like identity of the primary alignment.
    pub identity: Option<f64>,
}

impl ReadOutcome {
    fn reset(&mut self) {
        let mut adapters = std::mem::take(&mut self.adapters);
        let mut barcodes = std::mem::take(&mut self.barcodes);
        adapters.clear();
        barcodes.clear();
        *self = ReadOutcome {
            adapters,
            barcodes,
            ..ReadOutcome::default()
        };
    }
}

impl Stats {
    pub fn new(adapters: usize) -> Self {
        let bins = POSITIONS.len();
        Self {
            platform_reads: [0; 3],
            reads: 0,
            bases: 0,
            lengths: Vec::new(),
            read_q: vec![0; READ_Q_BINS],
            read_q_bases: vec![0; READ_Q_BINS],
            reads_at_q: [0; 3],
            bases_at_q20: 0,
            reads_with_q: 0,
            length_q: vec![0; bins::LENGTH_BINS * READ_Q_LEVELS],
            reads_with_seq: 0,
            reads_with_quals: 0,
            quals_from_start: vec![0; bins * QUAL_LEVELS],
            quals_from_end: vec![0; bins * QUAL_LEVELS],
            bases_from_start: vec![0; bins * 5],
            bases_from_end: vec![0; bins * 5],
            base_quals: [0; QUAL_LEVELS],
            gc: [0; 101],
            gc_long: [0; 101],
            pass: Yield::default(),
            fail: Yield::default(),
            adapters: (0..adapters).map(|_| AdapterCounts::new()).collect(),
            reads_adapter_at_end: 0,
            reads_adapter_internal: 0,
            poly_a_tail: 0,
            poly_t_head: 0,
            sample: EndSample::default(),
            name_hashes: HashSet::new(),
            duplicate_names: 0,
            pacbio_q_estimated: 0,
            time: HashMap::new(),
            channels: HashMap::new(),
            speed: vec![0; SPEED_BINS],
            end_reasons: HashMap::new(),
            passes: [0; 101],
            effective_coverage: [0; 101],
            barcode_score: [0; 101],
            with_kinetics: 0,
            reads_with_mods: 0,
            mods: HashMap::new(),
            mapping: MappingCounts::default(),
            read_q_tenths: Vec::new(),
            ont: OntCounts::default(),
            pacbio: PacBioCounts {
                passes_q: vec![0; 61 * READ_Q_LEVELS],
                zmw_reads: HashMap::new(),
            },
            cdna: CdnaCounts::default(),
            barcodes: BarcodeCounts::default(),
            screen: Vec::new(),
            reads_screened: 0,
            sparse_testable: 0,
            screen_testable: Vec::new(),
            kinnex: vec![0; KINNEX_MAX + 1],
            rna_reads: 0,
        }
    }

    /// Counts one record.
    pub fn add(&mut self, read: &Read, settings: &Settings<'_>, scratch: &mut Scratch) {
        scratch.outcome.reset();
        match read.kind {
            RecordKind::Secondary => {
                self.mapping.secondary += 1;
                return;
            }
            RecordKind::Supplementary => {
                self.mapping.supplementary += 1;
                return;
            }
            RecordKind::Primary => {}
        }
        scratch.outcome.counted = true;
        self.rna_reads += u64::from(read.molecule == Molecule::Rna);
        self.reads += 1;
        self.bases += read.len;
        self.platform_reads[platform_index(read.platform)] += 1;
        self.lengths
            .push(u32::try_from(read.len).unwrap_or(u32::MAX));
        if !self.name_hashes.insert(read.name_hash) {
            self.duplicate_names += 1;
        }

        let read_q = self.add_qualities(read);
        // Rounded down, so the filter grid's "read Q of at least" is exact at every tenth.
        self.read_q_tenths.push(read_q.map_or(NO_Q, |q| {
            (f64::from(q) * 10.0)
                .floor()
                .clamp(0.0, f64::from(NO_Q - 1)) as u16
        }));
        if let Some(q) = read_q {
            let bin = ((f64::from(q) * 10.0).round().max(0.0) as usize).min(READ_Q_BINS - 1);
            self.read_q[bin] += 1;
            self.read_q_bases[bin] += read.len;
            for (count, level) in self.reads_at_q.iter_mut().zip([10.0, 20.0, 30.0]) {
                *count += u64::from(q >= level);
            }
            if q >= 20.0 {
                self.bases_at_q20 += read.len;
            }
            self.reads_with_q += 1;
            let q_level = (q.round().max(0.0) as usize).min(READ_Q_LEVELS - 1);
            self.length_q[bins::length_bin(read.len) * READ_Q_LEVELS + q_level] += 1;
        }

        // A PacBio read without `rq` (SMRT Link's FASTQ export) has a read Q estimated from
        // binned base qualities, which runs several Q off `rq` either way. CCS output holds
        // only HiFi reads unless told otherwise, so its reads pass rather than fail on that
        // estimate; `--min-qscore` still applies when given.
        let q_estimated = read.platform == Platform::PacBio
            && read.reported_q.is_none()
            && read.accuracy.is_none()
            && read_q.is_some();
        self.pacbio_q_estimated += u64::from(q_estimated);
        let filter = read.filter.or_else(|| {
            if q_estimated && settings.min_q.is_none() && is_ccs_name(&read.name) {
                return Some(Filter::Pass);
            }
            read_q.map(|q| {
                if q >= settings.min_q_for(read.platform) {
                    Filter::Pass
                } else {
                    Filter::Fail
                }
            })
        });
        match filter {
            Some(Filter::Pass) => self.pass.add(read.len),
            Some(Filter::Fail) => self.fail.add(read.len),
            None => {}
        }

        scratch.outcome.read_q = read_q;
        scratch.outcome.filter = filter;
        if !read.bases.is_empty() {
            self.add_sequence(read, settings, scratch);
        }
        scratch.outcome.identity = self.add_alignment(read);
        self.add_run_metadata(read, read_q, filter);
        self.add_pacbio(read);
        if !read.mod_positions.is_empty() {
            self.add_mods(&read.mod_positions, &read.mod_probs, &read.bases);
        }
        self.add_ont(read, read_q);
        self.add_pacbio_extras(read, read_q);
    }

    fn add_ont(&mut self, read: &Read, read_q: Option<f32>) {
        let ont = &mut self.ont;
        if let Some(duplex) = read.duplex {
            let class = match duplex {
                DuplexStatus::Simplex => 0,
                DuplexStatus::Duplex => 1,
                DuplexStatus::SimplexWithDuplexOffspring => 2,
            };
            ont.duplex[class].add(read.len);
            if let Some(q) = read_q {
                let bin = ((f64::from(q) * 10.0).round().max(0.0) as usize).min(READ_Q_BINS - 1);
                ont.duplex_q[class * READ_Q_BINS + bin] += 1;
            }
        }
        if read.origin == Origin::SplitFromParent {
            ont.split_reads.add(read.len);
        }
        if let Some(mux) = read.mux {
            ont.mux.entry(mux).or_default().add(read.len);
        }
        let bump = |map: &mut HashMap<i64, u64>, value: f64, width: f64| {
            *map.entry((value / width).floor() as i64).or_default() += 1;
        };
        if let Some(shift) = read.signal_shift {
            bump(&mut ont.shift, f64::from(shift), SHIFT_BIN);
        }
        if let Some(scale) = read.signal_scale {
            bump(&mut ont.scale, f64::from(scale), SCALE_BIN);
        }
        if let (Some(samples), true) = (read.signal_samples, read.len > 0) {
            bump(
                &mut ont.samples_per_base,
                samples as f64 / read.len as f64,
                SAMPLES_PER_BASE_BIN,
            );
        }
        if let Some(length) = read.poly_a {
            ont.poly_a[((length / POLY_A_BIN) as usize).min(POLY_A_BINS - 1)] += 1;
        }
    }

    fn add_pacbio_extras(&mut self, read: &Read, read_q: Option<f32>) {
        if let (Some(passes), Some(q)) = (read.passes, read_q) {
            let q_level = (q.round().max(0.0) as usize).min(READ_Q_LEVELS - 1);
            self.pacbio.passes_q[(passes as usize).min(60) * READ_Q_LEVELS + q_level] += 1;
        }
        if read.platform == Platform::PacBio {
            if let Some(zmw) = read.zmw {
                *self.pacbio.zmw_reads.entry(zmw).or_default() += 1;
            }
        }
    }

    /// cDNA structure from the primers found at the read's ends.
    fn add_cdna(&mut self, hits: &[Hit], adapters: &Adapters, n: usize) {
        let at = |role: Role, end: bool| {
            hits.iter().any(|h| {
                adapters.role(h.adapter) == role
                    && if end {
                        n - h.end < END_WINDOW
                    } else {
                        h.start < END_WINDOW
                    }
            })
        };
        let ssp_start = at(Role::CdnaSsp, false);
        let ssp_end = at(Role::CdnaSsp, true);
        let vnp_start = at(Role::CdnaVnp, false);
        let vnp_end = at(Role::CdnaVnp, true);
        let c = &mut self.cdna;
        if ssp_start && vnp_end {
            c.full_length_forward += 1;
        } else if vnp_start && ssp_end {
            c.full_length_reverse += 1;
        } else if ssp_start || ssp_end {
            c.ssp_only += 1;
        } else if vnp_start || vnp_end {
            c.vnp_only += 1;
        }
    }

    /// Barcodes left in the read: its own, untrimmed or inside, and other samples'.
    fn add_barcodes(&mut self, read: &Read, index: &Adapters, hits: &[Hit]) {
        let n = read.bases.len();
        let b = &mut self.barcodes;
        b.reads_scanned += 1;
        let label = read.barcode.as_deref().unwrap_or_default();
        let (mut own_end, mut own_inside, mut foreign) = (false, false, false);
        for hit in hits {
            let name = &index.names()[hit.adapter];
            if barcodes::is_own(label, name) {
                if hit.start < END_WINDOW || n - hit.end < END_WINDOW {
                    own_end = true;
                } else {
                    own_inside = true;
                }
            } else {
                foreign = true;
                match b.foreign.get_mut(name) {
                    Some(count) => *count += 1,
                    None => {
                        b.foreign.insert(name.clone(), 1);
                    }
                }
            }
        }
        b.own_at_end += u64::from(own_end);
        b.own_internal += u64::from(own_inside);
        b.with_foreign += u64::from(foreign);
    }

    /// Per-position qualities; returns the read Q (reported, from `rq`, or from the bases).
    fn add_qualities(&mut self, read: &Read) -> Option<f32> {
        let quals = &read.quals;
        let computed = if quals.is_empty() {
            None
        } else {
            self.reads_with_quals += 1;
            let mut error_sum = 0.0;
            for &q in quals {
                let q = usize::from(q).min(QUAL_LEVELS - 1);
                self.base_quals[q] += 1;
                error_sum += ERROR_PROB[q];
            }
            let n = quals.len();
            POSITIONS.for_each_segment(n, |bin, start, end| {
                let row = &mut self.quals_from_start[bin * QUAL_LEVELS..(bin + 1) * QUAL_LEVELS];
                for &q in &quals[start..end] {
                    row[usize::from(q).min(QUAL_LEVELS - 1)] += 1;
                }
                let row = &mut self.quals_from_end[bin * QUAL_LEVELS..(bin + 1) * QUAL_LEVELS];
                for &q in &quals[n - end..n - start] {
                    row[usize::from(q).min(QUAL_LEVELS - 1)] += 1;
                }
            });
            Some(phred(error_sum / n as f64))
        };
        read.reported_q
            .or_else(|| read.accuracy.map(|rq| phred(1.0 - f64::from(rq))))
            .or(computed)
    }

    fn add_sequence(&mut self, read: &Read, settings: &Settings<'_>, scratch: &mut Scratch) {
        let adapters = settings.adapters;
        let bases = &read.bases;
        let n = bases.len();
        self.reads_with_seq += 1;
        let mut counts = [0u64; 5];
        POSITIONS.for_each_segment(n, |bin, start, end| {
            let row = &mut self.bases_from_start[bin * 5..bin * 5 + 5];
            for &b in &bases[start..end] {
                row[usize::from(b)] += 1;
            }
            let row = &mut self.bases_from_end[bin * 5..bin * 5 + 5];
            for &b in &bases[n - end..n - start] {
                row[usize::from(b)] += 1;
                counts[usize::from(b)] += 1;
            }
        });
        let acgt = counts[usize::from(BASE_A)]
            + counts[usize::from(BASE_C)]
            + counts[usize::from(BASE_G)]
            + counts[usize::from(BASE_T)];
        let gc = counts[usize::from(BASE_G)] + counts[usize::from(BASE_C)];
        if let Some(percent) = (gc * 100 + acgt / 2).checked_div(acgt) {
            self.gc[percent as usize] += 1;
            if n as u64 >= GC_PEAK_MIN_LENGTH {
                self.gc_long[percent as usize] += 1;
            }
            scratch.outcome.gc_pct = u8::try_from(percent).ok();
        }

        adapters.scan(bases, &mut scratch.adapter_scan, &mut scratch.hits);
        let mut any_at_end = false;
        let mut any_internal = false;
        for (id, counts) in self.adapters.iter_mut().enumerate() {
            let mut at_start = false;
            let mut at_end = false;
            let mut internal = false;
            let mut nearest_start = None;
            let mut nearest_end = None;
            for hit in scratch.hits.iter().filter(|h| h.adapter == id) {
                let from_end = n - hit.end;
                if hit.start < END_WINDOW {
                    at_start = true;
                } else if from_end < END_WINDOW {
                    at_end = true;
                } else {
                    internal = true;
                }
                nearest_start = Some(nearest_start.map_or(hit.start, |s: usize| s.min(hit.start)));
                nearest_end = Some(nearest_end.map_or(from_end, |e: usize| e.min(from_end)));
            }
            counts.reads_at_start += u64::from(at_start);
            counts.reads_at_end += u64::from(at_end);
            counts.reads_internal += u64::from(internal);
            any_at_end |= at_start || at_end;
            any_internal |= internal;
            if let Some(start) = nearest_start {
                counts.nearest_start[POSITIONS.bin_of(start as u64)] += 1;
            }
            if let Some(end) = nearest_end {
                counts.nearest_end[POSITIONS.bin_of(end as u64)] += 1;
            }
        }
        self.reads_adapter_at_end += u64::from(any_at_end);
        self.reads_adapter_internal += u64::from(any_internal);
        if longest_run(&bases[n.saturating_sub(POLY_WINDOW)..], BASE_A) >= POLY_MIN {
            self.poly_a_tail += 1;
        }
        if longest_run(&bases[..n.min(POLY_WINDOW)], BASE_T) >= POLY_MIN {
            self.poly_t_head += 1;
        }
        self.sample.offer(read.name_hash, bases);
        self.add_cdna(&scratch.hits, adapters, n);

        scratch
            .outcome
            .adapters
            .extend(scratch.hits.iter().map(|h| h.adapter));
        if let Some(index) = settings.barcodes.for_platform(read.platform) {
            index.scan(bases, &mut scratch.adapter_scan, &mut scratch.barcode_hits);
            self.add_barcodes(read, index, &scratch.barcode_hits);
            scratch.outcome.barcodes.extend(
                scratch
                    .barcode_hits
                    .iter()
                    .map(|h| index.names()[h.adapter].clone()),
            );
        }
        if let Some(screen) = settings.screen {
            if self.screen.is_empty() {
                self.screen = vec![Yield::default(); screen.references().len()];
                self.screen_testable = vec![0; screen.references().len()];
            }
            self.reads_screened += 1;
            scratch.outcome.screened = true;
            let found = screen.classify(bases, scratch.outcome.read_q, &mut scratch.screen);
            if let Some(reference) = found {
                self.screen[reference].add(read.len);
                scratch.outcome.screen = Some(reference);
            }
            if Screen::sparse_testable(&scratch.screen) {
                self.sparse_testable += 1;
                if let Some(reference) = found {
                    self.screen_testable[reference] += 1;
                }
            }
        }
        if let Some(kinnex) = settings.kinnex {
            kinnex.scan(bases, &mut scratch.adapter_scan, &mut scratch.kinnex_hits);
            self.kinnex[scratch.kinnex_hits.len().min(KINNEX_MAX)] += 1;
            scratch.outcome.kinnex_adapters = Some(scratch.kinnex_hits.len());
        }
    }

    fn add_run_metadata(&mut self, read: &Read, read_q: Option<f32>, filter: Option<Filter>) {
        if let Some(start) = read.start_secs {
            let bin = self.time.entry((start / 60.0).floor() as i64).or_default();
            bin.reads += 1;
            bin.bases += read.len;
            bin.pass_reads += u64::from(filter == Some(Filter::Pass));
            if let Some(q) = read_q {
                bin.q_tenths += (f64::from(q) * 10.0).round().max(0.0) as u64;
                bin.q_reads += 1;
            }
        }
        if let Some(channel) = read.channel {
            self.channels.entry(channel).or_default().add(read.len);
        }
        if let Some(duration) = read.duration_secs.filter(|&d| d > 0.0) {
            let speed = read.len as f64 / f64::from(duration);
            self.speed[((speed / SPEED_BIN_WIDTH) as usize).min(SPEED_BINS - 1)] += 1;
        }
        if let Some(reason) = &read.end_reason {
            match self.end_reasons.get_mut(reason.as_ref()) {
                Some(count) => *count += 1,
                None => {
                    self.end_reasons.insert(Arc::clone(reason), 1);
                }
            }
        }
    }

    fn add_pacbio(&mut self, read: &Read) {
        if let Some(passes) = read.passes {
            self.passes[(passes as usize).min(100)] += 1;
        }
        if let Some(coverage) = read.effective_coverage {
            self.effective_coverage[(coverage.max(0.0) as usize).min(100)] += 1;
        }
        if let Some(score) = read.barcode_score {
            self.barcode_score[(score as usize).min(100)] += 1;
        }
        if read.kinetics == Kinetics::Present {
            self.with_kinetics += 1;
        }
    }

    /// Counts the calls in an `MM`/`ML` pair. Each `MM` entry (`C+m?,0,3,1`) has one delta
    /// per call, counting skipped occurrences of its base in the read as sequenced, and, for
    /// combined codes such as `C+mh`, one `ML` value per code per call. With the bases at
    /// hand, calls on a C followed by G are also counted as CpG calls.
    ///
    /// In implicit mode (`C+m.`, or no mode, as PacBio's jasmine writes) the bases an entry
    /// skips are unmodified calls, and they are counted as such; in explicit mode (`C+m?`,
    /// as dorado writes) they were not called.
    fn add_mods(&mut self, positions: &[u8], probs: &[u8], bases: &[u8]) {
        self.reads_with_mods += 1;
        let mut probs = probs.iter().copied();
        for entry in positions.split(|&b| b == b';').filter(|e| !e.is_empty()) {
            let mut parts = entry.split(|&b| b == b',');
            let head = parts.next().unwrap_or_default();
            if head.len() < 3 {
                continue;
            }
            let (base_strand, codes) = head.split_at(2);
            let explicit = codes.ends_with(b"?");
            let codes = codes
                .strip_suffix(b"?")
                .or_else(|| codes.strip_suffix(b"."))
                .unwrap_or(codes);
            let codes: Vec<&[u8]> = if codes.iter().all(u8::is_ascii_digit) {
                vec![codes]
            } else {
                codes.chunks(1).collect()
            };
            // Counted locally and merged once per entry: a read carries thousands of calls.
            let mut local = vec![ModCounts::default(); codes.len()];
            // The read positions of the entry's base: deltas count them.
            let Some(target) = mod_base(base_strand) else {
                continue;
            };
            let mut occurrences = bases
                .iter()
                .enumerate()
                .filter(move |&(_, &b)| target.is_none_or(|t| b == t))
                .map(|(i, _)| i);
            let have_bases = !bases.is_empty();
            // CpG context is defined for calls on C of the sequenced strand.
            let cpg_context = base_strand == b"C+" && have_bases;
            let is_cpg = |i: usize| cpg_context.then(|| bases.get(i + 1) == Some(&BASE_G));
            let implicit = !explicit && have_bases;
            let unmodified = |local: &mut [ModCounts], cpg: Option<bool>| {
                for counts in local.iter_mut() {
                    counts.calls += 1;
                    counts.probs[0] += 1;
                    if cpg == Some(true) {
                        counts.cpg_calls += 1;
                    }
                }
            };
            for delta in parts.filter_map(parse_delta) {
                let cpg = if implicit {
                    for _ in 0..delta {
                        match occurrences.next() {
                            Some(i) => unmodified(&mut local, is_cpg(i)),
                            None => break,
                        }
                    }
                    occurrences.next().and_then(is_cpg)
                } else if have_bases {
                    occurrences.nth(delta).and_then(is_cpg)
                } else {
                    None
                };
                for counts in &mut local {
                    let Some(prob) = probs.next() else { break };
                    let modified = u64::from(prob >= MOD_THRESHOLD);
                    counts.calls += 1;
                    counts.modified += modified;
                    if cpg == Some(true) {
                        counts.cpg_calls += 1;
                        counts.cpg_modified += modified;
                    }
                    counts.probs[usize::from(prob) / 10] += 1;
                }
            }
            if implicit {
                for i in occurrences {
                    unmodified(&mut local, is_cpg(i));
                }
            }
            for (code, counts) in codes.iter().zip(local).filter(|(_, c)| c.calls > 0) {
                let key = format!(
                    "{}{}",
                    String::from_utf8_lossy(base_strand),
                    String::from_utf8_lossy(code)
                );
                let total = self.mods.entry(key).or_default();
                total.calls += counts.calls;
                total.modified += counts.modified;
                total.cpg_calls += counts.cpg_calls;
                total.cpg_modified += counts.cpg_modified;
                add_into(&mut total.probs, &counts.probs);
            }
        }
    }

    /// Counts the alignment and returns the read's identity, when it has one.
    fn add_alignment(&mut self, read: &Read) -> Option<f64> {
        let mut read_identity = None;
        match read.alignment {
            Alignment::NotAligned => {}
            Alignment::Unmapped => self.mapping.unmapped += 1,
            Alignment::Mapped(mapped) => {
                let m = &mut self.mapping;
                m.mapped += 1;
                m.mapped_bases += read.len;
                m.split += u64::from(mapped.split);
                m.mapq[usize::from(mapped.mapq).min(60)] += 1;
                let columns = mapped.aligned_columns + mapped.inserted + mapped.deleted;
                if let (Some(nm), true) = (mapped.edit_distance, columns > 0) {
                    let nm = u64::from(nm);
                    let identity = 1.0 - nm.min(columns) as f64 / columns as f64;
                    m.identity[identity_bin(identity)] += 1;
                    read_identity = Some(identity);
                    let mismatches = nm.saturating_sub(mapped.inserted + mapped.deleted);
                    let gaps = mapped.insertion_runs + mapped.deletion_runs;
                    let compressed = 1.0
                        - (mismatches + gaps) as f64
                            / (mapped.aligned_columns + gaps).max(1) as f64;
                    m.gap_compressed_identity[identity_bin(compressed)] += 1;
                }
                if read.len > 0 {
                    let aligned = (mapped.aligned_columns + mapped.inserted).min(read.len);
                    m.aligned_fraction[(aligned * 100 / read.len) as usize] += 1;
                }
            }
        }
        read_identity
    }

    /// Adds `other`'s counts to these.
    pub fn merge(&mut self, other: Stats) {
        for (a, b) in self.platform_reads.iter_mut().zip(other.platform_reads) {
            *a += b;
        }
        self.reads += other.reads;
        self.bases += other.bases;
        self.lengths.extend_from_slice(&other.lengths);
        self.read_q_tenths.extend_from_slice(&other.read_q_tenths);
        let (o, p) = (&mut self.ont, &other.ont);
        for (a, b) in o.duplex.iter_mut().zip(p.duplex) {
            a.merge(b);
        }
        add_into(&mut o.duplex_q, &p.duplex_q);
        o.split_reads.merge(p.split_reads);
        for (&mux, &y) in &p.mux {
            o.mux.entry(mux).or_default().merge(y);
        }
        for (mine, theirs) in [
            (&mut o.shift, &p.shift),
            (&mut o.scale, &p.scale),
            (&mut o.samples_per_base, &p.samples_per_base),
        ] {
            for (&bin, &count) in theirs {
                *mine.entry(bin).or_default() += count;
            }
        }
        add_into(&mut o.poly_a, &p.poly_a);
        add_into(&mut self.pacbio.passes_q, &other.pacbio.passes_q);
        for (&zmw, &count) in &other.pacbio.zmw_reads {
            *self.pacbio.zmw_reads.entry(zmw).or_default() += count;
        }
        let (c, d) = (&mut self.cdna, other.cdna);
        c.full_length_forward += d.full_length_forward;
        c.full_length_reverse += d.full_length_reverse;
        c.ssp_only += d.ssp_only;
        c.vnp_only += d.vnp_only;
        let (b, e) = (&mut self.barcodes, &other.barcodes);
        b.reads_scanned += e.reads_scanned;
        b.own_at_end += e.own_at_end;
        b.own_internal += e.own_internal;
        b.with_foreign += e.with_foreign;
        for (name, &count) in &e.foreign {
            *b.foreign.entry(name.clone()).or_default() += count;
        }
        if self.screen.is_empty() {
            self.screen.clone_from(&other.screen);
        } else {
            for (a, &y) in self.screen.iter_mut().zip(&other.screen) {
                a.merge(y);
            }
        }
        self.reads_screened += other.reads_screened;
        self.sparse_testable += other.sparse_testable;
        if self.screen_testable.is_empty() {
            self.screen_testable.clone_from(&other.screen_testable);
        } else {
            add_into(&mut self.screen_testable, &other.screen_testable);
        }
        add_into(&mut self.kinnex, &other.kinnex);
        self.rna_reads += other.rna_reads;
        add_into(&mut self.read_q, &other.read_q);
        add_into(&mut self.reads_at_q, &other.reads_at_q);
        self.bases_at_q20 += other.bases_at_q20;
        add_into(&mut self.read_q_bases, &other.read_q_bases);
        self.reads_with_q += other.reads_with_q;
        add_into(&mut self.length_q, &other.length_q);
        self.reads_with_seq += other.reads_with_seq;
        self.reads_with_quals += other.reads_with_quals;
        add_into(&mut self.quals_from_start, &other.quals_from_start);
        add_into(&mut self.quals_from_end, &other.quals_from_end);
        add_into(&mut self.bases_from_start, &other.bases_from_start);
        add_into(&mut self.bases_from_end, &other.bases_from_end);
        add_into(&mut self.base_quals, &other.base_quals);
        add_into(&mut self.gc, &other.gc);
        add_into(&mut self.gc_long, &other.gc_long);
        self.pass.merge(other.pass);
        self.fail.merge(other.fail);
        self.reads_adapter_at_end += other.reads_adapter_at_end;
        self.reads_adapter_internal += other.reads_adapter_internal;
        for (a, b) in self.adapters.iter_mut().zip(&other.adapters) {
            a.merge(b);
        }
        self.poly_a_tail += other.poly_a_tail;
        self.poly_t_head += other.poly_t_head;
        self.sample.merge(other.sample);
        self.duplicate_names += other.duplicate_names;
        self.pacbio_q_estimated += other.pacbio_q_estimated;
        // Fold the smaller set into the larger.
        let mine = std::mem::take(&mut self.name_hashes);
        let (mut larger, smaller) = if mine.len() >= other.name_hashes.len() {
            (mine, other.name_hashes)
        } else {
            (other.name_hashes, mine)
        };
        for hash in smaller {
            if !larger.insert(hash) {
                self.duplicate_names += 1;
            }
        }
        self.name_hashes = larger;
        for (minute, bin) in other.time {
            let mine = self.time.entry(minute).or_default();
            mine.reads += bin.reads;
            mine.bases += bin.bases;
            mine.pass_reads += bin.pass_reads;
            mine.q_tenths += bin.q_tenths;
            mine.q_reads += bin.q_reads;
        }
        for (channel, y) in other.channels {
            self.channels.entry(channel).or_default().merge(y);
        }
        add_into(&mut self.speed, &other.speed);
        for (reason, count) in other.end_reasons {
            *self.end_reasons.entry(reason).or_default() += count;
        }
        add_into(&mut self.passes, &other.passes);
        add_into(&mut self.effective_coverage, &other.effective_coverage);
        add_into(&mut self.barcode_score, &other.barcode_score);
        self.with_kinetics += other.with_kinetics;
        self.reads_with_mods += other.reads_with_mods;
        for (code, counts) in other.mods {
            let mine = self.mods.entry(code).or_default();
            mine.calls += counts.calls;
            mine.modified += counts.modified;
            mine.cpg_calls += counts.cpg_calls;
            mine.cpg_modified += counts.cpg_modified;
            add_into(&mut mine.probs, &counts.probs);
        }
        let (m, o) = (&mut self.mapping, other.mapping);
        m.mapped += o.mapped;
        m.unmapped += o.unmapped;
        m.secondary += o.secondary;
        m.supplementary += o.supplementary;
        m.split += o.split;
        m.mapped_bases += o.mapped_bases;
        add_into(&mut m.mapq, &o.mapq);
        add_into(&mut m.identity, &o.identity);
        add_into(&mut m.gap_compressed_identity, &o.gap_compressed_identity);
        add_into(&mut m.aligned_fraction, &o.aligned_fraction);
    }
}

/// The base an `MM` entry's deltas count, in the read as sequenced: the named base for
/// `X+`, its complement for `X-` (calls on the other strand), `None` for `N` (every base).
/// `None` overall for an unknown base letter.
fn mod_base(base_strand: &[u8]) -> Option<Option<u8>> {
    let base = match base_strand.first()?.to_ascii_uppercase() {
        b'A' => BASE_A,
        b'C' => BASE_C,
        b'G' => BASE_G,
        b'T' | b'U' => BASE_T,
        b'N' => return Some(None),
        _ => return None,
    };
    Some(Some(if base_strand.get(1) == Some(&b'-') {
        BASE_T - base
    } else {
        base
    }))
}

/// A PacBio CCS read name, `<movie>/<zmw>/ccs` (or `/ccs/fwd`, `/ccs/rev` by strand).
fn is_ccs_name(name: &[u8]) -> bool {
    name.ends_with(b"/ccs") || name.windows(5).any(|w| w == b"/ccs/")
}

/// An `MM` delta: decimal digits, perhaps padded with spaces.
fn parse_delta(text: &[u8]) -> Option<usize> {
    let digits = text.trim_ascii();
    if digits.is_empty() {
        return None;
    }
    digits.iter().try_fold(0usize, |n, &d| {
        d.is_ascii_digit()
            .then(|| n.checked_mul(10)?.checked_add(usize::from(d - b'0')))
            .flatten()
    })
}

fn platform_index(platform: Platform) -> usize {
    match platform {
        Platform::Ont => 0,
        Platform::PacBio => 1,
        Platform::Unknown => 2,
    }
}

/// The platform most of a group's reads came from.
pub fn majority_platform(platform_reads: &[u64; 3]) -> Platform {
    if platform_reads.iter().all(|&n| n == 0) {
        return Platform::Unknown;
    }
    let platforms = [Platform::Ont, Platform::PacBio, Platform::Unknown];
    let (index, _) = platform_reads
        .iter()
        .enumerate()
        .max_by_key(|&(i, &n)| (n, std::cmp::Reverse(i)))
        .expect("three platforms");
    platforms[index]
}

fn add_into(a: &mut [u64], b: &[u64]) {
    for (x, y) in a.iter_mut().zip(b) {
        *x += y;
    }
}

fn phred(error: f64) -> f32 {
    (-10.0 * error.max(1e-6).log10()) as f32
}

fn identity_bin(identity: f64) -> usize {
    ((identity.clamp(0.0, 1.0) * 1000.0).round() as usize).min(IDENTITY_BINS - 1)
}

fn longest_run(bases: &[u8], base: u8) -> usize {
    let mut best = 0;
    let mut run = 0;
    for &b in bases {
        run = if b == base { run + 1 } else { 0 };
        best = best.max(run);
    }
    best
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::read::Mapped;

    static BARCODES: LazyLock<Barcodes> = LazyLock::new(|| Barcodes::new(&[]));

    fn settings(adapters: &Adapters) -> Settings<'_> {
        Settings {
            min_q: None,
            adapters,
            barcodes: &BARCODES,
            screen: None,
            kinnex: None,
        }
    }

    fn read(bases: &[u8], quals: &[u8], name_hash: u64) -> Read {
        Read {
            name_hash,
            len: bases.len() as u64,
            bases: bases.to_vec(),
            quals: quals.to_vec(),
            ..Read::default()
        }
    }

    #[test]
    fn counts_q_thresholds_exactly() {
        let adapters = Adapters::load(None).expect("built-ins");
        let mut stats = Stats::new(adapters.names().len());
        let mut scratch = Scratch::default();
        for (hash, q) in [(1, 19.96f32), (2, 20.0), (3, 29.99)] {
            let mut r = read(&[0; 100], &[], hash);
            r.reported_q = Some(q);
            stats.add(&r, &settings(&adapters), &mut scratch);
        }
        // Q19.96 rounds to the Q20.0 bin, but is not Q20.
        assert_eq!(stats.reads_at_q, [3, 2, 0]);
        assert_eq!(stats.bases_at_q20, 200);
        assert_eq!(stats.read_q_tenths, vec![199, 200, 299]);
    }

    #[test]
    fn counts_qualities_and_composition() {
        let adapters = Adapters::load(None).expect("built-ins");
        let mut stats = Stats::new(adapters.names().len());
        let mut scratch = Scratch::default();
        stats.add(
            &read(&[0, 1, 2, 3], &[10, 10, 10, 10], 1),
            &settings(&adapters),
            &mut scratch,
        );
        assert_eq!(stats.reads, 1);
        assert_eq!(stats.base_quals[10], 4);
        assert_eq!(stats.quals_from_start[10], 1, "position 1 has Q10");
        assert_eq!(
            stats.quals_from_end[3 * QUAL_LEVELS + 10],
            1,
            "4th base from the end"
        );
        assert_eq!(stats.bases_from_start[..5], [1, 0, 0, 0, 0]);
        assert_eq!(stats.bases_from_end[..5], [0, 0, 0, 1, 0]);
        assert_eq!(stats.gc[50], 1);
        assert_eq!(stats.read_q[100], 1, "all Q10 bases make a Q10 read");
        assert_eq!(stats.pass.reads, 1, "Q10 passes the ONT default");
    }

    #[test]
    fn keeps_short_reads_out_of_the_gc_peak_histogram() {
        let adapters = Adapters::load(None).expect("built-ins");
        let mut stats = Stats::new(adapters.names().len());
        let mut scratch = Scratch::default();
        // A 60-base adapter fragment and a 600-base read, both 50% GC.
        stats.add(
            &read(&[1, 0].repeat(30), &[], 1),
            &settings(&adapters),
            &mut scratch,
        );
        stats.add(
            &read(&[2, 3].repeat(300), &[], 2),
            &settings(&adapters),
            &mut scratch,
        );
        assert_eq!(stats.gc[50], 2);
        assert_eq!(stats.gc_long[50], 1);
    }

    #[test]
    fn counts_duplicate_names_across_merges() {
        let adapters = Adapters::load(None).expect("built-ins");
        let mut a = Stats::new(adapters.names().len());
        let mut b = Stats::new(adapters.names().len());
        let mut scratch = Scratch::default();
        a.add(&read(&[0; 10], &[], 7), &settings(&adapters), &mut scratch);
        a.add(&read(&[0; 10], &[], 8), &settings(&adapters), &mut scratch);
        b.add(&read(&[0; 10], &[], 7), &settings(&adapters), &mut scratch);
        a.merge(b);
        assert_eq!(a.reads, 3);
        assert_eq!(a.duplicate_names, 1);
    }

    #[test]
    fn counts_mod_calls() {
        let adapters = Adapters::load(None).expect("built-ins");
        let mut stats = Stats::new(adapters.names().len());
        // Read CGACCG: C at 0 (CpG), 3 (not), 4 (CpG); A at 2.
        let bases = [BASE_C, BASE_G, BASE_A, BASE_C, BASE_C, BASE_G];
        stats.add_mods(b"C+mh?,0,1;A+a.,0;", &[200, 10, 100, 100, 255], &bases);
        assert_eq!(stats.mods["C+m"].calls, 2);
        assert_eq!(stats.mods["C+m"].modified, 1);
        assert_eq!(stats.mods["C+h"].calls, 2);
        assert_eq!(stats.mods["A+a"].modified, 1);
        // The second call skips one C (at 3) and lands on the C at 4, a CpG.
        assert_eq!(stats.mods["C+m"].cpg_calls, 2);
        assert_eq!(stats.mods["C+m"].cpg_modified, 1);
        stats.add_mods(b"C+m?,1;", &[255], &bases);
        assert_eq!(stats.mods["C+m"].cpg_calls, 2, "the C at 3 is not a CpG");
    }

    #[test]
    fn counts_implicit_mod_skips_as_unmodified() {
        let adapters = Adapters::load(None).expect("built-ins");
        let mut stats = Stats::new(adapters.names().len());
        // Read CGACCGCG: C at 0 (CpG), 3 (not), 4 (CpG), 6 (CpG).
        let bases = [
            BASE_C, BASE_G, BASE_A, BASE_C, BASE_C, BASE_G, BASE_C, BASE_G,
        ];
        // Implicit: one call, on the C at 4; the other three C's are unmodified calls.
        stats.add_mods(b"C+m.,2;", &[255], &bases);
        let m = &stats.mods["C+m"];
        assert_eq!((m.calls, m.modified), (4, 1));
        assert_eq!((m.cpg_calls, m.cpg_modified), (3, 1));
        // No mode is implicit too (jasmine); `?` is explicit, skipping nothing.
        let mut stats = Stats::new(adapters.names().len());
        stats.add_mods(b"C+m,2;", &[255], &bases);
        assert_eq!(stats.mods["C+m"].calls, 4);
        let mut stats = Stats::new(adapters.names().len());
        stats.add_mods(b"C+m?,2;", &[255], &bases);
        assert_eq!(stats.mods["C+m"].calls, 1);
        assert_eq!(stats.mods["C+m"].cpg_modified, 1);
    }

    #[test]
    fn measures_identity() {
        let adapters = Adapters::load(None).expect("built-ins");
        let mut stats = Stats::new(adapters.names().len());
        let mut r = read(&[0; 100], &[], 1);
        r.alignment = Alignment::Mapped(Mapped {
            mapq: 60,
            edit_distance: Some(5),
            aligned_columns: 98,
            inserted: 2,
            deleted: 0,
            insertion_runs: 1,
            deletion_runs: 0,
            split: false,
        });
        stats.add(&r, &settings(&adapters), &mut Scratch::default());
        // 5 edits over 100 columns: 95.0%.
        assert_eq!(stats.mapping.identity[950], 1);
        // 3 mismatches plus one gap over 99: 95.96%.
        assert_eq!(stats.mapping.gap_compressed_identity[960], 1);
        assert_eq!(stats.mapping.aligned_fraction[100], 1);
    }
}
