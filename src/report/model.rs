//! Turning a group's accumulated counts into summary numbers and plot data.

use std::collections::HashMap;

use serde::{Deserialize, Serialize};

use super::checks::{self, Module};
use super::thresholds::Thresholds;
use super::GroupKind;
use crate::adapters::Adapters;
use crate::bins::{self, POSITIONS};
use crate::read::{Molecule, Platform, BASE_N};
use crate::sampler::{self, Ends};
use crate::screen::Reference;
use crate::stats::{self, Stats, QUAL_LEVELS, READ_Q_LEVELS};

/// A per-position bin is plotted while at least this fraction of reads reach it.
const MIN_POSITION_COVERAGE: f64 = 0.005;

/// An end sequence is overrepresented at this fraction of the sample (FastQC's 0.1%).
pub const OVERREPRESENTED_FRACTION: f64 = 0.001;

/// Overrepresented sequences listed per end.
const OVERREPRESENTED_LISTED: usize = 20;

/// Run-wide inputs to building a group.
#[derive(Debug)]
pub struct Context<'run> {
    pub adapters: &'run [String],
    pub adapter_index: &'run Adapters,
    pub min_q: Option<f32>,
    pub total_reads: u64,
    pub thresholds: &'run Thresholds,
    /// Haploid genome size, for coverage.
    pub genome_size: Option<u64>,
    /// Sample names by barcode, from a sample sheet.
    pub samples: &'run HashMap<String, String>,
    pub screen: Option<&'run [Reference]>,
    /// The first minute of each run (keyed as `run_of` names runs), so every group's yield
    /// over time starts at its run's start, not at its own first read.
    pub run_starts: &'run HashMap<String, i64>,
}

/// One group of reads: all reads, or one barcode.
#[derive(Debug, Serialize, Deserialize)]
pub struct Group {
    pub name: String,
    /// Sample name from the sample sheet.
    pub sample: Option<String>,
    pub kind: GroupKind,
    pub platform: Platform,
    pub molecule: Molecule,
    pub summary: Summary,
    pub modules: Vec<Module>,
    pub data: Data,
    /// Read Q warn and fail levels, as the checks and the charts' bands use them.
    #[serde(default)]
    pub quality_levels: Option<(f64, f64)>,
}

/// Headline numbers of a group.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct Summary {
    pub reads: u64,
    pub bases: u64,
    /// Percent of the run's reads.
    pub share_pct: f64,
    pub mean_length: f64,
    /// The middle length; for an even number of reads, the mean of the two middle ones.
    pub median_length: f64,
    pub n50: u64,
    pub n90: u64,
    pub shortest: u64,
    pub longest: u64,
    pub median_read_q: Option<f64>,
    pub mean_read_q: Option<f64>,
    pub reads_q10_pct: Option<f64>,
    pub reads_q20_pct: Option<f64>,
    pub reads_q30_pct: Option<f64>,
    /// Bases in reads of Q20 or more.
    pub bases_q20: Option<u64>,
    /// Base Q from the mean error probability of every base, as read Q is computed.
    pub base_q: Option<f64>,
    pub gc_pct: Option<f64>,
    pub pass_reads: u64,
    pub pass_bases: u64,
    pub pass_pct: Option<f64>,
    pub mapped_pct: Option<f64>,
    pub median_identity_pct: Option<f64>,
    pub duplicate_names: u64,
    pub mod_reads_pct: Option<f64>,
    /// Percent of 5mC calls (`C+m`) at CpGs called modified.
    pub methylated_pct: Option<f64>,
    /// Bases over the genome size, when `--genome-size` is given.
    pub coverage: Option<f64>,
    pub pass_coverage: Option<f64>,
}

/// Plot data of a group.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct Data {
    pub length_bins: Vec<LengthBin>,
    pub read_q: Vec<QBin>,
    /// Reads by log length bin and whole read Q, nonzero cells only.
    pub length_q: Vec<(usize, usize, u64)>,
    pub quality_from_start: Vec<PositionQuality>,
    pub quality_from_end: Vec<PositionQuality>,
    pub content_from_start: Vec<PositionContent>,
    pub content_from_end: Vec<PositionContent>,
    /// Bases by Phred value.
    pub base_quals: Vec<u64>,
    /// Reads by GC percent, 0 to 100.
    pub gc: Vec<u64>,
    /// A normal distribution with the observed mean and deviation, for comparison.
    pub gc_expected: Vec<f64>,
    /// Summed departure from that curve, as a percent of reads (FastQC's measure).
    pub gc_deviation_pct: Option<f64>,
    /// A separate GC peak beside the main one, if any.
    pub gc_second_peak: Option<GcPeak>,
    pub duplication: Option<Duplication>,
    pub overrepresented: Vec<Overrepresented>,
    pub adapters: Vec<AdapterSummary>,
    /// Reads with any adapter or primer near either end, and inside, as percents.
    #[serde(default)]
    pub adapter_at_ends_pct: f64,
    #[serde(default)]
    pub adapter_internal_pct: f64,
    pub poly_a_tail_pct: Option<f64>,
    pub poly_t_head_pct: Option<f64>,
    pub time: Vec<TimePoint>,
    pub channels: Vec<ChannelYield>,
    /// Reads by translocation speed, 10 bases per second per bin.
    pub speed: Vec<u64>,
    pub end_reasons: Vec<(String, u64)>,
    pub pacbio: Option<PacBio>,
    pub mods: Vec<ModSummary>,
    pub mapping: Option<MappingSummary>,
    pub filter_grid: Option<FilterGrid>,
    pub contaminants: Vec<ScreenHit>,
    pub reads_screened: u64,
    pub barcode_scan: Option<BarcodeScan>,
    pub ont: Option<OntExtras>,
    pub cdna: Option<CdnaSummary>,
    /// Reads by dorado poly(A) estimate, in bins of 5 bases.
    pub poly_a: Vec<u64>,
    /// Reads by Kinnex adapters found.
    pub kinnex: Vec<u64>,
}

/// Reads and bases left after filtering by minimum length (rows) and read Q (columns).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FilterGrid {
    pub min_lengths: Vec<u64>,
    pub min_qs: Vec<f64>,
    pub reads: Vec<Vec<u64>>,
    pub bases: Vec<Vec<u64>>,
}

/// Reads assigned to one screen reference.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ScreenHit {
    pub reference: Reference,
    pub reads: u64,
    pub bases: u64,
    /// Share of reads screened, or for a sparse reference with `testable`, of those reads.
    pub reads_pct: f64,
    pub bases_pct: f64,
    /// For a sparse reference (human) whose share is of the reads long enough to test: how
    /// many there were.
    #[serde(default)]
    pub testable: Option<u64>,
    /// For a sparse reference whose share is of all reads because too few were long enough
    /// to test: a lower bound, as the short reads were not tested.
    #[serde(default)]
    pub lower_bound: bool,
}

/// A sparse reference's share is of the reads long enough to test only when there are at
/// least this many of them...
const MIN_TESTABLE_READS: u64 = 100;

/// ...and they are at least this fraction of the reads screened. Their share assumes the
/// sparse reference's reads are no longer than the rest: in a 1.5 kb amplicon run with 2%
/// human reads of 10 kb, the human reads are the only ones long enough to test, and their
/// share of those (100%) would call the sample human.
const MIN_TESTABLE_FRACTION: f64 = 0.2;

/// Barcodes found inside reads.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BarcodeScan {
    pub scanned: u64,
    pub own_at_end_pct: f64,
    pub own_internal_pct: f64,
    pub foreign_pct: f64,
    /// The most frequent foreign barcodes: name, reads, percent of reads.
    pub foreign: Vec<(String, u64, f64)>,
}

/// ONT read properties beyond the basics.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OntExtras {
    pub duplex: Option<DuplexSummary>,
    pub split_reads: u64,
    pub split_pct: f64,
    /// Reads and bases by mux.
    pub mux: Vec<(u8, u64, u64)>,
    /// Reads by signal shift, signal scale and samples per base: bin start and count.
    pub shift: Vec<(f64, u64)>,
    pub scale: Vec<(f64, u64)>,
    pub samples_per_base: Vec<(f64, u64)>,
}

/// dorado duplex classes.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DuplexSummary {
    /// Class name, reads, bases, median read Q.
    pub classes: Vec<(String, u64, u64, Option<f64>)>,
    /// Bases in duplex reads, twice over (each covers two strands), as a percent of
    /// all simplex bases (with and without duplex offspring): ONT's duplex rate.
    pub duplex_rate_pct: f64,
}

/// cDNA primer structure.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct CdnaSummary {
    pub full_length_pct: f64,
    pub forward_pct: f64,
    pub reverse_pct: f64,
    pub ssp_only_pct: f64,
    pub vnp_only_pct: f64,
}

/// Reads and bases in one log length bin.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct LengthBin {
    pub start: f64,
    pub end: f64,
    pub reads: u64,
    pub bases: u64,
}

/// Reads and bases at one read Q (in tenths).
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct QBin {
    pub q: f64,
    pub reads: u64,
    pub bases: u64,
}

/// Quality distribution of one position bin.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct PositionQuality {
    pub bin: usize,
    /// Reads reaching the bin, estimated from its bases.
    pub reads: f64,
    pub mean: f64,
    pub p10: f64,
    pub p25: f64,
    pub median: f64,
    pub p75: f64,
    pub p90: f64,
}

/// Base composition of one position bin, in percent.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct PositionContent {
    pub bin: usize,
    /// Bases counted in the bin.
    pub bases: u64,
    pub a: f64,
    pub c: f64,
    pub g: f64,
    pub t: f64,
    pub n: f64,
}

/// A secondary peak of the per-read GC distribution.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct GcPeak {
    pub gc: usize,
    pub main_gc: usize,
    /// Percent of reads within 5 GC points of the peak.
    pub share_pct: f64,
    /// The peak was judged on reads at least this long; 0 for all reads.
    pub min_length: u64,
}

/// Duplication of sampled reads, judged by both 50-base ends.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Duplication {
    pub sampled: usize,
    /// Percent of sampled reads left after removing duplicates.
    pub distinct_pct: f64,
    /// Percent of sampled reads at each duplication level (FastQC's levels).
    pub levels: Vec<(String, f64)>,
}

/// A 50-base read end seen in more than 0.1% of sampled reads.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Overrepresented {
    pub end: ReadEnd,
    pub sequence: String,
    pub count: u64,
    pub pct: f64,
    pub source: String,
}

/// Which end of a read.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ReadEnd {
    Start,
    End,
}

/// Where one adapter was found.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AdapterSummary {
    pub name: String,
    pub start_pct: f64,
    pub end_pct: f64,
    pub internal_pct: f64,
    /// Percent of reads with the adapter at or before each position bin, from each end.
    pub cumulative_from_start: Vec<f64>,
    pub cumulative_from_end: Vec<f64>,
}

/// Yield in one slice of run time.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct TimePoint {
    pub hours: f64,
    pub reads: u64,
    pub bases: u64,
    pub cumulative_bases: u64,
    pub mean_q: Option<f64>,
    /// `None` for a slice without reads (or without pass/fail).
    pub pass_pct: Option<f64>,
    pub mean_length: Option<f64>,
}

/// Yield of one channel.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct ChannelYield {
    pub channel: u32,
    pub reads: u64,
    pub bases: u64,
}

/// PacBio-specific distributions.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PacBio {
    pub passes: Vec<u64>,
    pub effective_coverage: Vec<u64>,
    pub barcode_score: Vec<u64>,
    pub median_passes: Option<u64>,
    pub kinetics_pct: f64,
    /// Reads by passes and whole read Q, nonzero cells only.
    pub passes_q: Vec<(usize, usize, u64)>,
    /// ZMWs by reads they gave (1, 2, ... with the last bin holding more).
    pub reads_per_zmw: Vec<u64>,
    pub zmws: u64,
}

/// One base modification code.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ModSummary {
    pub code: String,
    pub label: String,
    pub calls: u64,
    pub modified_pct: f64,
    /// Calls on a CpG, and the percent of those called modified.
    pub cpg_calls: u64,
    pub cpg_modified_pct: Option<f64>,
    pub calls_per_read: f64,
    /// Calls by ML probability bin.
    pub probs: Vec<u64>,
}

/// Alignment summary of aligned input.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MappingSummary {
    pub mapped: u64,
    pub unmapped: u64,
    pub secondary: u64,
    pub supplementary: u64,
    pub mapped_pct: f64,
    pub split_pct: f64,
    pub mapq: Vec<u64>,
    /// Reads by identity, in 0.1% bins from `identity_floor`.
    pub identity: Vec<u64>,
    pub gap_compressed_identity: Vec<u64>,
    pub identity_floor: f64,
    pub median_identity_pct: Option<f64>,
    pub median_gap_compressed_pct: Option<f64>,
    pub aligned_fraction: Vec<u64>,
}

impl Group {
    /// Summarizes a group's statistics and runs the checks.
    pub fn build(name: &str, kind: GroupKind, stats: &Stats, context: &Context<'_>) -> Self {
        let platform = stats::majority_platform(&stats.platform_reads);
        let molecule = if stats.rna_reads * 2 > stats.reads {
            Molecule::Rna
        } else {
            Molecule::Dna
        };
        let mut summary = summarize(stats, context.total_reads, context.genome_size);
        let mut data = Data::build(stats, context);
        let run_start = context.run_starts.get(super::run_of(name)).copied();
        data.time = time_points(stats, run_start);
        if molecule == Molecule::Rna {
            // CpG methylation is a DNA measure; in RNA, C+m is m5C in any context.
            summary.methylated_pct = None;
            for m in &mut data.mods {
                m.cpg_calls = 0;
                m.cpg_modified_pct = None;
            }
        }
        for m in &mut data.mods {
            m.label = mod_label(&m.code, molecule);
        }
        let modules = checks::evaluate(
            stats,
            &summary,
            &data,
            platform,
            molecule,
            context.min_q,
            context.thresholds,
        );
        // Group names are a barcode, or `run · barcode` when runs are split.
        let barcode = name.rsplit(" · ").next().unwrap_or(name);
        Self {
            name: name.to_owned(),
            sample: sample_name(context.samples, barcode),
            kind,
            platform,
            molecule,
            summary,
            modules,
            data,
            quality_levels: Some(checks::quality_levels(
                platform,
                context.min_q,
                context.thresholds,
            )),
        }
    }

    /// The status of one module.
    pub fn module(&self, id: checks::ModuleId) -> Option<&Module> {
        self.modules.iter().find(|m| m.id == id)
    }
}

/// A group's sample name. Sample sheets name PacBio barcodes singly (`bc2052`) while lima's
/// `bc` tags name the pair (`bc2052--bc2052`), so a symmetric pair also matches its barcode.
fn sample_name(samples: &HashMap<String, String>, barcode: &str) -> Option<String> {
    samples
        .get(barcode)
        .or_else(|| {
            let (forward, reverse) = barcode.split_once("--")?;
            (forward == reverse).then(|| samples.get(forward)).flatten()
        })
        .cloned()
}

/// Phred Q of the mean error probability over a histogram of base Qs, or `None` when empty.
///
/// Averaging Phred values directly overstates quality: it is a log scale, so the many
/// confident bases outweigh the few likely errors.
fn error_averaged_q(counts: &[u64]) -> Option<f64> {
    let total: u64 = counts.iter().sum();
    if total == 0 {
        return None;
    }
    let errors: f64 = counts
        .iter()
        .enumerate()
        .map(|(q, &c)| c as f64 * 10f64.powf(-(q as f64) / 10.0))
        .sum();
    Some(-10.0 * (errors / total as f64).log10())
}

fn pct(part: u64, whole: u64) -> f64 {
    if whole == 0 {
        0.0
    } else {
        100.0 * part as f64 / whole as f64
    }
}

fn summarize(stats: &Stats, total_reads: u64, genome_size: Option<u64>) -> Summary {
    let mut lengths: Vec<u32> = stats.lengths.clone();
    lengths.sort_unstable();
    let n = lengths.len();
    let (n50, n90) = nx_values(&lengths, stats.bases);
    let read_q_total = stats.reads_with_q;
    let reads_at_least = |q: usize| stats.reads_at_q[q / 10 - 1];
    let q_summary = |value: f64| (read_q_total > 0).then_some(value);
    let gc_bases: u64 = stats.bases_from_end.chunks(5).map(|c| c[1] + c[2]).sum();
    let acgt_bases: u64 = stats
        .bases_from_end
        .chunks(5)
        .map(|c| c[..4].iter().sum::<u64>())
        .sum();
    let filtered = stats.pass.reads + stats.fail.reads;
    let mapping = &stats.mapping;
    let mapped_total = mapping.mapped + mapping.unmapped;
    let five_mc = stats.mods.get("C+m");

    let coverage = |bases: u64| genome_size.map(|g| bases as f64 / g as f64);
    Summary {
        coverage: coverage(stats.bases),
        pass_coverage: coverage(stats.pass.bases),
        reads: stats.reads,
        bases: stats.bases,
        share_pct: pct(stats.reads, total_reads),
        mean_length: if n == 0 {
            0.0
        } else {
            stats.bases as f64 / n as f64
        },
        median_length: match n {
            0 => 0.0,
            _ if n % 2 == 1 => f64::from(lengths[n / 2]),
            _ => (f64::from(lengths[n / 2 - 1]) + f64::from(lengths[n / 2])) / 2.0,
        },
        n50,
        n90,
        shortest: lengths.first().map_or(0, |&l| u64::from(l)),
        longest: lengths.last().map_or(0, |&l| u64::from(l)),
        median_read_q: q_summary(histogram_median(&stats.read_q) / 10.0),
        mean_read_q: q_summary(
            stats
                .read_q
                .iter()
                .enumerate()
                .map(|(i, &c)| i as f64 * c as f64)
                .sum::<f64>()
                / 10.0
                / read_q_total.max(1) as f64,
        ),
        reads_q10_pct: q_summary(pct(reads_at_least(10), read_q_total)),
        reads_q20_pct: q_summary(pct(reads_at_least(20), read_q_total)),
        reads_q30_pct: q_summary(pct(reads_at_least(30), read_q_total)),
        bases_q20: (read_q_total > 0).then_some(stats.bases_at_q20),
        base_q: error_averaged_q(&stats.base_quals),
        gc_pct: (acgt_bases > 0).then(|| pct(gc_bases, acgt_bases)),
        pass_reads: stats.pass.reads,
        pass_bases: stats.pass.bases,
        pass_pct: (filtered > 0).then(|| pct(stats.pass.reads, filtered)),
        mapped_pct: (mapped_total > 0).then(|| pct(mapping.mapped, mapped_total)),
        median_identity_pct: (mapping.identity.iter().sum::<u64>() > 0)
            .then(|| percentile(&mapping.identity, 0.5) / 10.0),
        duplicate_names: stats.duplicate_names,
        mod_reads_pct: (stats.reads_with_mods > 0).then(|| pct(stats.reads_with_mods, stats.reads)),
        methylated_pct: five_mc
            .filter(|m| m.cpg_calls > 0)
            .map(|m| pct(m.cpg_modified, m.cpg_calls)),
    }
}

/// The standard median of the values a histogram counts (its bin indexes): the middle value,
/// or for an even count the mean of the two middle ones.
fn histogram_median(hist: &[u64]) -> f64 {
    let total: u64 = hist.iter().sum();
    if total == 0 {
        return 0.0;
    }
    // The value at 1-based rank `rank`.
    let nth = |rank: u64| {
        let mut seen = 0;
        for (value, &count) in hist.iter().enumerate() {
            seen += count;
            if seen >= rank {
                return value as f64;
            }
        }
        (hist.len() - 1) as f64
    };
    if total % 2 == 1 {
        nth(total / 2 + 1)
    } else {
        (nth(total / 2) + nth(total / 2 + 1)) / 2.0
    }
}

/// N50 and N90 of sorted lengths.
fn nx_values(sorted: &[u32], total: u64) -> (u64, u64) {
    let mut n50 = 0;
    let mut n90 = 0;
    let mut cumulative = 0u64;
    for &len in sorted.iter().rev() {
        cumulative += u64::from(len);
        if n50 == 0 && cumulative * 2 >= total {
            n50 = u64::from(len);
        }
        if cumulative * 10 >= total * 9 {
            n90 = u64::from(len);
            break;
        }
    }
    (n50, n90)
}

/// The first bin at which a histogram's cumulative count reaches `fraction` of its total,
/// without interpolation: read Q histograms are in tenths, so medians are to 0.1.
pub fn percentile(hist: &[u64], fraction: f64) -> f64 {
    let total: u64 = hist.iter().sum();
    if total == 0 {
        return 0.0;
    }
    let target = fraction * total as f64;
    let mut cumulative = 0.0;
    for (i, &count) in hist.iter().enumerate() {
        let next = cumulative + count as f64;
        if next >= target && count > 0 {
            return i as f64;
        }
        cumulative = next;
    }
    (hist.len() - 1) as f64
}

impl Data {
    fn build(stats: &Stats, context: &Context<'_>) -> Self {
        let mut data = Data {
            length_bins: length_bins(stats),
            read_q: read_q_bins(stats),
            length_q: stats
                .length_q
                .iter()
                .enumerate()
                .filter(|&(_, &c)| c > 0)
                .map(|(i, &c)| (i / READ_Q_LEVELS, i % READ_Q_LEVELS, c))
                .collect(),
            quality_from_start: position_quality(&stats.quals_from_start, stats.reads_with_quals),
            quality_from_end: position_quality(&stats.quals_from_end, stats.reads_with_quals),
            content_from_start: position_content(&stats.bases_from_start, stats.reads_with_seq),
            content_from_end: position_content(&stats.bases_from_end, stats.reads_with_seq),
            base_quals: trim_trailing_zeros(&stats.base_quals),
            gc: stats.gc.to_vec(),
            speed: trim_trailing_zeros(&stats.speed),
            ..Data::default()
        };
        if stats.reads_with_seq > 0 {
            let (expected, deviation) = gc_expected(&stats.gc);
            data.gc_expected = expected;
            data.gc_deviation_pct = Some(deviation);
            // Peaks are judged on reads of 500 bases or more when there are enough of them,
            // so adapter fragments do not pass for a second organism.
            let long: u64 = stats.gc_long.iter().sum();
            data.gc_second_peak = if long >= 100 {
                second_gc_peak(&stats.gc_long).map(|p| GcPeak {
                    min_length: stats::GC_PEAK_MIN_LENGTH,
                    ..p
                })
            } else {
                second_gc_peak(&stats.gc)
            };
            data.adapters = adapter_summaries(stats, context.adapters);
            data.adapter_at_ends_pct = pct(stats.reads_adapter_at_end, stats.reads_with_seq);
            data.adapter_internal_pct = pct(stats.reads_adapter_internal, stats.reads_with_seq);
            data.poly_a_tail_pct = Some(pct(stats.poly_a_tail, stats.reads_with_seq));
            data.poly_t_head_pct = Some(pct(stats.poly_t_head, stats.reads_with_seq));
        }
        if stats.sample.len() > 0 {
            data.duplication = Some(duplication(stats.sample.ends()));
            data.overrepresented = overrepresented(
                stats,
                context.adapter_index,
                listing_pct(context.thresholds),
            );
        }
        let mut channels: Vec<ChannelYield> = stats
            .channels
            .iter()
            .map(|(&channel, y)| ChannelYield {
                channel,
                reads: y.reads,
                bases: y.bases,
            })
            .collect();
        channels.sort_by_key(|c| c.channel);
        data.channels = channels;
        let mut reasons: Vec<(String, u64)> = stats
            .end_reasons
            .iter()
            .map(|(r, &c)| (r.to_string(), c))
            .collect();
        reasons.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
        data.end_reasons = reasons;
        data.pacbio = pacbio(stats);
        data.mods = mod_summaries(stats);
        data.mapping = mapping(stats);
        data.filter_grid = filter_grid(&stats.lengths, &stats.read_q_tenths);
        if let Some(references) = context.screen {
            data.reads_screened = stats.reads_screened;
            data.contaminants = references
                .iter()
                .zip(&stats.screen)
                .zip(&stats.screen_testable)
                .map(|((reference, y), &testable_called)| {
                    // A sparse reference calls only long reads. When they are a good part of
                    // the run, its share is of them; otherwise it is of all reads, a lower
                    // bound.
                    let n = stats.sparse_testable;
                    let of_testable = reference.sparse
                        && n >= MIN_TESTABLE_READS
                        && n as f64 >= MIN_TESTABLE_FRACTION * stats.reads_screened as f64;
                    ScreenHit {
                        reference: reference.clone(),
                        reads: y.reads,
                        bases: y.bases,
                        reads_pct: if of_testable {
                            pct(testable_called, n)
                        } else {
                            pct(y.reads, stats.reads_screened)
                        },
                        bases_pct: pct(y.bases, stats.bases),
                        testable: of_testable.then_some(n),
                        lower_bound: reference.sparse && !of_testable,
                    }
                })
                .collect();
        }
        data.barcode_scan = barcode_scan(stats);
        data.ont = ont_extras(stats);
        data.cdna = cdna(stats);
        data.poly_a = trim_trailing_zeros(&stats.ont.poly_a);
        data.kinnex = if stats.kinnex.iter().skip(1).any(|&c| c > 0) {
            trim_trailing_zeros(&stats.kinnex)
        } else {
            Vec::new()
        };
        data
    }
}

fn length_bins(stats: &Stats) -> Vec<LengthBin> {
    let mut reads = vec![0u64; bins::LENGTH_BINS];
    let mut bases = vec![0u64; bins::LENGTH_BINS];
    for &len in &stats.lengths {
        let bin = bins::length_bin(u64::from(len));
        reads[bin] += 1;
        bases[bin] += u64::from(len);
    }
    let Some(first) = reads.iter().position(|&c| c > 0) else {
        return Vec::new();
    };
    let last = reads.iter().rposition(|&c| c > 0).unwrap_or(first);
    (first..=last)
        .map(|bin| LengthBin {
            start: bins::length_bin_start(bin),
            end: bins::length_bin_start(bin + 1),
            reads: reads[bin],
            bases: bases[bin],
        })
        .collect()
}

fn read_q_bins(stats: &Stats) -> Vec<QBin> {
    let Some(first) = stats.read_q.iter().position(|&c| c > 0) else {
        return Vec::new();
    };
    let last = stats.read_q.iter().rposition(|&c| c > 0).unwrap_or(first);
    (first..=last)
        .map(|i| QBin {
            q: i as f64 / 10.0,
            reads: stats.read_q[i],
            bases: stats.read_q_bases[i],
        })
        .collect()
}

fn position_quality(counts: &[u64], reads: u64) -> Vec<PositionQuality> {
    let mut out = Vec::new();
    for (bin, row) in counts.as_chunks::<QUAL_LEVELS>().0.iter().enumerate() {
        let bases: u64 = row.iter().sum();
        let covering = bases as f64 / POSITIONS.width(bin) as f64;
        if bases == 0 || covering < (reads as f64 * MIN_POSITION_COVERAGE).max(1.0) {
            break;
        }
        let mean = row
            .iter()
            .enumerate()
            .map(|(q, &c)| q as f64 * c as f64)
            .sum::<f64>()
            / bases as f64;
        out.push(PositionQuality {
            bin,
            reads: covering,
            mean,
            p10: percentile(row, 0.10),
            p25: percentile(row, 0.25),
            median: percentile(row, 0.50),
            p75: percentile(row, 0.75),
            p90: percentile(row, 0.90),
        });
    }
    out
}

fn position_content(counts: &[u64], reads: u64) -> Vec<PositionContent> {
    let mut out = Vec::new();
    for (bin, row) in counts.as_chunks::<5>().0.iter().enumerate() {
        let bases: u64 = row.iter().sum();
        let covering = bases as f64 / POSITIONS.width(bin) as f64;
        if bases == 0 || covering < (reads as f64 * MIN_POSITION_COVERAGE).max(1.0) {
            break;
        }
        let p = |i: usize| pct(row[i], bases);
        out.push(PositionContent {
            bin,
            bases,
            a: p(0),
            c: p(1),
            g: p(2),
            t: p(3),
            n: p(usize::from(BASE_N)),
        });
    }
    out
}

fn trim_trailing_zeros(values: &[u64]) -> Vec<u64> {
    let end = values.iter().rposition(|&v| v > 0).map_or(0, |i| i + 1);
    values[..end].to_vec()
}

/// A normal curve with the observed GC mean and deviation, scaled to the read count, and
/// the percent of reads by which the observed distribution departs from it.
fn gc_expected(gc: &[u64; 101]) -> (Vec<f64>, f64) {
    let total: u64 = gc.iter().sum();
    if total == 0 {
        return (vec![0.0; 101], 0.0);
    }
    let n = total as f64;
    let mean = gc
        .iter()
        .enumerate()
        .map(|(i, &c)| i as f64 * c as f64)
        .sum::<f64>()
        / n;
    let variance = gc
        .iter()
        .enumerate()
        .map(|(i, &c)| (i as f64 - mean).powi(2) * c as f64)
        .sum::<f64>()
        / n;
    let sd = variance.sqrt().max(0.5);
    let density: Vec<f64> = (0..=100)
        .map(|i| (-0.5 * ((i as f64 - mean) / sd).powi(2)).exp())
        .collect();
    let scale = n / density.iter().sum::<f64>();
    let expected: Vec<f64> = density.iter().map(|d| d * scale).collect();
    let deviation = gc
        .iter()
        .zip(&expected)
        .map(|(&o, &e)| (o as f64 - e).abs())
        .sum::<f64>()
        / n
        * 100.0;
    (expected, deviation)
}

/// A second peak of the GC distribution: a local maximum at least 8 points from the main
/// peak, a tenth of its height or more, and separated from it by a valley.
fn second_gc_peak(gc: &[u64; 101]) -> Option<GcPeak> {
    let total: u64 = gc.iter().sum();
    if total < 100 {
        return None;
    }
    let smooth: Vec<f64> = (0..=100)
        .map(|i: usize| {
            let window = &gc[i.saturating_sub(2)..=(i + 2).min(100)];
            window.iter().sum::<u64>() as f64 / window.len() as f64
        })
        .collect();
    let main = (0..=100).max_by(|&a, &b| smooth[a].total_cmp(&smooth[b]))?;
    let is_local_max =
        |i: usize| (i.saturating_sub(3)..=(i + 3).min(100)).all(|j| smooth[j] <= smooth[i]);
    (0..=100)
        .filter(|&i: &usize| {
            i.abs_diff(main) >= 8 && smooth[i] >= 0.1 * smooth[main] && is_local_max(i)
        })
        .filter(|&i| {
            let (lo, hi) = (i.min(main), i.max(main));
            let valley = smooth[lo..=hi]
                .iter()
                .copied()
                .fold(f64::INFINITY, f64::min);
            valley < 0.7 * smooth[i]
        })
        .map(|i| {
            let share: u64 = gc[i.saturating_sub(5)..=(i + 5).min(100)].iter().sum();
            GcPeak {
                gc: i,
                main_gc: main,
                share_pct: pct(share, total),
                min_length: 0,
            }
        })
        .max_by(|a, b| a.share_pct.total_cmp(&b.share_pct))
}

/// FastQC's duplication levels, over reads keyed by both ends.
fn duplication(ends: impl Iterator<Item = Ends>) -> Duplication {
    let mut counts: HashMap<Ends, u64> = HashMap::new();
    let mut sampled = 0usize;
    for e in ends {
        *counts.entry(e).or_default() += 1;
        sampled += 1;
    }
    const LEVELS: [(&str, u64, u64); 16] = [
        ("1", 1, 1),
        ("2", 2, 2),
        ("3", 3, 3),
        ("4", 4, 4),
        ("5", 5, 5),
        ("6", 6, 6),
        ("7", 7, 7),
        ("8", 8, 8),
        ("9", 9, 9),
        (">10", 10, 49),
        (">50", 50, 99),
        (">100", 100, 499),
        (">500", 500, 999),
        (">1k", 1_000, 4_999),
        (">5k", 5_000, 9_999),
        (">10k", 10_000, u64::MAX),
    ];
    let mut levels = vec![0u64; LEVELS.len()];
    for &count in counts.values() {
        if let Some(i) = LEVELS
            .iter()
            .position(|&(_, lo, hi)| (lo..=hi).contains(&count))
        {
            levels[i] += count;
        }
    }
    let total = sampled as u64;
    Duplication {
        sampled,
        distinct_pct: pct(counts.len() as u64, total),
        levels: LEVELS
            .iter()
            .zip(levels)
            .map(|(&(name, _, _), reads)| (name.to_owned(), pct(reads, total)))
            .collect(),
    }
}

/// The share above which end sequences are listed: FastQC's 0.1%, or the warn level
/// when that is lower, so every sequence the check grades is listed.
pub fn listing_pct(thresholds: &Thresholds) -> f64 {
    (OVERREPRESENTED_FRACTION * 100.0).min(thresholds.overrepresented_warn)
}

fn overrepresented(stats: &Stats, adapters: &Adapters, listing_pct: f64) -> Vec<Overrepresented> {
    let sampled = stats.sample.len() as u64;
    // Above the listing share, and at least twice: one read is not a repeat.
    let listed = |count: u64| count >= 2 && count as f64 * 100.0 > sampled as f64 * listing_pct;
    let mut out = Vec::new();
    for end in [ReadEnd::Start, ReadEnd::End] {
        let mut counts: HashMap<u128, u64> = HashMap::new();
        for ends in stats.sample.ends() {
            let key = match end {
                ReadEnd::Start => ends.start,
                ReadEnd::End => ends.end,
            };
            *counts.entry(key).or_default() += 1;
        }
        let mut frequent: Vec<(u128, u64)> =
            counts.into_iter().filter(|&(_, c)| listed(c)).collect();
        frequent.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
        for (packed, count) in frequent.into_iter().take(OVERREPRESENTED_LISTED) {
            let sequence = sampler::unpack(packed);
            out.push(Overrepresented {
                end,
                source: possible_source(&sequence, adapters),
                sequence,
                count,
                pct: pct(count, sampled),
            });
        }
    }
    out
}

fn possible_source(sequence: &str, adapters: &Adapters) -> String {
    let codes: Vec<u8> = sequence
        .bytes()
        .map(|b| crate::read::ASCII_TO_CODE[usize::from(b)])
        .collect();
    let mut hits = Vec::new();
    adapters.scan(
        &codes,
        &mut crate::adapters::ScanBuffer::default(),
        &mut hits,
    );
    if let Some(hit) = hits.first() {
        return adapters.names()[hit.adapter].clone();
    }
    let max_run = sequence
        .as_bytes()
        .chunk_by(|a, b| a == b)
        .map(<[u8]>::len)
        .max()
        .unwrap_or(0);
    if max_run >= 20 {
        "Homopolymer".to_owned()
    } else if let Some(period) = tandem_period(sequence.as_bytes()) {
        format!("Simple repeat (period {period})")
    } else {
        "No hit".to_owned()
    }
}

/// The shortest period, up to 6, at which a sequence repeats itself almost exactly.
fn tandem_period(seq: &[u8]) -> Option<usize> {
    (1..=6).find(|&period| {
        let mismatches = seq
            .iter()
            .zip(&seq[period..])
            .filter(|(a, b)| a != b)
            .count();
        mismatches * 10 <= seq.len() - period
    })
}

fn adapter_summaries(stats: &Stats, names: &[String]) -> Vec<AdapterSummary> {
    let reads = stats.reads_with_seq;
    let cumulative = |counts: &[u64]| -> Vec<f64> {
        let mut running = 0;
        counts
            .iter()
            .map(|&c| {
                running += c;
                pct(running, reads)
            })
            .collect()
    };
    stats
        .adapters
        .iter()
        .zip(names)
        .filter(|(counts, _)| {
            counts.reads_at_start + counts.reads_at_end + counts.reads_internal > 0
        })
        .map(|(counts, name)| AdapterSummary {
            name: name.clone(),
            start_pct: pct(counts.reads_at_start, reads),
            end_pct: pct(counts.reads_at_end, reads),
            internal_pct: pct(counts.reads_internal, reads),
            cumulative_from_start: cumulative(&counts.nearest_start),
            cumulative_from_end: cumulative(&counts.nearest_end),
        })
        .collect()
}

/// Run time in slices of whole minutes, at most about 300 slices.
/// Yield over time in at most 300 slices, from `origin` (the run's first minute) or else
/// the group's own first read.
fn time_points(stats: &Stats, origin: Option<i64>) -> Vec<TimePoint> {
    let Some(&own_first) = stats.time.keys().min() else {
        return Vec::new();
    };
    let first = origin.map_or(own_first, |o| o.min(own_first));
    let last = stats.time.keys().max().copied().unwrap_or(first);
    let span = last - first + 1;
    let width = (span + 299) / 300;
    let slices = usize::try_from((last - first) / width + 1).unwrap_or(1);
    let mut merged = vec![stats::TimeBin::default(); slices];
    for (&minute, bin) in &stats.time {
        let slot = &mut merged[usize::try_from((minute - first) / width).unwrap_or(0)];
        slot.reads += bin.reads;
        slot.bases += bin.bases;
        slot.pass_reads += bin.pass_reads;
        slot.q_tenths += bin.q_tenths;
        slot.q_reads += bin.q_reads;
    }
    let mut cumulative = 0;
    merged
        .iter()
        .enumerate()
        .map(|(i, bin)| {
            cumulative += bin.bases;
            TimePoint {
                hours: (i as i64 * width) as f64 / 60.0,
                reads: bin.reads,
                bases: bin.bases,
                cumulative_bases: cumulative,
                mean_q: (bin.q_reads > 0).then(|| bin.q_tenths as f64 / 10.0 / bin.q_reads as f64),
                pass_pct: (bin.reads > 0).then(|| pct(bin.pass_reads, bin.reads)),
                mean_length: (bin.reads > 0).then(|| bin.bases as f64 / bin.reads as f64),
            }
        })
        .collect()
}

fn pacbio(stats: &Stats) -> Option<PacBio> {
    let with_passes: u64 = stats.passes.iter().sum();
    let with_coverage: u64 = stats.effective_coverage.iter().sum();
    let with_score: u64 = stats.barcode_score.iter().sum();
    if with_passes + with_coverage + with_score + stats.with_kinetics == 0 {
        return None;
    }
    Some(PacBio {
        passes: trim_trailing_zeros(&stats.passes),
        effective_coverage: trim_trailing_zeros(&stats.effective_coverage),
        barcode_score: trim_trailing_zeros(&stats.barcode_score),
        median_passes: (with_passes > 0).then(|| percentile(&stats.passes, 0.5) as u64),
        kinetics_pct: pct(stats.with_kinetics, stats.reads),
        passes_q: stats
            .pacbio
            .passes_q
            .iter()
            .enumerate()
            .filter(|&(_, &c)| c > 0)
            .map(|(i, &c)| (i / READ_Q_LEVELS, i % READ_Q_LEVELS, c))
            .collect(),
        reads_per_zmw: {
            let mut counts = vec![0u64; 5];
            for &n in stats.pacbio.zmw_reads.values() {
                counts[(n as usize).clamp(1, 5) - 1] += 1;
            }
            counts
        },
        zmws: stats.pacbio.zmw_reads.len() as u64,
    })
}

/// Names of the modification codes of the SAM specification and dorado's ChEBI codes. In
/// RNA, `C+m` and `A+a` are m5C and m6A.
fn mod_label(code: &str, molecule: Molecule) -> String {
    let rna = molecule == Molecule::Rna;
    let name = match code.get(2..) {
        Some("m") if rna => "m5C",
        Some("a") if rna => "m6A",
        Some("m") => "5mC",
        Some("h") => "5hmC",
        Some("f") => "5fC",
        Some("c") => "5caC",
        Some("a") => "6mA",
        Some("o") => "8oxoG",
        Some("e") => "5hmU",
        Some("21839") => "4mC",
        Some("17596") => "inosine",
        Some("17802") => "pseudouridine",
        Some("69426") => "2'-O-methyladenosine (Am)",
        Some("19228") => "2'-O-methylcytidine (Cm)",
        Some("19229") => "2'-O-methylguanosine (Gm)",
        Some("19227") => "2'-O-methyluridine (Um)",
        _ => return code.to_owned(),
    };
    format!("{name} ({code})")
}

fn mod_summaries(stats: &Stats) -> Vec<ModSummary> {
    let mut out: Vec<ModSummary> = stats
        .mods
        .iter()
        .map(|(code, m)| ModSummary {
            code: code.clone(),
            label: mod_label(code, Molecule::Dna),
            calls: m.calls,
            modified_pct: pct(m.modified, m.calls),
            cpg_calls: m.cpg_calls,
            cpg_modified_pct: (m.cpg_calls > 0).then(|| pct(m.cpg_modified, m.cpg_calls)),
            calls_per_read: m.calls as f64 / stats.reads_with_mods.max(1) as f64,
            probs: m.probs.to_vec(),
        })
        .collect();
    out.sort_by(|a, b| b.calls.cmp(&a.calls).then(a.code.cmp(&b.code)));
    out
}

fn mapping(stats: &Stats) -> Option<MappingSummary> {
    let m = &stats.mapping;
    if m.mapped + m.unmapped == 0 {
        return None;
    }
    // Plot identity from the lowest occupied bin, or from 70% if that is lower.
    let first = |h: &[u64]| h.iter().position(|&c| c > 0).unwrap_or(700);
    let floor_bin = first(&m.identity)
        .min(first(&m.gap_compressed_identity))
        .min(700);
    let with_identity: u64 = m.identity.iter().sum();
    Some(MappingSummary {
        mapped: m.mapped,
        unmapped: m.unmapped,
        secondary: m.secondary,
        supplementary: m.supplementary,
        mapped_pct: pct(m.mapped, m.mapped + m.unmapped),
        split_pct: pct(m.split, m.mapped),
        mapq: m.mapq.to_vec(),
        identity: m.identity[floor_bin..].to_vec(),
        gap_compressed_identity: m.gap_compressed_identity[floor_bin..].to_vec(),
        identity_floor: floor_bin as f64 / 10.0,
        median_identity_pct: (with_identity > 0).then(|| percentile(&m.identity, 0.5) / 10.0),
        median_gap_compressed_pct: (with_identity > 0)
            .then(|| percentile(&m.gap_compressed_identity, 0.5) / 10.0),
        aligned_fraction: m.aligned_fraction.to_vec(),
    })
}

/// Minimum lengths and read Qs of the filter grid.
const GRID_LENGTHS: [u64; 8] = [0, 1_000, 5_000, 10_000, 15_000, 20_000, 50_000, 100_000];
const GRID_QS: [u16; 7] = [0, 70, 100, 150, 200, 250, 300];

/// Reads and bases passing each pair of length and read Q minimums, counted exactly: each
/// read lands in the cell of the highest minimums it meets, then cells are summed over
/// everything at or above them.
fn filter_grid(lengths: &[u32], qs: &[u16]) -> Option<FilterGrid> {
    if lengths.is_empty() {
        return None;
    }
    let has_q = qs.iter().any(|&q| q != stats::NO_Q);
    let q_levels = if has_q { GRID_QS.len() } else { 1 };
    let (rows, cols) = (GRID_LENGTHS.len(), q_levels);
    let mut reads = vec![vec![0u64; cols]; rows];
    let mut bases = vec![vec![0u64; cols]; rows];
    for (&len, &q) in lengths.iter().zip(qs) {
        let row = GRID_LENGTHS
            .iter()
            .rposition(|&l| u64::from(len) >= l)
            .unwrap_or(0);
        let col = if has_q && q != stats::NO_Q {
            GRID_QS[..cols].iter().rposition(|&g| q >= g).unwrap_or(0)
        } else {
            0
        };
        reads[row][col] += 1;
        bases[row][col] += u64::from(len);
    }
    for grid in [&mut reads, &mut bases] {
        for r in (0..rows).rev() {
            for c in (0..cols).rev() {
                let below = if r + 1 < rows { grid[r + 1][c] } else { 0 };
                let right = if c + 1 < cols { grid[r][c + 1] } else { 0 };
                let diagonal = if r + 1 < rows && c + 1 < cols {
                    grid[r + 1][c + 1]
                } else {
                    0
                };
                grid[r][c] += below + right - diagonal;
            }
        }
    }
    Some(FilterGrid {
        min_lengths: GRID_LENGTHS.to_vec(),
        min_qs: GRID_QS[..cols]
            .iter()
            .map(|&q| f64::from(q) / 10.0)
            .collect(),
        reads,
        bases,
    })
}

fn barcode_scan(stats: &Stats) -> Option<BarcodeScan> {
    let b = &stats.barcodes;
    if b.reads_scanned == 0 {
        return None;
    }
    let mut foreign: Vec<(String, u64, f64)> = b
        .foreign
        .iter()
        .map(|(name, &n)| (name.clone(), n, pct(n, b.reads_scanned)))
        .collect();
    foreign.sort_by(|x, y| y.1.cmp(&x.1).then(x.0.cmp(&y.0)));
    foreign.truncate(10);
    Some(BarcodeScan {
        scanned: b.reads_scanned,
        own_at_end_pct: pct(b.own_at_end, b.reads_scanned),
        own_internal_pct: pct(b.own_internal, b.reads_scanned),
        foreign_pct: pct(b.with_foreign, b.reads_scanned),
        foreign,
    })
}

/// A sparse histogram as sorted bin starts and counts.
fn sparse_histogram(map: &HashMap<i64, u64>, width: f64) -> Vec<(f64, u64)> {
    let mut bins: Vec<(f64, u64)> = map.iter().map(|(&b, &c)| (b as f64 * width, c)).collect();
    bins.sort_by(|a, b| a.0.total_cmp(&b.0));
    bins
}

fn ont_extras(stats: &Stats) -> Option<OntExtras> {
    let o = &stats.ont;
    let classified: u64 = o.duplex.iter().map(|y| y.reads).sum();
    let any = classified > 0
        || o.split_reads.reads > 0
        || !o.mux.is_empty()
        || !o.shift.is_empty()
        || !o.samples_per_base.is_empty();
    if !any {
        return None;
    }
    let duplex = (classified > 0).then(|| {
        let names = ["simplex", "duplex", "simplex with duplex offspring"];
        let classes = names
            .iter()
            .enumerate()
            .map(|(i, name)| {
                let hist = &o.duplex_q[i * stats::READ_Q_BINS..(i + 1) * stats::READ_Q_BINS];
                let with_q: u64 = hist.iter().sum();
                (
                    (*name).to_owned(),
                    o.duplex[i].reads,
                    o.duplex[i].bases,
                    (with_q > 0).then(|| percentile(hist, 0.5) / 10.0),
                )
            })
            .collect();
        let duplex_bases = o.duplex[1].bases;
        let parent_bases = o.duplex[2].bases;
        let simplex_bases = o.duplex[0].bases + parent_bases;
        DuplexSummary {
            classes,
            duplex_rate_pct: pct(2 * duplex_bases, simplex_bases),
        }
    });
    let mut mux: Vec<(u8, u64, u64)> = o.mux.iter().map(|(&m, y)| (m, y.reads, y.bases)).collect();
    mux.sort_unstable();
    Some(OntExtras {
        duplex,
        split_reads: o.split_reads.reads,
        split_pct: pct(o.split_reads.reads, stats.reads),
        mux,
        shift: sparse_histogram(&o.shift, stats::SHIFT_BIN),
        scale: sparse_histogram(&o.scale, stats::SCALE_BIN),
        samples_per_base: sparse_histogram(&o.samples_per_base, stats::SAMPLES_PER_BASE_BIN),
    })
}

/// cDNA structure, shown when at least 1% of reads carry a cDNA primer at an end.
fn cdna(stats: &Stats) -> Option<CdnaSummary> {
    let c = &stats.cdna;
    let reads = stats.reads_with_seq;
    let with_primer = c.full_length_forward + c.full_length_reverse + c.ssp_only + c.vnp_only;
    (reads > 0 && with_primer * 100 >= reads).then(|| CdnaSummary {
        full_length_pct: pct(c.full_length_forward + c.full_length_reverse, reads),
        forward_pct: pct(c.full_length_forward, reads),
        reverse_pct: pct(c.full_length_reverse, reads),
        ssp_only_pct: pct(c.ssp_only, reads),
        vnp_only_pct: pct(c.vnp_only, reads),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lists_every_end_sequence_above_the_listing_share() {
        assert!(listing_pct(&Thresholds::default()) <= 0.1);
        let low = Thresholds {
            overrepresented_warn: 0.05,
            ..Thresholds::default()
        };
        assert!((listing_pct(&low) - 0.05).abs() < 1e-12);
    }

    #[test]
    fn takes_the_middle_of_a_histogram() {
        assert!((histogram_median(&[0, 2, 0, 1]) - 1.0).abs() < 1e-9);
        assert!((histogram_median(&[0, 1, 0, 1]) - 2.0).abs() < 1e-9);
        assert!((histogram_median(&[0, 1, 1, 1, 1]) - 2.5).abs() < 1e-9);
        assert!(histogram_median(&[0, 0]).abs() < 1e-9);
    }

    #[test]
    fn takes_the_middle_of_an_even_number_of_lengths() {
        let mut stats = Stats::new(0);
        stats.lengths = vec![100, 400, 200, 300];
        stats.bases = 1_000;
        assert!((summarize(&stats, 4, None).median_length - 250.0).abs() < 1e-9);
        stats.lengths.push(500);
        stats.bases += 500;
        assert!((summarize(&stats, 5, None).median_length - 300.0).abs() < 1e-9);
    }

    #[test]
    fn base_q_averages_error_probabilities() {
        // Half Q10 (error 0.1) and half Q50 (0.00001): mean error 0.050005, about Q13,
        // where averaging the Phred values would give Q30.
        let mut counts = vec![0; 51];
        counts[10] = 1;
        counts[50] = 1;
        let q = super::error_averaged_q(&counts).expect("bases");
        assert!((q - 13.01).abs() < 0.01, "{q}");
        assert_eq!(super::error_averaged_q(&[0; 51]), None);
    }

    #[test]
    fn filter_grid_counts_reads_meeting_both_minimums() {
        let lengths = [500, 1_500, 12_000, 30_000];
        let qs = [80, 120, 210, stats::NO_Q];
        let grid = filter_grid(&lengths, &qs).expect("reads");
        assert_eq!(grid.reads[0][0], 4, "no filter keeps everything");
        // At least 1 kb and Q10: the 1.5 kb Q12 and 12 kb Q21 reads.
        assert_eq!(grid.reads[1][2], 2);
        assert_eq!(grid.bases[1][2], 13_500);
        // At least 10 kb and Q20: only the 12 kb read; the 30 kb read has no Q.
        assert_eq!(grid.reads[3][4], 1);
        assert_eq!(grid.reads[3][0], 2);
    }

    #[test]
    fn labels_rna_and_dna_modifications() {
        assert_eq!(mod_label("C+m", Molecule::Dna), "5mC (C+m)");
        assert_eq!(mod_label("C+m", Molecule::Rna), "m5C (C+m)");
        assert_eq!(mod_label("A+a", Molecule::Rna), "m6A (A+a)");
        assert_eq!(
            mod_label("A+69426", Molecule::Rna),
            "2'-O-methyladenosine (Am) (A+69426)"
        );
        assert_eq!(
            mod_label("T+19227", Molecule::Rna),
            "2'-O-methyluridine (Um) (T+19227)"
        );
    }

    #[test]
    fn names_symmetric_pacbio_pairs_from_single_barcodes() {
        let samples: HashMap<String, String> =
            [("bc2052".to_owned(), "B. subtilis".to_owned())].into();
        assert_eq!(
            sample_name(&samples, "bc2052--bc2052").as_deref(),
            Some("B. subtilis")
        );
        assert_eq!(
            sample_name(&samples, "bc2052").as_deref(),
            Some("B. subtilis")
        );
        assert_eq!(sample_name(&samples, "bc2052--bc2053"), None);
    }

    #[test]
    fn computes_nx() {
        let lengths = [1, 2, 3, 4, 10];
        // From the longest: 10 reaches half of 20; 10+4+3+2 = 19 reaches 90%.
        assert_eq!(nx_values(&lengths, 20), (10, 2));
    }

    #[test]
    fn percentiles_of_histograms() {
        let hist = [0, 5, 0, 5];
        assert_eq!(percentile(&hist, 0.5), 1.0);
        assert_eq!(percentile(&hist, 0.9), 3.0);
        assert_eq!(percentile(&[0, 0], 0.5), 0.0);
    }

    #[test]
    fn gc_of_a_normal_sample_barely_deviates() {
        let mut gc = [0u64; 101];
        for (i, slot) in gc.iter_mut().enumerate() {
            *slot = (1000.0 * (-0.5 * ((i as f64 - 41.0) / 4.0).powi(2)).exp()) as u64;
        }
        let (_, deviation) = gc_expected(&gc);
        assert!(deviation < 2.0, "deviation {deviation}");
    }

    #[test]
    fn finds_a_contaminating_gc_peak() {
        let curve = |center: f64, height: f64| -> [u64; 101] {
            let mut gc = [0u64; 101];
            for (i, slot) in gc.iter_mut().enumerate() {
                *slot = (height * (-0.5 * ((i as f64 - center) / 3.0).powi(2)).exp()) as u64;
            }
            gc
        };
        let human = curve(41.0, 1000.0);
        assert!(second_gc_peak(&human).is_none());
        let mut mixed = human;
        for (slot, extra) in mixed.iter_mut().zip(curve(65.0, 300.0)) {
            *slot += extra;
        }
        let peak = second_gc_peak(&mixed).expect("a second peak");
        assert_eq!((peak.gc, peak.main_gc), (65, 41));
        assert!((peak.share_pct - 23.0).abs() < 2.0, "{}", peak.share_pct);
    }

    #[test]
    fn duplication_levels_sum_to_all_reads() {
        let one = Ends { start: 1, end: 1 };
        let two = Ends { start: 2, end: 2 };
        let d = duplication([one, one, one, two].into_iter());
        assert_eq!(d.distinct_pct, 50.0);
        let total: f64 = d.levels.iter().map(|(_, p)| p).sum();
        assert!((total - 100.0).abs() < 1e-9);
        assert_eq!(d.levels[2].1, 75.0);
    }
}
