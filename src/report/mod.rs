//! The report: per-group summaries, plot data and pass/warn/fail checks, and its renderings.

pub mod chart;
pub mod checks;
pub mod html;
pub mod model;
pub mod multiqc;
pub mod runs;
pub mod samples;
pub mod thresholds;
pub mod tsv;

use std::collections::HashMap;

use serde::{Deserialize, Serialize};

use crate::barcode::{NO_BARCODE, UNCLASSIFIED};
use crate::bins::POSITIONS;
use crate::pipeline::{Collected, FileInfo, Options};
use crate::screen::Reference;
use crate::stats::Stats;

pub use model::Group;
pub use runs::RunInfo;
pub use thresholds::Thresholds;

/// Everything the HTML, TSV and JSON outputs show.
#[derive(Debug, Serialize, Deserialize)]
pub struct Report {
    pub tool: String,
    pub version: String,
    pub title: String,
    pub command: String,
    pub inputs: Vec<FileInfo>,
    /// Read groups and programs from BAM headers.
    pub runs: RunInfo,
    /// Position bins of the per-position plots: first base (0-based) and end.
    pub position_bins: Vec<(u64, Option<u64>)>,
    pub adapters: Vec<String>,
    pub screen: Vec<Reference>,
    pub genome_size: Option<u64>,
    /// All reads (or a one-barcode run's barcode) first, then barcodes by read count, then
    /// unclassified and unbarcoded reads by read count.
    pub groups: Vec<Group>,
    /// There are groups besides all reads: barcodes, or runs with --split-runs.
    pub barcoded: bool,
    /// Seconds between refreshes of a report written in watch mode.
    pub refresh_secs: Option<u64>,
    pub thresholds: Thresholds,
}

/// What kind of group of reads a [`Group`] is.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum GroupKind {
    AllReads,
    Barcode,
    Unclassified,
    NoBarcode,
}

/// Report settings that do not affect how reads are counted.
#[derive(Debug)]
pub struct ReportOptions {
    pub title: String,
    pub command: String,
    pub thresholds: Thresholds,
    pub genome_size: Option<u64>,
    /// Sample names by barcode.
    pub samples: HashMap<String, String>,
    pub refresh_secs: Option<u64>,
}

/// The barcode part of a group name (`PAY87954 · barcode01` when runs are split).
pub fn barcode_of(name: &str) -> &str {
    name.rsplit(" · ").next().unwrap_or(name)
}

fn kind_of(name: &str) -> GroupKind {
    match barcode_of(name) {
        UNCLASSIFIED => GroupKind::Unclassified,
        NO_BARCODE => GroupKind::NoBarcode,
        _ => GroupKind::Barcode,
    }
}

/// Builds the report from the pipeline's statistics.
pub fn build(collected: &Collected, options: &Options, report: &ReportOptions) -> Report {
    let groups = &collected.groups;
    let total_reads: u64 = groups.values().map(|s| s.reads).sum();
    // One group, barcoded or not, is the whole run: a separate "All reads" would repeat it.
    let barcoded = groups.len() > 1;

    let mut all = Stats::new(options.adapters.names().len());
    for stats in groups.values() {
        all.merge(stats.clone());
    }

    let mut ordered: Vec<(&String, &Stats)> = groups.iter().collect();
    ordered.sort_by(|(a_name, a), (b_name, b)| {
        let rank = |name: &str| u8::from(kind_of(name) != GroupKind::Barcode);
        rank(a_name)
            .cmp(&rank(b_name))
            .then(b.reads.cmp(&a.reads))
            .then(a_name.cmp(b_name))
    });

    let screen = options.screen.as_ref().map(|s| s.references());
    let mut run_starts: HashMap<String, i64> = HashMap::new();
    for (name, stats) in groups {
        if let Some(&first) = stats.time.keys().min() {
            let start = run_starts.entry(run_of(name).to_owned()).or_insert(first);
            *start = (*start).min(first);
        }
    }
    let context = model::Context {
        adapters: options.adapters.names(),
        adapter_index: &options.adapters,
        min_q: options.min_q,
        total_reads,
        thresholds: &report.thresholds,
        genome_size: report.genome_size,
        samples: &report.samples,
        screen,
        run_starts: &run_starts,
    };
    // A run of one barcode (one MinKNOW barcode folder, one lima file) keeps that barcode's
    // name, and so its sample name; only a run without barcodes is "All reads".
    let single = (!barcoded)
        .then(|| ordered.first())
        .flatten()
        .filter(|(name, _)| kind_of(name) != GroupKind::NoBarcode);
    let mut built = vec![match single {
        Some(&(name, stats)) => model::Group::build(name, kind_of(name), stats, &context),
        None => model::Group::build("All reads", GroupKind::AllReads, &all, &context),
    }];
    if barcoded {
        // With --split-runs, barcodes are balanced, and reads unclassified, within each run.
        let mut run_reads: HashMap<&str, u64> = HashMap::new();
        for &(name, stats) in &ordered {
            *run_reads.entry(run_of(name)).or_default() += stats.reads;
        }
        for &(name, stats) in &ordered {
            let mut group = model::Group::build(name, kind_of(name), stats, &context);
            let run = run_of(name);
            let in_run: Vec<(&String, &Stats)> = ordered
                .iter()
                .filter(|(other, _)| run_of(other) == run)
                .copied()
                .collect();
            checks::barcode_balance(
                &mut group,
                median_barcode_reads(&in_run),
                &report.thresholds,
            );
            let share = 100.0 * stats.reads as f64 / run_reads[run].max(1) as f64;
            checks::unclassified_share(&mut group, share, &report.thresholds);
            built.push(group);
        }
        // Only a demultiplexed run pools samples; runs split from one sample do not.
        let barcodes: std::collections::HashSet<&str> = ordered
            .iter()
            .filter(|(name, _)| kind_of(name) == GroupKind::Barcode)
            .map(|(name, _)| barcode_of(name))
            .collect();
        if barcodes.len() > 1 {
            checks::pools_not_graded_for_mixtures(&mut built);
        }
    }

    Report {
        tool: "lrqc".to_owned(),
        version: env!("CARGO_PKG_VERSION").to_owned(),
        title: report.title.clone(),
        command: report.command.clone(),
        runs: RunInfo::from_inputs(&collected.files),
        inputs: collected.files.clone(),
        position_bins: (0..POSITIONS.len()).map(|b| POSITIONS.range(b)).collect(),
        adapters: options.adapters.names().to_vec(),
        screen: screen.map(<[Reference]>::to_vec).unwrap_or_default(),
        genome_size: report.genome_size,
        groups: built,
        barcoded,
        refresh_secs: report.refresh_secs,
        thresholds: report.thresholds.clone(),
    }
}

/// The run part of a group name (`PAY87954` of `PAY87954 · barcode01`); empty without
/// `--split-runs`.
pub(crate) fn run_of(name: &str) -> &str {
    name.rsplit_once(" · ").map_or("", |(run, _)| run)
}

fn median_barcode_reads(groups: &[(&String, &Stats)]) -> u64 {
    let mut reads: Vec<u64> = groups
        .iter()
        .filter(|(name, _)| kind_of(name) == GroupKind::Barcode)
        .map(|(_, s)| s.reads)
        .collect();
    reads.sort_unstable();
    reads.get(reads.len() / 2).copied().unwrap_or(0)
}

/// Parses a genome size such as `3.1g`, `4.6m`, `12k` or `48502`.
///
/// # Errors
///
/// Fails on anything else.
pub fn parse_genome_size(text: &str) -> anyhow::Result<u64> {
    let text = text.trim().to_ascii_lowercase();
    // `4.6mb`, `4.6mbp`, `48502bp`: the units people write.
    let text = text
        .strip_suffix("bp")
        .or_else(|| text.strip_suffix('b'))
        .unwrap_or(&text)
        .trim_end()
        .to_owned();
    let (number, scale) = match text.strip_suffix(['g', 'm', 'k']) {
        Some(number) => {
            let scale = match text.chars().last() {
                Some('g') => 1e9,
                Some('m') => 1e6,
                _ => 1e3,
            };
            (number, scale)
        }
        None => (text.as_str(), 1.0),
    };
    let value: f64 = number.trim().parse().map_err(|_| {
        anyhow::anyhow!("genome size {text:?} is not a number like 3.1g, 4.6m or 48502")
    })?;
    anyhow::ensure!(
        value.is_finite() && value > 0.0,
        "genome size must be a positive number"
    );
    Ok((value * scale).round() as u64)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_genome_sizes() {
        assert_eq!(parse_genome_size("3.1g").expect("valid"), 3_100_000_000);
        assert_eq!(parse_genome_size("4.6M").expect("valid"), 4_600_000);
        assert_eq!(parse_genome_size("48502").expect("valid"), 48_502);
        assert!(parse_genome_size("big").is_err());
        assert_eq!(parse_genome_size("4.6Mb").expect("valid"), 4_600_000);
        assert_eq!(parse_genome_size("3.1 Gbp").expect("valid"), 3_100_000_000);
        assert_eq!(parse_genome_size("48502bp").expect("valid"), 48_502);
        assert!(parse_genome_size("inf").is_err());
    }

    #[test]
    fn splits_group_names() {
        assert_eq!(barcode_of("PAY87954 · barcode01"), "barcode01");
        assert_eq!(barcode_of("barcode01"), "barcode01");
        assert_eq!(kind_of("PAY87954 · unclassified"), GroupKind::Unclassified);
    }
}
