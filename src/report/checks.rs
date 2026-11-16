//! Pass, warn and fail verdicts, FastQC-style, with long-read thresholds.
//!
//! The levels come from [`Thresholds`]; the reasoning behind each is next to its check.

use serde::{Deserialize, Serialize};

use super::model::{Data, Group, ScreenHit, Summary};
use super::thresholds::Thresholds;
use super::GroupKind;
use crate::read::{Molecule, Platform};
use crate::screen::Kind;
use crate::stats::Stats;

/// Outcome of a check.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
pub enum Status {
    Pass,
    Warn,
    Fail,
    /// Shown without a verdict: what is normal depends on the experiment.
    Info,
    /// The input has no data for this module.
    Absent,
}

impl Status {
    pub fn label(self) -> &'static str {
        match self {
            Status::Pass => "pass",
            Status::Warn => "warn",
            Status::Fail => "fail",
            Status::Info => "info",
            Status::Absent => "n/a",
        }
    }
}

macro_rules! modules {
    ($($id:ident => $key:literal, $title:literal;)*) => {
        /// The report's modules, in display order.
        #[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
        pub enum ModuleId {
            $(#[serde(rename = $key)] $id,)*
        }

        impl ModuleId {
            pub const ALL: &[ModuleId] = &[$(ModuleId::$id,)*];

            pub fn title(self) -> &'static str {
                match self {
                    $(ModuleId::$id => $title,)*
                }
            }

            /// Stable identifier for anchors, JSON and TSV columns.
            pub fn key(self) -> &'static str {
                match self {
                    $(ModuleId::$id => $key,)*
                }
            }
        }
    };
}

modules! {
    BasicStatistics => "basic_statistics", "Basic statistics";
    FilterYield => "filter_yield", "Yield after filtering";
    ReadLength => "read_length", "Read length distribution";
    ReadQuality => "read_quality", "Read quality distribution";
    LengthVsQuality => "length_vs_quality", "Read length vs quality";
    BaseQuality => "base_quality", "Per-position base quality";
    BaseContent => "base_content", "Per-position base content";
    GcContent => "gc_content", "Per-read GC content";
    NContent => "n_content", "Per-position N content";
    Duplication => "duplication", "Duplication levels";
    Overrepresented => "overrepresented", "Overrepresented end sequences";
    Adapters => "adapters", "Adapter and primer content";
    Contaminants => "contaminants", "Controls and contaminants";
    BarcodeCrosstalk => "barcode_crosstalk", "Barcodes inside reads";
    ReadIds => "read_ids", "Duplicate read IDs";
    BarcodeBalance => "barcode_balance", "Barcode balance";
    Transcripts => "transcripts", "cDNA and poly(A)";
    YieldOverTime => "yield_over_time", "Yield over time";
    Channels => "channels", "Channel activity";
    TranslocationSpeed => "translocation_speed", "Translocation speed";
    Duplex => "duplex", "Duplex";
    Signal => "signal", "Signal, muxes and split reads";
    EndReasons => "end_reasons", "Read end reasons";
    PacBio => "pacbio", "PacBio HiFi metrics";
    Kinnex => "kinnex", "Kinnex arrays";
    Modifications => "modifications", "Base modifications";
    Alignment => "alignment", "Alignment";
}

/// A module's verdict and the reason for it.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Module {
    pub id: ModuleId,
    pub status: Status,
    pub note: String,
}

fn module(id: ModuleId, status: Status, note: impl Into<String>) -> Module {
    Module {
        id,
        status,
        note: note.into(),
    }
}

/// Base content is graded from this position bin (position 21) on: every long-read run has
/// biased first and last bases, and the adapter check grades adapters left there.
const GRADED_FROM_POSITION: usize = 20;

/// Position bins with fewer bases are not graded: at 2,000 bases a 25% base share varies by
/// about one point by chance, far below the warn level.
const MIN_GRADED_BASES: u64 = 2_000;

/// A sample is taken to be the organism a contaminant reference holds when most reads match.
const MAJORITY_PCT: f64 = 50.0;

/// Warn and fail levels for read and base quality.
pub fn quality_levels(platform: Platform, min_q: Option<f32>, t: &Thresholds) -> (f64, f64) {
    let warn = min_q.map_or(
        match platform {
            Platform::PacBio => t.read_q_warn_pacbio,
            Platform::Ont | Platform::Unknown => t.read_q_warn_ont,
        },
        f64::from,
    );
    (warn, warn * t.read_q_fail_ratio)
}

/// Grades `value` against warn and fail levels, in the given direction.
fn grade(value: f64, warn: f64, fail: f64, direction: Direction) -> Status {
    let (warn_hit, fail_hit) = match direction {
        Direction::HigherIsWorse => (value > warn, value > fail),
        Direction::LowerIsWorse => (value < warn, value < fail),
    };
    if fail_hit {
        Status::Fail
    } else if warn_hit {
        Status::Warn
    } else {
        Status::Pass
    }
}

#[derive(Debug, Clone, Copy)]
enum Direction {
    HigherIsWorse,
    LowerIsWorse,
}

/// Runs every check on a group.
pub fn evaluate(
    stats: &Stats,
    summary: &Summary,
    data: &Data,
    platform: Platform,
    molecule: Molecule,
    min_q: Option<f32>,
    t: &Thresholds,
) -> Vec<Module> {
    let (q_warn, q_fail) = quality_levels(platform, min_q, t);
    let mut modules = evaluate_dna(stats, summary, data, q_warn, q_fail, t);
    if molecule == Molecule::Rna {
        relax_for_rna(&mut modules);
    }
    modules
}

/// Direct RNA breaks the assumptions of four checks: single-stranded transcripts need not
/// balance A with T or G with C, transcriptome GC has several peaks, and highly expressed
/// transcripts repeat by design. Those are shown without a verdict.
fn relax_for_rna(modules: &mut [Module]) {
    for m in modules.iter_mut().filter(|m| {
        matches!(
            m.id,
            ModuleId::BaseContent
                | ModuleId::GcContent
                | ModuleId::Duplication
                | ModuleId::Overrepresented
        ) && m.status != Status::Absent
    }) {
        m.status = Status::Info;
        m.note = format!("direct RNA: shown, not graded ({})", m.note);
    }
}

fn evaluate_dna(
    stats: &Stats,
    summary: &Summary,
    data: &Data,
    q_warn: f64,
    q_fail: f64,
    t: &Thresholds,
) -> Vec<Module> {
    let has_seq = stats.reads_with_seq > 0;
    let absent = |id: ModuleId, why: &str| module(id, Status::Absent, why);
    let info = |id: ModuleId, note: String| module(id, Status::Info, note);
    let mut modules = Vec::with_capacity(ModuleId::ALL.len());

    modules.push(if summary.reads == 0 {
        module(ModuleId::BasicStatistics, Status::Fail, "no reads")
    } else {
        let coverage = summary
            .coverage
            .map_or_else(String::new, |c| format!(", {c:.1}x coverage"));
        module(
            ModuleId::BasicStatistics,
            Status::Pass,
            format!("{} reads{coverage}", summary.reads),
        )
    });

    modules.push(match &data.filter_grid {
        Some(_) => info(ModuleId::FilterYield, String::new()),
        None => absent(ModuleId::FilterYield, "no reads"),
    });

    modules.push(if stats.lengths.contains(&0) {
        module(
            ModuleId::ReadLength,
            Status::Warn,
            "some reads have no bases",
        )
    } else {
        info(
            ModuleId::ReadLength,
            format!("read N50 {} bases", summary.n50),
        )
    });

    modules.push(match summary.median_read_q {
        // Most read Qs estimated from binned PacBio base qualities: not graded.
        Some(q) if stats.pacbio_q_estimated * 2 > stats.reads_with_q => info(
            ModuleId::ReadQuality,
            format!(
                "median read Q{q:.1}, estimated from base qualities: the input lacks PacBio's predicted accuracy (rq), and binned HiFi base qualities put the estimate several Q off, so it is not graded"
            ),
        ),
        Some(q) => module(
            ModuleId::ReadQuality,
            grade(q, q_warn, q_fail, Direction::LowerIsWorse),
            format!("median read Q{q:.1}; warn below Q{q_warn:.1}, fail below Q{q_fail:.1}"),
        ),
        None => absent(ModuleId::ReadQuality, "no qualities"),
    });

    modules.push(if summary.median_read_q.is_some() {
        info(ModuleId::LengthVsQuality, String::new())
    } else {
        absent(ModuleId::LengthVsQuality, "no qualities")
    });

    // Base quality: the median of every plotted position bin against the read Q levels,
    // and the lower quartile against the fail level, as FastQC uses its quartiles.
    let positions = data.quality_from_start.iter().chain(&data.quality_from_end);
    modules.push(match positions.clone().map(|p| p.median).reduce(f64::min) {
        Some(worst_median) => {
            let worst_quartile = positions.map(|p| p.p25).fold(f64::INFINITY, f64::min);
            let status = grade(worst_median, q_warn, q_fail, Direction::LowerIsWorse).max(
                grade(worst_quartile, q_fail, q_fail * 0.5, Direction::LowerIsWorse)
                    .min(Status::Warn),
            );
            module(
                ModuleId::BaseQuality,
                status,
                format!(
                    "lowest median Q{worst_median:.0}, lowest lower quartile Q{worst_quartile:.0}; warn when a median is below Q{q_warn:.1} or a lower quartile below Q{q_fail:.1}, fail when a median is below Q{q_fail:.1}"
                ),
            )
        }
        None => absent(ModuleId::BaseQuality, "no qualities"),
    });

    let graded = data
        .content_from_start
        .iter()
        .chain(&data.content_from_end)
        .filter(|c| c.bin >= GRADED_FROM_POSITION && c.bases >= MIN_GRADED_BASES);
    modules.push(
        match graded
            .map(|c| (c.a - c.t).abs().max((c.g - c.c).abs()))
            .reduce(f64::max)
        {
            Some(worst) => module(
                ModuleId::BaseContent,
                grade(worst, t.base_content_warn, t.base_content_fail, Direction::HigherIsWorse),
                format!(
                    "largest A/T or G/C imbalance {worst:.1}% beyond the first and last 20 bases; warn above {}%, fail above {}%",
                    t.base_content_warn, t.base_content_fail
                ),
            ),
            None if has_seq => info(
                ModuleId::BaseContent,
                "too few bases beyond the read ends to grade".to_owned(),
            ),
            None => absent(ModuleId::BaseContent, "no sequence"),
        },
    );

    // GC: long reads average over isochores, so their GC is skewed, not normal, and
    // FastQC's deviation measure flags clean human samples. A second GC peak is what
    // contamination or a mixed sample looks like.
    modules.push(match data.gc_deviation_pct {
        Some(d) => match &data.gc_second_peak {
            Some(peak) => module(
                ModuleId::GcContent,
                grade(
                    peak.share_pct,
                    t.gc_second_peak_warn,
                    t.gc_second_peak_fail,
                    Direction::HigherIsWorse,
                ),
                format!(
                    "second GC peak at {}% holding {:.1}% of reads{} (main peak {}%): contamination or a mixed sample; warn above {}%, fail above {}%",
                    peak.gc,
                    peak.share_pct,
                    if peak.min_length > 0 { format!(" of {} bases or more", peak.min_length) } else { String::new() },
                    peak.main_gc,
                    t.gc_second_peak_warn,
                    t.gc_second_peak_fail
                ),
            ),
            None => module(
                ModuleId::GcContent,
                Status::Pass,
                format!("one GC peak; {d:.1}% departure from a normal curve (shown, not graded)"),
            ),
        },
        None => absent(ModuleId::GcContent, "no sequence"),
    });

    let n_worst = data
        .content_from_start
        .iter()
        .chain(&data.content_from_end)
        .map(|c| c.n)
        .reduce(f64::max);
    modules.push(match n_worst {
        Some(n) => module(
            ModuleId::NContent,
            grade(
                n,
                t.n_content_warn,
                t.n_content_fail,
                Direction::HigherIsWorse,
            ),
            format!(
                "highest N content {n:.2}%; warn above {}%, fail above {}%",
                t.n_content_warn, t.n_content_fail
            ),
        ),
        None => absent(ModuleId::NContent, "no sequence"),
    });

    modules.push(match &data.duplication {
        Some(d) => {
            let duplicated = 100.0 - d.distinct_pct;
            module(
                ModuleId::Duplication,
                grade(duplicated, t.duplication_warn, t.duplication_fail, Direction::HigherIsWorse),
                format!(
                    "{duplicated:.1}% of {} sampled reads are copies of another read's 50-base ends (FastQC's measure: copies beyond the first); warn above {}%, fail above {}% (expected for amplicons)",
                    d.sampled, t.duplication_warn, t.duplication_fail
                ),
            )
        }
        None => absent(ModuleId::Duplication, "no reads of 50 bases or more without N"),
    });

    modules.push(if data.duplication.is_none() {
        absent(
            ModuleId::Overrepresented,
            "no reads of 50 bases or more without N",
        )
    } else {
        let top = data
            .overrepresented
            .iter()
            .map(|o| o.pct)
            .fold(0.0, f64::max);
        let status = if data.overrepresented.is_empty() {
            Status::Pass
        } else {
            grade(
                top,
                t.overrepresented_warn,
                t.overrepresented_fail,
                Direction::HigherIsWorse,
            )
        };
        module(
            ModuleId::Overrepresented,
            status,
            format!(
                "{} end sequences above {}% of sampled reads; warn above {}%, fail above {}%",
                data.overrepresented.len(),
                super::model::listing_pct(t),
                t.overrepresented_warn,
                t.overrepresented_fail
            ),
        )
    });

    modules.push(if !has_seq {
        absent(ModuleId::Adapters, "no sequence")
    } else {
        // Reads with any adapter or primer: one read with adapter A at its start and
        // another with adapter B at its end are two untrimmed reads.
        let at_ends = data.adapter_at_ends_pct;
        let internal = data.adapter_internal_pct;
        let status = grade(at_ends, t.adapter_end_warn, t.adapter_end_fail, Direction::HigherIsWorse)
            .max(grade(
                internal,
                t.adapter_internal_warn,
                t.adapter_internal_fail,
                Direction::HigherIsWorse,
            ));
        module(
            ModuleId::Adapters,
            status,
            format!(
                "{at_ends:.2}% of reads with an adapter or primer within 150 bases of an end, {internal:.2}% with one inside (chimeras); warn above {}% and {}%, fail above {}% and {}%",
                t.adapter_end_warn, t.adapter_internal_warn, t.adapter_end_fail, t.adapter_internal_fail
            ),
        )
    });

    modules.push(contaminants(data, t));
    modules.push(barcodes_inside(data, t));

    modules.push(if summary.duplicate_names == 0 {
        module(ModuleId::ReadIds, Status::Pass, "every read ID is unique")
    } else {
        let p = 100.0 * summary.duplicate_names as f64 / summary.reads.max(1) as f64;
        module(
            ModuleId::ReadIds,
            if p > t.duplicate_ids_fail { Status::Fail } else { Status::Warn },
            format!(
                "{} reads repeat an earlier read ID ({p:.2}%): files given twice, or overlapping inputs; warn on any, fail above {}%",
                summary.duplicate_names, t.duplicate_ids_fail
            ),
        )
    });

    // Barcode balance is set per barcode once all groups are built.
    modules.push(absent(ModuleId::BarcodeBalance, "not demultiplexed"));

    modules.push(match (&data.cdna, data.poly_a.is_empty()) {
        (Some(c), _) => info(
            ModuleId::Transcripts,
            format!(
                "{:.1}% of reads full-length (both cDNA primers)",
                c.full_length_pct
            ),
        ),
        (None, false) => info(
            ModuleId::Transcripts,
            "poly(A) lengths from dorado".to_owned(),
        ),
        (None, true) => absent(
            ModuleId::Transcripts,
            "no cDNA primers or poly(A) estimates",
        ),
    });

    modules.push(if data.time.is_empty() {
        absent(
            ModuleId::YieldOverTime,
            "no read start times (ONT st tags or summaries)",
        )
    } else {
        info(ModuleId::YieldOverTime, String::new())
    });
    modules.push(if data.channels.is_empty() {
        absent(
            ModuleId::Channels,
            "no channel numbers (ONT ch tags or summaries)",
        )
    } else {
        info(
            ModuleId::Channels,
            format!("{} active channels", data.channels.len()),
        )
    });
    modules.push(if data.speed.is_empty() {
        absent(
            ModuleId::TranslocationSpeed,
            "no read durations (ONT du tags or summaries)",
        )
    } else {
        info(ModuleId::TranslocationSpeed, String::new())
    });
    let ont = data.ont.as_ref();
    modules.push(match ont.and_then(|o| o.duplex.as_ref()) {
        Some(d) => info(
            ModuleId::Duplex,
            format!("duplex rate {:.1}%", d.duplex_rate_pct),
        ),
        None => absent(ModuleId::Duplex, "no dx tags (dorado duplex)"),
    });
    modules.push(match ont {
        Some(o) if !o.shift.is_empty() || !o.mux.is_empty() || o.split_reads > 0 => info(
            ModuleId::Signal,
            format!("{:.2}% of reads split by the basecaller", o.split_pct),
        ),
        _ => absent(ModuleId::Signal, "no ONT signal, mux or split tags"),
    });
    modules.push(if data.end_reasons.is_empty() {
        absent(
            ModuleId::EndReasons,
            "no end reasons (dorado er tags or sequencing summaries)",
        )
    } else {
        info(ModuleId::EndReasons, String::new())
    });
    modules.push(match &data.pacbio {
        Some(p) => info(
            ModuleId::PacBio,
            p.median_passes
                .map_or_else(String::new, |n| format!("median {n} passes")),
        ),
        None => absent(ModuleId::PacBio, "no PacBio tags (np, ec, bq, kinetics)"),
    });
    modules.push(if data.kinnex.is_empty() {
        absent(ModuleId::Kinnex, "no --kinnex-adapters, or no arrays found")
    } else {
        info(ModuleId::Kinnex, String::new())
    });
    modules.push(match summary.mod_reads_pct {
        Some(p) => info(
            ModuleId::Modifications,
            format!("{p:.1}% of reads carry MM/ML calls"),
        ),
        None => absent(ModuleId::Modifications, "no MM/ML tags"),
    });
    modules.push(match &data.mapping {
        Some(m) => module(
            ModuleId::Alignment,
            grade(
                m.mapped_pct,
                t.mapped_warn,
                t.mapped_fail,
                Direction::LowerIsWorse,
            ),
            format!(
                "{:.2}% of reads mapped; warn below {}%, fail below {}%",
                m.mapped_pct, t.mapped_warn, t.mapped_fail
            ),
        ),
        None => absent(ModuleId::Alignment, "unaligned input"),
    });
    debug_assert_eq!(modules.len(), ModuleId::ALL.len(), "one result per module");
    modules
}

/// Contaminants are graded on their share of reads; controls and the mitochondrion are
/// reported. A sample that mostly matches a "contaminant" is taken to be that organism.
fn contaminants(data: &Data, t: &Thresholds) -> Module {
    if data.reads_screened == 0 {
        return module(
            ModuleId::Contaminants,
            Status::Absent,
            "no sequence, or --no-screen",
        );
    }
    let hits = &data.contaminants;
    // Controls, organelles and vectors are named, not graded.
    let named: Vec<String> = [Kind::Control, Kind::Organelle, Kind::Vector]
        .into_iter()
        .flat_map(|kind| {
            hits.iter()
                .filter(move |h| h.reference.kind == kind && h.reads > 0)
                .map(|h| format!("{} {:.2}%", h.reference.name, h.reads_pct))
        })
        .collect();
    let others: String = named.iter().map(|n| format!("; {n}")).collect();
    let (sample, worst) = split_contaminants(hits);
    // What a sparse reference's share is of.
    let of_reads = |w: &ScreenHit| match w.testable {
        Some(n) => format!("of the {n} reads long enough to test"),
        None if w.lower_bound => {
            "of all reads (a lower bound: too few reads were long enough to test)".to_owned()
        }
        None => "of reads".to_owned(),
    };
    // A sample that mostly matches a "contaminant" is that organism; the others are still
    // graded, so E. coli does not hide phiX, nor human hide Mycoplasma.
    let organism = sample.map(|s| {
        format!(
            "{:.1}% {} match {}: the sample is likely {}",
            s.reads_pct,
            of_reads(s),
            s.reference.name,
            s.reference.name
        )
    });
    match worst {
        Some(w) if w.reads > 0 => {
            let status = grade(
                w.reads_pct,
                t.contaminant_warn,
                t.contaminant_fail,
                Direction::HigherIsWorse,
            );
            // A lower bound cannot tell contamination from a sample of that organism whose
            // short reads went untested (human RNA), so it warns at most.
            let (status, caveat) = if w.lower_bound && status == Status::Fail {
                (
                    Status::Warn,
                    "; a warning only: the untested short reads may be the sample's own",
                )
            } else {
                (status, "")
            };
            let prefix = organism.map(|o| format!("{o}; ")).unwrap_or_default();
            module(
                ModuleId::Contaminants,
                status,
                format!(
                    "{prefix}{:.2}% {} match {}; warn above {}%, fail above {}%{caveat}{others}",
                    w.reads_pct,
                    of_reads(w),
                    w.reference.name,
                    t.contaminant_warn,
                    t.contaminant_fail
                ),
            )
        }
        _ => match organism {
            Some(o) => module(ModuleId::Contaminants, Status::Info, format!("{o}{others}")),
            None => module(
                ModuleId::Contaminants,
                Status::Pass,
                format!("no reads match a contaminant{others}"),
            ),
        },
    }
}

/// Splits the contaminant hits into the organism the sample is (one matching most reads),
/// if any, and the worst of the rest, which is graded.
pub fn split_contaminants(hits: &[ScreenHit]) -> (Option<&ScreenHit>, Option<&ScreenHit>) {
    let by_share = |a: &&ScreenHit, b: &&ScreenHit| a.reads_pct.total_cmp(&b.reads_pct);
    let contaminants = || {
        hits.iter()
            .filter(|h| h.reference.kind == Kind::Contaminant)
    };
    let sample = contaminants()
        .max_by(by_share)
        .filter(|h| h.reads_pct > MAJORITY_PCT);
    let worst = contaminants()
        .filter(|h| !sample.is_some_and(|s| std::ptr::eq(*h, s)))
        .max_by(by_share);
    (sample, worst)
}

/// Share of reads from a contaminant other than the organism the sample is, for the
/// summary tables.
pub fn foreign_contaminant_pct(hits: &[ScreenHit]) -> Option<f64> {
    match split_contaminants(hits) {
        (_, Some(w)) => Some(w.reads_pct),
        (Some(_), None) => Some(0.0),
        (None, None) => None,
    }
}

/// Barcodes left inside reads: the read's own at an end means barcodes were not trimmed;
/// another sample's means cross-talk or chimeras across samples.
fn barcodes_inside(data: &Data, t: &Thresholds) -> Module {
    let Some(scan) = &data.barcode_scan else {
        return module(
            ModuleId::BarcodeCrosstalk,
            Status::Absent,
            "no barcode sequences for this platform (PacBio: barcodes named in file names or the sample sheet, or --barcode-names)",
        );
    };
    let status = grade(
        scan.own_at_end_pct,
        t.untrimmed_barcode_warn,
        t.untrimmed_barcode_fail,
        Direction::HigherIsWorse,
    )
    .max(grade(
        scan.foreign_pct,
        t.foreign_barcode_warn,
        t.foreign_barcode_fail,
        Direction::HigherIsWorse,
    ));
    module(
        ModuleId::BarcodeCrosstalk,
        status,
        format!(
            "{:.2}% of reads keep their barcode at an end, {:.2}% hold another sample's barcode; warn above {}% and {}%, fail above {}% and {}%",
            scan.own_at_end_pct,
            scan.foreign_pct,
            t.untrimmed_barcode_warn,
            t.foreign_barcode_warn,
            t.untrimmed_barcode_fail,
            t.foreign_barcode_fail
        ),
    )
}

/// Flags barcodes with far fewer reads than the median barcode.
pub fn barcode_balance(group: &mut Group, median_reads: u64, t: &Thresholds) {
    if group.kind != GroupKind::Barcode || median_reads == 0 {
        return;
    }
    let ratio = group.summary.reads as f64 / median_reads as f64;
    let status = grade(
        ratio,
        t.barcode_ratio_warn,
        t.barcode_ratio_fail,
        Direction::LowerIsWorse,
    );
    set(
        group,
        ModuleId::BarcodeBalance,
        status,
        format!(
            "{ratio:.2}x the median barcode's reads; warn below {}x, fail below {}x",
            t.barcode_ratio_warn, t.barcode_ratio_fail
        ),
    );
}

/// Grades unclassified groups by their share of the run, and makes foreign barcodes in
/// them informational: unassigned reads holding barcodes are what unclassified means.
pub fn unclassified_share(group: &mut Group, share: f64, t: &Thresholds) {
    if group.kind == GroupKind::Unclassified {
        set(
            group,
            ModuleId::BarcodeBalance,
            grade(
                share,
                t.unclassified_warn,
                t.unclassified_fail,
                Direction::HigherIsWorse,
            ),
            format!(
                "{share:.1}% of reads unclassified; warn above {}%, fail above {}%",
                t.unclassified_warn, t.unclassified_fail
            ),
        );
        if let Some(scan) = &group.data.barcode_scan {
            let note = format!(
                "{:.2}% of unclassified reads hold a barcode the demultiplexer did not assign",
                scan.foreign_pct
            );
            set(group, ModuleId::BarcodeCrosstalk, Status::Info, note);
        }
    }
}

/// In a demultiplexed run, all reads and unclassified reads pool every sample: each sample's
/// organism shows up there as a "contaminant", and a pool of different organisms has several
/// GC peaks by design. Both are graded per barcode; in the pools they are shown without a
/// verdict.
pub fn pools_not_graded_for_mixtures(groups: &mut [Group]) {
    for group in groups
        .iter_mut()
        .filter(|g| matches!(g.kind, GroupKind::AllReads | GroupKind::Unclassified))
    {
        for id in [ModuleId::Contaminants, ModuleId::GcContent] {
            if let Some(m) = group.module(id) {
                if matches!(m.status, Status::Pass | Status::Warn | Status::Fail) {
                    let note = format!(
                        "pooled reads of every sample: shown, not graded ({})",
                        m.note
                    );
                    set(group, id, Status::Info, note);
                }
            }
        }
    }
}

fn set(group: &mut Group, id: ModuleId, status: Status, note: String) {
    if let Some(m) = group.modules.iter_mut().find(|m| m.id == id) {
        m.status = status;
        m.note = note;
    }
}

/// The most severe verdict among graded modules of the given groups.
pub fn worst<'a>(groups: impl IntoIterator<Item = &'a Group>) -> Option<(Status, String)> {
    groups
        .into_iter()
        .flat_map(|g| {
            g.modules
                .iter()
                .filter(|m| matches!(m.status, Status::Warn | Status::Fail))
                .map(move |m| {
                    (
                        m.status,
                        format!("{}: {} ({})", g.name, m.id.title(), m.note),
                    )
                })
        })
        .max_by_key(|(status, _)| *status)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn grades_contaminants_beside_the_sample_organism() {
        let mut hits = vec![
            hit("E. coli", Kind::Contaminant, 9_090, 90.9, None),
            hit("phiX", Kind::Contaminant, 910, 9.1, None),
        ];
        let m = screened(hits.clone());
        assert_eq!(
            m.status,
            Status::Fail,
            "the majority must not hide phiX: {m:?}"
        );
        assert!(m.note.contains("likely E. coli") && m.note.contains("phiX"));
        assert_eq!(foreign_contaminant_pct(&hits), Some(9.1));

        hits.truncate(1);
        assert_eq!(screened(hits.clone()).status, Status::Info);
        assert_eq!(foreign_contaminant_pct(&hits), Some(0.0));
    }

    #[test]
    fn grades_both_directions() {
        assert_eq!(
            grade(5.0, 10.0, 20.0, Direction::HigherIsWorse),
            Status::Pass
        );
        assert_eq!(
            grade(15.0, 10.0, 20.0, Direction::HigherIsWorse),
            Status::Warn
        );
        assert_eq!(
            grade(25.0, 10.0, 20.0, Direction::HigherIsWorse),
            Status::Fail
        );
        assert_eq!(grade(9.0, 10.0, 7.5, Direction::LowerIsWorse), Status::Warn);
        assert_eq!(grade(7.0, 10.0, 7.5, Direction::LowerIsWorse), Status::Fail);
    }

    #[test]
    fn quality_levels_follow_platform() {
        let t = Thresholds::default();
        assert_eq!(quality_levels(Platform::Ont, None, &t), (10.0, 7.5));
        assert_eq!(quality_levels(Platform::PacBio, None, &t), (20.0, 15.0));
        assert_eq!(
            quality_levels(Platform::PacBio, Some(30.0), &t),
            (30.0, 22.5)
        );
    }

    #[test]
    fn relaxes_dna_assumptions_for_rna() {
        let mut modules = vec![
            module(ModuleId::BaseContent, Status::Fail, "imbalance"),
            module(ModuleId::Overrepresented, Status::Fail, "ends"),
            module(ModuleId::ReadQuality, Status::Fail, "low Q"),
            module(ModuleId::GcContent, Status::Absent, "no sequence"),
        ];
        relax_for_rna(&mut modules);
        let statuses: Vec<Status> = modules.iter().map(|m| m.status).collect();
        assert_eq!(
            statuses,
            [Status::Info, Status::Info, Status::Fail, Status::Absent]
        );
    }

    /// A screen hit; `basis` is `Some(Ok(n))` for a share of `n` testable reads, and
    /// `Some(Err(()))` for a lower bound of all reads.
    fn hit(
        key: &str,
        kind: Kind,
        reads: u64,
        reads_pct: f64,
        basis: Option<Result<u64, ()>>,
    ) -> ScreenHit {
        ScreenHit {
            reference: crate::screen::Reference {
                key: key.to_owned(),
                name: key.to_owned(),
                kind,
                accession: String::new(),
                length: 0,
                sparse: basis.is_some(),
            },
            reads,
            bases: 0,
            reads_pct,
            bases_pct: 0.0,
            testable: basis.and_then(Result::ok),
            lower_bound: basis.is_some_and(|b| b.is_err()),
        }
    }

    fn screened(hits: Vec<ScreenHit>) -> Module {
        let data = Data {
            reads_screened: 10_000,
            contaminants: hits,
            ..Data::default()
        };
        contaminants(&data, &Thresholds::default())
    }

    #[test]
    fn grades_sparse_references_by_what_their_share_is_of() {
        // Long human reads (HG002 ONT): nearly every read long enough to test is human.
        let human = screened(vec![hit(
            "human",
            Kind::Contaminant,
            9_000,
            98.8,
            Some(Ok(7_000)),
        )]);
        assert_eq!(human.status, Status::Info);
        assert!(
            human.note.contains("of the 7000 reads long enough to test"),
            "{}",
            human.note
        );
        // A 1.5 kb amplicon run with 2% long human reads: few reads could be tested, so the
        // share is of all reads, and the contamination is flagged, not taken for the sample.
        let amplicon = screened(vec![hit(
            "human",
            Kind::Contaminant,
            200,
            2.0,
            Some(Err(())),
        )]);
        assert_eq!(amplicon.status, Status::Warn);
        assert!(amplicon.note.contains("a lower bound"), "{}", amplicon.note);
        // Human RNA, mostly short reads: 10% is past the fail level, but the untested short
        // reads may be human too, so it only warns.
        let rna = screened(vec![hit(
            "human",
            Kind::Contaminant,
            1_000,
            10.0,
            Some(Err(())),
        )]);
        assert_eq!(rna.status, Status::Warn);
        assert!(rna.note.contains("may be the sample's own"), "{}", rna.note);
        // A lower bound past half the reads is still a majority.
        let most = screened(vec![hit(
            "human",
            Kind::Contaminant,
            6_000,
            60.0,
            Some(Err(())),
        )]);
        assert_eq!(most.status, Status::Info);
        // A dense contaminant is graded as before, and vectors are named, not graded.
        let ecoli = screened(vec![
            hit("ecoli", Kind::Contaminant, 300, 3.0, None),
            hit("univec", Kind::Vector, 900, 9.0, None),
        ]);
        assert_eq!(ecoli.status, Status::Warn);
        assert!(ecoli.note.contains("; univec 9.00%"), "{}", ecoli.note);
        let dense_fail = screened(vec![hit("ecoli", Kind::Contaminant, 1_000, 10.0, None)]);
        assert_eq!(dense_fail.status, Status::Fail);
    }

    #[test]
    fn module_keys_are_unique() {
        let mut keys: Vec<&str> = ModuleId::ALL.iter().map(|m| m.key()).collect();
        keys.sort_unstable();
        keys.dedup();
        assert_eq!(keys.len(), ModuleId::ALL.len());
    }
}
