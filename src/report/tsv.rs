//! One row per group: headline numbers and every module's status.

use std::fmt::Write as _;

use super::checks::ModuleId;
use super::Report;

/// Renders `summary.tsv`.
pub fn render(report: &Report) -> String {
    let mut out = String::new();
    let columns = [
        "group",
        "sample",
        "platform",
        "reads",
        "bases",
        "share_pct",
        "mean_length",
        "median_length",
        "n50",
        "n90",
        "longest",
        "median_read_q",
        "mean_read_q",
        "reads_q10_pct",
        "reads_q20_pct",
        "reads_q30_pct",
        "bases_q20",
        "base_q",
        "gc_pct",
        "pass_reads",
        "pass_bases",
        "pass_pct",
        "mapped_pct",
        "median_identity_pct",
        "duplicate_read_ids",
        "mod_reads_pct",
        "cpg_5mc_pct",
        "coverage",
        "contaminant_pct",
    ];
    out.push_str(&columns.join("\t"));
    for &id in ModuleId::ALL {
        let _ = write!(out, "\t{}", id.key());
    }
    out.push('\n');

    let opt = |v: Option<f64>| v.map_or_else(|| "NA".to_owned(), |v| format!("{v:.2}"));
    for group in &report.groups {
        let s = &group.summary;
        let fields = [
            group.name.clone(),
            group.sample.clone().unwrap_or_default(),
            group.platform.label().to_owned(),
            s.reads.to_string(),
            s.bases.to_string(),
            format!("{:.2}", s.share_pct),
            format!("{:.1}", s.mean_length),
            format!("{:.1}", s.median_length),
            s.n50.to_string(),
            s.n90.to_string(),
            s.longest.to_string(),
            opt(s.median_read_q),
            opt(s.mean_read_q),
            opt(s.reads_q10_pct),
            opt(s.reads_q20_pct),
            opt(s.reads_q30_pct),
            s.bases_q20
                .map_or_else(|| "NA".to_owned(), |b| b.to_string()),
            opt(s.base_q),
            opt(s.gc_pct),
            s.pass_reads.to_string(),
            s.pass_bases.to_string(),
            opt(s.pass_pct),
            opt(s.mapped_pct),
            opt(s.median_identity_pct),
            s.duplicate_names.to_string(),
            opt(s.mod_reads_pct),
            opt(s.methylated_pct),
            opt(s.coverage),
            opt(super::checks::foreign_contaminant_pct(
                &group.data.contaminants,
            )),
        ];
        out.push_str(&fields.join("\t"));
        for &id in ModuleId::ALL {
            let status = group.module(id).map_or("NA", |m| m.status.label());
            let _ = write!(out, "\t{status}");
        }
        out.push('\n');
    }
    out
}
