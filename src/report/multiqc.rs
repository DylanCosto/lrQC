//! MultiQC custom content: lrqc's numbers as general statistics and its checks as a table.
//!
//! MultiQC picks up files ending in `_mqc.tsv` and reads the `#` header lines as settings.

use std::fmt::Write as _;

use super::checks::{ModuleId, Status};
use super::{GroupKind, Report};

/// MultiQC's names for the groups, in order: the report title, then the sample name if
/// known, else the group name. MultiQC keeps one row per name, so a sample name shared by
/// several groups (one alias for two barcodes, or one sample over runs with --split-runs)
/// falls back to the group name with the sample in brackets.
fn sample_names(report: &Report) -> Vec<String> {
    let short = |group: &super::Group| match (group.kind, &group.sample) {
        (GroupKind::AllReads, _) => report.title.clone(),
        (_, Some(sample)) => format!("{} {sample}", report.title),
        (_, None) => format!("{} {}", report.title, group.name),
    };
    let names: Vec<String> = report.groups.iter().map(short).collect();
    report
        .groups
        .iter()
        .zip(&names)
        .map(|(group, name)| {
            let shared = names.iter().filter(|n| *n == name).count() > 1;
            match &group.sample {
                Some(sample) if shared => format!("{} {} ({sample})", report.title, group.name),
                _ => name.clone(),
            }
        })
        .collect()
}

/// `lrqc_general_stats_mqc.tsv`: one row per group in MultiQC's general statistics table.
pub fn general_stats(report: &Report) -> String {
    let mut out = String::from(
        // MultiQC shows only the general statistics columns its pconfig names.
        "# id: 'lrqc'\n# section_name: 'lrQC'\n# plot_type: 'generalstats'\n# pconfig:\n#     reads:\n#         title: 'Reads'\n#         format: '{:,.0f}'\n#     bases:\n#         title: 'Bases'\n#         format: '{:,.0f}'\n#     n50:\n#         title: 'Read N50'\n#         format: '{:,.0f}'\n#     median_read_q:\n#         title: 'Median Q'\n#         format: '{:,.1f}'\n#     pass_pct:\n#         title: '% Pass'\n#         suffix: '%'\n#         max: 100\n#         format: '{:,.1f}'\n#     gc_pct:\n#         title: '% GC'\n#         suffix: '%'\n#         max: 100\n#         format: '{:,.1f}'\nSample\treads\tbases\tn50\tmedian_read_q\tpass_pct\tgc_pct\n",
    );
    let opt = |v: Option<f64>| v.map_or_else(String::new, |v| format!("{v:.2}"));
    for (group, name) in report.groups.iter().zip(sample_names(report)) {
        let s = &group.summary;
        let _ = writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            name,
            s.reads,
            s.bases,
            s.n50,
            opt(s.median_read_q),
            opt(s.pass_pct),
            opt(s.gc_pct)
        );
    }
    out
}

/// `lrqc_checks_mqc.tsv`: every check's verdict per group, as a MultiQC table.
pub fn checks(report: &Report) -> String {
    let shown: Vec<ModuleId> = ModuleId::ALL
        .iter()
        .copied()
        .filter(|&id| {
            report.groups.iter().any(|g| {
                g.module(id)
                    .is_some_and(|m| matches!(m.status, Status::Pass | Status::Warn | Status::Fail))
            })
        })
        .collect();
    let mut out = String::from(
        "# id: 'lrqc_checks'\n# section_name: 'lrQC checks'\n# description: 'pass, warn or fail for each lrQC check'\n# plot_type: 'table'\nSample",
    );
    for id in &shown {
        let _ = write!(out, "\t{}", id.title());
    }
    out.push('\n');
    for (group, name) in report.groups.iter().zip(sample_names(report)) {
        out.push_str(&name);
        for &id in &shown {
            let status = group.module(id).map_or("", |m| m.status.label());
            let _ = write!(out, "\t{status}");
        }
        out.push('\n');
    }
    out
}
