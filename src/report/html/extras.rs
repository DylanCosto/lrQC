//! Module bodies for filtering, contaminants, barcodes inside reads, transcripts, ONT duplex
//! and signal, Kinnex arrays, and the run information panel.

use std::fmt::Write as _;

use super::{figure, histogram_chart, long_table_wrap, opt, thousands, Page};
use crate::report::chart::{escape, si, Axis, Bar, Chart, Layer, Units};
use crate::report::model::Group;
use crate::report::Report;
use crate::screen::Kind;

pub(super) fn filter_yield(page: &mut Page, report: &Report, group: &Group) {
    let Some(grid) = &group.data.filter_grid else {
        return;
    };
    for (title, values) in [("Reads", &grid.reads), ("Bases", &grid.bases)] {
        let _ = write!(
            page,
            r#"<div class="table-wrap"><table><tr><th class="l">{title}: min length ↓ / min Q →</th>"#
        );
        for q in &grid.min_qs {
            let _ = write!(page, "<th>Q{q:.0}</th>");
        }
        page.push_str("</tr>");
        let total = values[0][0].max(1) as f64;
        for (row, &min_length) in values.iter().zip(&grid.min_lengths) {
            let _ = write!(
                page,
                r#"<tr><td class="l">≥ {}</td>"#,
                si(min_length as f64)
            );
            for &value in row {
                let shown = match (title, report.genome_size) {
                    ("Bases", Some(genome)) => format!(
                        "{} ({:.1}x)",
                        si(value as f64),
                        value as f64 / genome as f64
                    ),
                    ("Bases", None) => si(value as f64),
                    _ => thousands(value),
                };
                let share = 100.0 * value as f64 / total;
                let _ = write!(page, r#"<td title="{share:.1}% kept">{shown}</td>"#);
            }
            page.push_str("</tr>");
        }
        page.push_str("</table></div>");
    }
}

pub(super) fn contaminants(page: &mut Page, group: &Group) {
    let data = &group.data;
    if data.contaminants.is_empty() {
        return;
    }
    let _ = write!(
        page,
        r#"<div class="table-wrap"><table><tr><th class="l">Reference</th><th class="l">Kind</th><th class="l">Accession</th><th>Reads</th><th>% of {} reads screened</th><th>Bases</th><th>% of bases</th></tr>"#,
        thousands(data.reads_screened)
    );
    for hit in &data.contaminants {
        let kind = match hit.reference.kind {
            Kind::Control => "spike-in control",
            Kind::Contaminant => "contaminant",
            Kind::Organelle => "organelle",
            Kind::Vector => "vector or plasmid backbone",
        };
        let mark = match (hit.testable, hit.lower_bound) {
            (Some(_), _) => "*",
            (None, true) => "†",
            (None, false) => "",
        };
        let _ = write!(
            page,
            r#"<tr><td class="l">{}</td><td class="l">{kind}</td><td class="l"><code>{}</code></td><td>{}</td><td>{:.3}%{mark}</td><td>{}</td><td>{:.3}%</td></tr>"#,
            escape(&hit.reference.name),
            escape(&hit.reference.accession),
            thousands(hit.reads),
            hit.reads_pct,
            si(hit.bases as f64),
            hit.bases_pct
        );
    }
    page.push_str("</table></div>");
    if let Some(n) = data.contaminants.iter().find_map(|h| h.testable) {
        let _ = write!(
            page,
            r#"<p class="note">* Of the {} reads long enough to test for human (about 7 kb and up): the human genome is sampled sparsely, so shorter reads are not called either way.</p>"#,
            thousands(n)
        );
    }
    if data.contaminants.iter().any(|h| h.lower_bound) {
        page.push_str(r#"<p class="note">† Of all reads, at least: too few reads were long enough (about 7 kb and up) to test for human, and the shorter ones were not tested, so more of them may be human.</p>"#);
    }
}

pub(super) fn barcodes_inside(page: &mut Page, group: &Group) {
    let Some(scan) = &group.data.barcode_scan else {
        return;
    };
    let _ = write!(
        page,
        r#"<div class="table-wrap"><table><tr><th class="l">Barcodes inside {} reads</th><th>% of reads</th></tr><tr><td class="l">Own barcode within 150 bases of an end (untrimmed)</td><td>{:.3}%</td></tr><tr><td class="l">Own barcode further inside (chimera within the sample)</td><td>{:.3}%</td></tr><tr><td class="l">Another sample's barcode</td><td>{:.3}%</td></tr></table></div>"#,
        thousands(scan.scanned),
        scan.own_at_end_pct,
        scan.own_internal_pct,
        scan.foreign_pct
    );
    if !scan.foreign.is_empty() {
        page.push_str(r#"<div class="table-wrap"><table><tr><th class="l">Other barcode found</th><th>Reads</th><th>%</th></tr>"#);
        for (name, reads, share) in &scan.foreign {
            let _ = write!(
                page,
                r#"<tr><td class="l">{}</td><td>{}</td><td>{share:.3}%</td></tr>"#,
                escape(name),
                thousands(*reads)
            );
        }
        page.push_str("</table></div>");
    }
}

pub(super) fn transcripts(page: &mut Page, group: &Group) {
    let data = &group.data;
    if let Some(c) = &data.cdna {
        let _ = write!(
            page,
            r#"<div class="table-wrap"><table><tr><th class="l">cDNA structure</th><th>% of reads</th></tr><tr><td class="l">Full length, SSP first (read from the transcript's 5' end)</td><td>{:.2}%</td></tr><tr><td class="l">Full length, VNP first (read from the poly(A) end)</td><td>{:.2}%</td></tr><tr><td class="l">SSP only (truncated at the 3' end)</td><td>{:.2}%</td></tr><tr><td class="l">VNP only (truncated at the 5' end)</td><td>{:.2}%</td></tr><tr><td class="l"><b>Full length</b></td><td><b>{:.2}%</b></td></tr></table></div>"#,
            c.forward_pct, c.reverse_pct, c.ssp_only_pct, c.vnp_only_pct, c.full_length_pct
        );
    }
    let total: u64 = data.poly_a.iter().sum();
    if total > 0 {
        let width = f64::from(crate::stats::POLY_A_BIN);
        let pcts: Vec<f64> = data
            .poly_a
            .iter()
            .map(|&c| 100.0 * c as f64 / total as f64)
            .collect();
        let top = pcts.iter().copied().fold(0.0, f64::max) * 1.08;
        let bars = pcts
            .iter()
            .enumerate()
            .map(|(i, &p)| {
                let lo = i as f64 * width;
                Bar {
                    x0: lo,
                    x1: lo + width,
                    y: p,
                    tip: format!("{lo:.0}–{:.0} A: {p:.2}% of reads", lo + width),
                }
            })
            .collect();
        let chart = Chart::new(
            "Poly(A) length",
            Axis::linear(
                "poly(A) length (bases; the last bin holds 600 and more)",
                0.0,
                pcts.len() as f64 * width,
                Units::Plain,
            ),
            Axis::linear("% of reads with an estimate", 0.0, top, Units::Percent),
        )
        .layer(Layer::Bars { series: 4, bars });
        figure(page, "dorado poly(A) estimates", chart, &[]);
    }
}

pub(super) fn duplex(page: &mut Page, group: &Group) {
    let Some(d) = group.data.ont.as_ref().and_then(|o| o.duplex.as_ref()) else {
        return;
    };
    page.push_str(r#"<div class="table-wrap"><table><tr><th class="l">Class</th><th>Reads</th><th>Bases</th><th>Median read Q</th></tr>"#);
    for (name, reads, bases, q) in &d.classes {
        let _ = write!(
            page,
            r#"<tr><td class="l">{}</td><td>{}</td><td>{}</td><td>{}</td></tr>"#,
            escape(name),
            thousands(*reads),
            si(*bases as f64),
            opt(*q, 1, "")
        );
    }
    let _ = write!(
        page,
        r#"</table></div><p class="note">Duplex rate {:.2}%. Yield above counts duplex bases and their simplex parents both; for unique bases, use the simplex reads.</p>"#,
        d.duplex_rate_pct
    );
}

pub(super) fn signal(page: &mut Page, group: &Group) {
    let Some(o) = &group.data.ont else { return };
    let _ = write!(
        page,
        r#"<p class="note">{} reads ({:.2}%) were split from a longer read by the basecaller.</p>"#,
        thousands(o.split_reads),
        o.split_pct
    );
    if !o.mux.is_empty() {
        let total: u64 = o.mux.iter().map(|m| m.1).sum();
        page.push_str(r#"<div class="table-wrap"><table><tr><th class="l">Mux</th><th>Reads</th><th>Bases</th><th>% of reads</th></tr>"#);
        for &(mux, reads, bases) in &o.mux {
            let _ = write!(
                page,
                r#"<tr><td class="l">{mux}</td><td>{}</td><td>{}</td><td>{:.1}%</td></tr>"#,
                thousands(reads),
                si(bases as f64),
                100.0 * reads as f64 / total.max(1) as f64
            );
        }
        page.push_str("</table></div>");
    }
    for (title, label, values, series) in [
        ("Signal shift", "shift (sm)", &o.shift, 0),
        ("Signal scale", "scale (sd)", &o.scale, 1),
        (
            "Samples per base",
            "samples per base (ns / length)",
            &o.samples_per_base,
            2,
        ),
    ] {
        if let Some(chart) = sparse_chart(title, label, values, series) {
            figure(page, title, chart, &[]);
        }
    }
}

/// A histogram of sparse bins, trimmed to the central 99.8% of reads.
fn sparse_chart(title: &str, label: &str, bins: &[(f64, u64)], series: usize) -> Option<Chart> {
    let total: u64 = bins.iter().map(|b| b.1).sum();
    if total == 0 || bins.len() < 2 {
        return None;
    }
    let width = bins
        .windows(2)
        .map(|w| w[1].0 - w[0].0)
        .fold(f64::INFINITY, f64::min);
    let mut cumulative = 0;
    let lo = total / 1000;
    let hi = total - lo;
    let kept: Vec<&(f64, u64)> = bins
        .iter()
        .filter(|b| {
            let before = cumulative;
            cumulative += b.1;
            cumulative > lo && before < hi
        })
        .collect();
    let (first, last) = (kept.first()?.0, kept.last()?.0 + width);
    let top = kept
        .iter()
        .map(|b| 100.0 * b.1 as f64 / total as f64)
        .fold(0.0, f64::max)
        * 1.08;
    let bars = kept
        .iter()
        .map(|&&(start, count)| {
            let p = 100.0 * count as f64 / total as f64;
            Bar {
                x0: start,
                x1: start + width,
                y: p,
                tip: format!("{start:.2}–{:.2}: {p:.2}% of reads", start + width),
            }
        })
        .collect();
    Some(
        Chart::new(
            title,
            Axis::linear(label, first, last, Units::Plain),
            Axis::linear("% of reads", 0.0, top, Units::Percent),
        )
        .layer(Layer::Bars { series, bars }),
    )
}

pub(super) fn kinnex(page: &mut Page, group: &Group) {
    if let Some(chart) = histogram_chart(
        "Kinnex adapters per read",
        "segmentation adapters found",
        &group.data.kinnex,
        5,
    ) {
        figure(page, "Reads by Kinnex adapters found", chart, &[]);
    }
}

/// Run information from BAM headers: flow cells, devices, models, kits and programs.
pub(super) fn run_info(page: &mut Page, report: &Report) {
    let runs = &report.runs;
    if runs.read_groups.is_empty() && runs.programs.is_empty() {
        return;
    }
    page.push_str(r#"<div class="panel"><h3>Runs</h3>"#);
    if !runs.read_groups.is_empty() {
        page.push_str(long_table_wrap(runs.read_groups.len()));
        page.push_str(r#"<table><tr><th class="l">Read group</th><th class="l">Platform</th><th class="l">Sample</th><th class="l">Flow cell / movie</th><th class="l">Device</th><th class="l">Started</th><th class="l">Basecall model</th><th class="l">Modification models</th><th class="l">PacBio</th></tr>"#);
        let cell = |v: &Option<String>| v.as_deref().map_or_else(String::new, escape);
        for g in &runs.read_groups {
            let _ = write!(
                page,
                r#"<tr><td class="l"><code>{}</code></td><td class="l">{}</td><td class="l">{}</td><td class="l">{}</td><td class="l">{}</td><td class="l">{}</td><td class="l">{}</td><td class="l">{}</td><td class="l">{}</td></tr>"#,
                escape(&g.id),
                cell(&g.platform),
                cell(&g.sample),
                cell(&g.flow_cell),
                cell(&g.device),
                cell(&g.started),
                cell(&g.basecall_model),
                cell(&g.modbase_models),
                cell(&g.pacbio)
            );
        }
        page.push_str("</table></div>");
    }
    if !runs.programs.is_empty() {
        let programs: Vec<String> = runs
            .programs
            .iter()
            .map(|(name, version)| format!("{} {}", escape(name), escape(version)))
            .collect();
        let _ = write!(
            page,
            r#"<p class="explain">Programs: {}</p>"#,
            programs.join(" · ")
        );
    }
    page.push_str("</div>");
}
