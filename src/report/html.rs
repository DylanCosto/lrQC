//! The self-contained HTML report: one section per group, an overview comparing groups,
//! charts as inline SVG, and navigation by CSS `:target` alone.

use std::collections::BTreeMap;
use std::fmt::Write as _;

use super::chart::{
    self, escape, si, Axis, Bar, Cell, Chart, Detail, Hit, Layer, LineStyle, Units, PALETTE,
};
use super::checks::{ModuleId, Status};
use super::model::{Group, PositionContent, PositionQuality, ReadEnd};
use super::{GroupKind, Report};

mod extras;
use crate::bins;
use crate::read::Platform;
use crate::stats::{READ_Q_LEVELS, SPEED_BIN_WIDTH};

const CSS: &str = include_str!("../../assets/report.css");

/// Reports with more groups than this draw compact charts and summarize the barcodes in
/// the overview instead of drawing each one.
const MANY_GROUPS: usize = 12;

/// The page being written, with the chart detail its size calls for.
///
/// Charts are held back until other markup follows them, so that each run of charts can
/// be laid out as one grid: two to a row, the odd one out spanning the row.
struct Page {
    html: String,
    detail: Detail,
    figures: Vec<Figure>,
}

/// A chart waiting to be placed in the chart grid.
struct Figure {
    caption: String,
    chart: Chart,
    /// Legend markup, or empty.
    legend: String,
    span: Span,
}

/// How much of the chart grid's row a figure takes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Span {
    /// One column, or the whole row when the figure is the odd one out.
    Column,
    Row,
}

impl Page {
    fn new(groups: usize) -> Self {
        Self {
            html: String::with_capacity(1 << 20),
            detail: if groups > MANY_GROUPS {
                Detail::Compact
            } else {
                Detail::Full
            },
            figures: Vec::new(),
        }
    }

    fn push_str(&mut self, text: &str) {
        self.place_figures();
        self.html.push_str(text);
    }

    /// Writes the charts held back as one grid.
    fn place_figures(&mut self) {
        if self.figures.is_empty() {
            return;
        }
        let mut figures = std::mem::take(&mut self.figures);
        // The last of an odd run of single-column figures takes the whole row.
        let mut run = 0;
        for at in 0..=figures.len() {
            match figures.get(at).map(|f| f.span) {
                Some(Span::Column) => run += 1,
                Some(Span::Row) | None => {
                    if run % 2 == 1 {
                        figures[at - 1].span = Span::Row;
                    }
                    run = 0;
                }
            }
        }
        self.html.push_str(r#"<div class="charts">"#);
        for figure in figures {
            let (class, chart) = match figure.span {
                Span::Column => ("", figure.chart),
                Span::Row => (r#" class="wide""#, figure.chart.wide()),
            };
            let _ = write!(
                self.html,
                "<figure{class}><figcaption>{}</figcaption>{}{}</figure>",
                escape(&figure.caption),
                chart.render(self.detail),
                figure.legend
            );
        }
        self.html.push_str("</div>");
    }

    /// The finished page.
    fn finish(mut self) -> String {
        self.place_figures();
        self.html
    }
}

impl std::fmt::Write for Page {
    fn write_str(&mut self, text: &str) -> std::fmt::Result {
        self.push_str(text);
        Ok(())
    }
}

/// Renders the whole report page.
pub fn render(report: &Report) -> String {
    let mut page = Page::new(report.groups.len());
    let _ = write!(
        page,
        r#"<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">{refresh}<title>lrQC: {title}</title><style>{CSS}</style></head><body>"#,
        title = escape(&report.title),
        refresh = report.refresh_secs.map_or_else(String::new, |secs| format!(
            r#"<meta http-equiv="refresh" content="{secs}">"#
        )),
    );
    let _ = write!(
        page,
        r#"<header class="top"><h1>{}</h1><span class="meta">lrQC {}</span></header><div class="layout">"#,
        escape(&report.title),
        report.version
    );
    nav(&mut page, report);
    page.push_str("<main>");
    for (index, group) in report.groups.iter().enumerate() {
        group_section(&mut page, report, index, group);
    }
    overview(&mut page, report);
    page.push_str("</main></div></body></html>");
    page.finish()
}

/// One page comparing the groups of several reports.
pub fn render_comparison(reports: &[Report]) -> String {
    let mut page = Page::new(reports.iter().map(|r| r.groups.len()).sum());
    let title = format!("lrQC: {} reports compared", reports.len());
    let titles: Vec<&str> = reports.iter().map(|r| r.title.as_str()).collect();
    let _ = write!(
        page,
        r#"<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>{title}</title><style>{CSS}</style></head><body><header class="top"><h1>{title}</h1><span class="meta">{}</span></header><main><section class="group" id="overview">"#,
        escape(&titles.join(" · "))
    );
    let rows: Vec<Row<'_>> = reports
        .iter()
        .flat_map(|r| {
            r.groups.iter().map(move |g| Row {
                label: format!("{} · {}", r.title, display_name(g)),
                link: None,
                group: g,
            })
        })
        .collect();
    status_matrix(&mut page, &rows);
    summary_table(&mut page, &rows);
    if let Some(first) = reports.first() {
        let totals: Vec<&Group> = reports.iter().map(|r| &r.groups[0]).collect();
        comparisons(
            &mut page,
            first,
            &totals,
            &titles.iter().map(|t| (*t).to_owned()).collect::<Vec<_>>(),
            "All reads, by report",
        );
        let barcode_rows: Vec<&Row<'_>> = rows
            .iter()
            .filter(|r| r.group.kind != GroupKind::AllReads)
            .collect();
        if barcode_rows.len() > 1 {
            let groups: Vec<&Group> = barcode_rows.iter().map(|r| r.group).collect();
            let labels: Vec<String> = barcode_rows.iter().map(|r| r.label.clone()).collect();
            comparisons(&mut page, first, &groups, &labels, "Barcodes, by report");
        }
    }
    page.push_str("</section></main></body></html>");
    page.finish()
}

/// A group's name with its sample name, if the sample sheet gives one.
fn display_name(group: &Group) -> String {
    match &group.sample {
        Some(sample) => format!("{} · {sample}", group.name),
        None => group.name.clone(),
    }
}

fn group_id(index: usize, group: &Group) -> String {
    let slug: String = group
        .name
        .chars()
        .map(|c| {
            if c.is_ascii_alphanumeric() || c == '-' || c == '_' {
                c
            } else {
                '-'
            }
        })
        .collect();
    format!("g{index}-{slug}")
}

/// Most severe verdict of a group, for its navigation dot.
fn worst(group: &Group) -> Status {
    group
        .modules
        .iter()
        .map(|m| m.status)
        .filter(|s| matches!(s, Status::Pass | Status::Warn | Status::Fail))
        .max()
        .unwrap_or(Status::Info)
}

fn status_class(status: Status) -> &'static str {
    match status {
        Status::Pass => "pass",
        Status::Warn => "warn",
        Status::Fail => "fail",
        Status::Info => "info",
        Status::Absent => "absent",
    }
}

fn nav(page: &mut Page, report: &Report) {
    page.push_str(r##"<nav class="side"><a href="#overview">Overview</a><h3>Groups</h3>"##);
    for (index, group) in report.groups.iter().enumerate() {
        let _ = write!(
            page,
            r##"<a href="#{}"><span class="dot {}"></span>{}<span class="n">{}</span></a>"##,
            group_id(index, group),
            status_class(worst(group)),
            escape(&display_name(group)),
            si(group.summary.reads as f64)
        );
    }
    page.push_str("</nav>");
}

fn thousands(n: u64) -> String {
    let digits = n.to_string();
    let mut out = String::with_capacity(digits.len() + digits.len() / 3);
    for (i, c) in digits.chars().enumerate() {
        if i > 0 && (digits.len() - i).is_multiple_of(3) {
            out.push(',');
        }
        out.push(c);
    }
    out
}

fn opt(value: Option<f64>, digits: usize, suffix: &str) -> String {
    value.map_or_else(|| "–".to_owned(), |v| format!("{v:.digits$}{suffix}"))
}

fn tile(page: &mut Page, value: &str, key: &str) {
    let _ = write!(
        page,
        r#"<div class="tile"><div class="v">{value}</div><div class="k">{key}</div></div>"#
    );
}

/// Queues a chart for the chart grid, with a legend of its series.
fn figure(page: &mut Page, caption: &str, chart: Chart, legend: &[(String, usize, LineStyle)]) {
    let legend = if legend.is_empty() {
        String::new()
    } else {
        let items: String = legend
            .iter()
            .map(|(name, series, style)| legend_item(name, *series, *style))
            .collect();
        format!(r#"<div class="legend">{items}</div>"#)
    };
    page.figures.push(Figure {
        caption: caption.to_owned(),
        chart,
        legend,
        span: Span::Column,
    });
}

fn legend_item(name: &str, series: usize, style: LineStyle) -> String {
    let class = match style {
        LineStyle::Dashed => "dashed",
        LineStyle::Solid | LineStyle::Thin => "",
    };
    format!(
        r#"<span class="s{}"><i class="{class}"></i>{}</span>"#,
        series % PALETTE,
        escape(name)
    )
}

/// What each module shows and how to read it.
fn explanation(id: ModuleId) -> &'static str {
    match id {
        ModuleId::BasicStatistics => "Read counts, yield and length and quality summaries. Read Q is the basecaller's (ONT qs tag or summary column), PacBio's predicted accuracy (rq), or else the mean error probability of the bases. Base Q is the mean error probability of all bases, as a Phred score.",
        ModuleId::ReadLength => "Reads and bases by read length, on a log scale. The bases curve shows where the yield is; a peak of short reads often means degraded DNA or short fragments.",
        ModuleId::ReadQuality => "Reads by read Q. Shaded bands mark the fail (red), warn (amber) and pass (green) ranges.",
        ModuleId::LengthVsQuality => "Reads by length and read Q together; color is the read count on a log scale. Long, low-quality reads or short, high-quality ones stand out as separate clouds.",
        ModuleId::BaseQuality => "Base quality by position, from the start (5') and from the end (3') of reads: median line, mean dashed, 25-75% dark band and 10-90% light band. Positions are single bases up to 100, then wider bins; bins reached by fewer than 0.5% of reads are left out.",
        ModuleId::BaseContent => "Percent of A, C, G and T by position from each end. Lines should run flat and parallel; divergence at the ends shows adapters, barcodes or primers left on the reads.",
        ModuleId::GcContent => "Reads by GC percent, with a normal curve of the same mean and spread for reference. Long reads average GC over tens of kilobases, so a skewed single peak is normal; a separate second peak means contamination or a mixed sample, and is what this check grades.",
        ModuleId::NContent => "Percent of N (or other ambiguous) bases by position from each end.",
        ModuleId::Duplication => "Duplication of a reproducible sample of reads (those with the smallest name hashes), counting reads as duplicates when their first and last 50 bases both match. Long genomic reads are almost never duplicates; amplicons and PCR libraries are.",
        ModuleId::Overrepresented => "First and last 50 bases seen in more than 0.1% of sampled reads, with the adapter or primer they match. Recurring read starts or ends are untrimmed adapters, primers, barcodes or a dominant amplicon.",
        ModuleId::Adapters => "Reads holding each known adapter or primer, found by shared 12-mers on either strand: within 150 bases of the start or the end, or inside the read. Adapters inside reads mark chimeras (ONT ligation) or missed adapter splits. Curves show the cumulative percent of reads with the adapter at or before each position.",
        ModuleId::ReadIds => "Reads whose ID appears more than once, as when a file is given twice or overlapping inputs are combined.",
        ModuleId::BarcodeBalance => "Reads of this barcode against the median barcode's, and the unclassified share of the run.",
        ModuleId::YieldOverTime => "Bases sequenced over the run (from read start times), read Q over time and the share of passing reads. Read Q dropping or yield flattening early points to pore blocking, low library or the need for a flush.",
        ModuleId::Channels => "Reads per channel, laid out in channel number order (not the flow cell's physical layout), and the spread of reads per channel. Many quiet channels mean few active pores.",
        ModuleId::TranslocationSpeed => "Read length over read duration, in bases per second. R10.4.1 kits run near 400 b/s (5 kHz) or 260 b/s; a slow or broad peak suggests motor protein or temperature trouble.",
        ModuleId::EndReasons => "Why the device ended each read, from dorado's er tags or sequencing summaries. Many unblock or mux-change ends point to pore blocking; data_service_unblock ends are adaptive sampling rejections.",
        ModuleId::PacBio => "PacBio HiFi read properties: CCS passes (np), effective coverage (ec), lima barcode score (bq) and kinetics tags.",
        ModuleId::Modifications => "Modified-base calls from MM/ML tags, by modification code: calls per read, percent called modified (probability 0.5 or more), the same for calls on a CpG of the read, and the call probability spread. Bases an implicit-mode entry (C+m., as dorado and jasmine write) skips count as unmodified calls, as the SAM specification has it. Models that call every C (5mC_5hmC) give low all-context rates; human CpG 5mC runs about 70-80% in blood and tissues, and often lower in cell lines. A peak near 0.5 means uncertain calls.",
        ModuleId::FilterYield => "Reads, bases and (with --genome-size) coverage left after filtering by minimum read length and minimum read Q, counted exactly from every read. Use it to choose filters before assembly or variant calling.",
        ModuleId::Contaminants => "Reads matched to 14 built-in references by sampled 19-mers: spike-in controls (lambda DCS, yeast ENO2, ERCC, SIRV), the human mitochondrion, cloning vectors (UniVec), and contaminants (phiX174, E. coli K-12, five Mycoplasma and Acholeplasma species, the human nuclear genome). A read goes to the reference holding the largest share of its sampled k-mers, when that is at least 3 k-mers and at least half the share a read of its Q keeps (never below 8%). Only contaminants are graded: a contaminant matching most reads is taken to be the sample, and the others are still graded. Human is judged on the reads long enough to test (about 7 kb), or as a lower bound over all reads, which only warns.",
        ModuleId::BarcodeCrosstalk => "Barcode sequences found inside reads (every ONT kit dorado knows, matched within their flanks; PacBio's from --barcode-names, or PacBio's built-in indexes for the barcodes the file names or sample sheet name). A read's own barcode at an end means barcodes were not trimmed; another sample's barcode means cross-talk, misassignment or chimeras joining two samples' molecules.",
        ModuleId::Transcripts => "cDNA structure from ONT's strand-switching (SSP) and VN primers: full-length reads have one primer at each end, and which comes first gives the strand read. Poly(A) lengths come from dorado's estimate (pt tag, --estimate-poly-a).",
        ModuleId::Duplex => "dorado duplex classes (dx tag): duplex reads, simplex reads whose strand pair made a duplex read, and other simplex reads, with their yield and read Q. The duplex rate counts each duplex base twice (it covers both strands) against all simplex bases.",
        ModuleId::Signal => "Raw-signal properties dorado records: signal shift and scale (sm, sd), samples per base (ns over length, about 12.5 at 5 kHz and 400 b/s), yield by mux (mx) and reads split from a longer read (pi). Broad or shifted signal distributions point to pore or chemistry problems.",
        ModuleId::Kinnex => "Kinnex/MAS segmentation adapters found per read (--kinnex-adapters, skera's FASTA). A full array holds every adapter of the kit; reads with fewer are partial arrays, which give fewer segmented reads.",
        ModuleId::Alignment => "Mapping of aligned BAMs: primary alignments mapped, unmapped, secondary and supplementary records, reads split across loci (SA tag), MAPQ, identity (BLAST-like, mismatches and indels over alignment columns; gap-compressed counts each indel once) and the fraction of each read aligned.",
    }
}

fn group_section(page: &mut Page, report: &Report, index: usize, group: &Group) {
    let id = group_id(index, group);
    let _ = write!(
        page,
        r#"<section class="group" id="{id}"><h2>{}</h2>"#,
        escape(&display_name(group))
    );
    for module in &group.modules {
        if module.status == Status::Absent {
            continue;
        }
        let _ = write!(
            page,
            r#"<article class="module" id="{id}-{}"><h3>{}</h3>"#,
            module.id.key(),
            module.id.title()
        );
        if !module.note.is_empty() {
            let _ = write!(page, r#"<p class="note">{}</p>"#, escape(&module.note));
        }
        let _ = write!(page, r#"<p class="explain">{}</p>"#, explanation(module.id));
        module_body(page, report, group, module.id);
        page.push_str("</article>");
    }
    let absent: Vec<&str> = group
        .modules
        .iter()
        .filter(|m| m.status == Status::Absent)
        .map(|m| m.id.title())
        .collect();
    if !absent.is_empty() {
        let _ = write!(
            page,
            r#"<p class="explain">Not available for this input: {}.</p>"#,
            escape(&absent.join(", "))
        );
    }
    page.push_str("</section>");
}

fn module_body(page: &mut Page, report: &Report, group: &Group, id: ModuleId) {
    match id {
        ModuleId::BasicStatistics => basic_table(page, group),
        ModuleId::ReadLength => read_length(page, group),
        ModuleId::ReadQuality => read_quality(page, group),
        ModuleId::LengthVsQuality => length_vs_quality(page, group),
        ModuleId::BaseQuality => base_quality(page, report, group),
        ModuleId::BaseContent => base_content(page, report, group),
        ModuleId::GcContent => gc_content(page, group),
        ModuleId::NContent => n_content(page, report, group),
        ModuleId::Duplication => duplication(page, group),
        ModuleId::Overrepresented => overrepresented(page, group),
        ModuleId::Adapters => adapters(page, report, group),
        ModuleId::ReadIds | ModuleId::BarcodeBalance => {}
        ModuleId::YieldOverTime => yield_over_time(page, group),
        ModuleId::Channels => channels(page, group),
        ModuleId::TranslocationSpeed => speed(page, group),
        ModuleId::EndReasons => end_reasons(page, group),
        ModuleId::PacBio => pacbio(page, group),
        ModuleId::Modifications => modifications(page, group),
        ModuleId::Alignment => alignment(page, group),
        ModuleId::FilterYield => extras::filter_yield(page, report, group),
        ModuleId::Contaminants => extras::contaminants(page, group),
        ModuleId::BarcodeCrosstalk => extras::barcodes_inside(page, group),
        ModuleId::Transcripts => extras::transcripts(page, group),
        ModuleId::Duplex => extras::duplex(page, group),
        ModuleId::Signal => extras::signal(page, group),
        ModuleId::Kinnex => extras::kinnex(page, group),
    }
}

fn basic_table(page: &mut Page, group: &Group) {
    let s = &group.summary;
    let rows: Vec<(&str, String)> = vec![
        ("Platform", group.platform.label().to_owned()),
        ("Reads", thousands(s.reads)),
        ("Bases", thousands(s.bases)),
        ("Share of run reads", format!("{:.2}%", s.share_pct)),
        ("Mean length", format!("{:.0}", s.mean_length)),
        ("Median length", thousands(s.median_length.round() as u64)),
        ("Read N50", thousands(s.n50)),
        ("Read N90", thousands(s.n90)),
        (
            "Shortest / longest",
            format!("{} / {}", thousands(s.shortest), thousands(s.longest)),
        ),
        (
            "Median / mean read Q",
            format!(
                "{} / {}",
                opt(s.median_read_q, 2, ""),
                opt(s.mean_read_q, 2, "")
            ),
        ),
        (
            "Reads ≥ Q10 / Q20 / Q30",
            format!(
                "{} / {} / {}",
                opt(s.reads_q10_pct, 1, "%"),
                opt(s.reads_q20_pct, 1, "%"),
                opt(s.reads_q30_pct, 1, "%")
            ),
        ),
        (
            "Bases in reads ≥ Q20",
            s.bases_q20.map_or_else(|| "–".to_owned(), thousands),
        ),
        ("Base Q", opt(s.base_q, 2, "")),
        ("GC", opt(s.gc_pct, 2, "%")),
        ("Duplicate read IDs", thousands(s.duplicate_names)),
    ];
    page.push_str(r#"<dl class="kv">"#);
    for (key, value) in rows {
        let _ = write!(page, "<div><dt>{key}</dt><dd>{value}</dd></div>");
    }
    page.push_str("</dl>");
}

/// Fail, warn and pass bands for a chart with Q along y, up to `top`.
fn quality_zones(group: &Group, top: f64) -> Layer {
    Layer::Zones(quality_bands(group, top))
}

/// The same bands for a chart with read Q along x, up to `right`.
fn quality_zones_x(group: &Group, right: f64) -> Layer {
    Layer::XZones(quality_bands(group, right))
}

fn quality_bands(group: &Group, end: f64) -> Vec<(f64, f64, &'static str)> {
    // The levels the checks used (--min-qscore, --thresholds); older reports lack them.
    let (warn, fail) = group.quality_levels.unwrap_or(match group.platform {
        Platform::PacBio => (20.0, 15.0),
        Platform::Ont | Platform::Unknown => (10.0, 7.5),
    });
    vec![
        (0.0, fail, "zone-bad"),
        (fail, warn, "zone-warn"),
        (warn, end, "zone-good"),
    ]
}

fn read_length(page: &mut Page, group: &Group) {
    let bins = &group.data.length_bins;
    if bins.is_empty() {
        return;
    }
    let reads: u64 = bins.iter().map(|b| b.reads).sum();
    let bases: u64 = bins.iter().map(|b| b.bases).sum();
    let x = Axis::log(
        "read length (bases)",
        bins[0].start,
        bins[bins.len() - 1].end,
        Units::Si,
    );
    for (title, total, value, series) in [
        ("Reads by length", reads, 0usize, 0usize),
        ("Bases by length", bases, 1, 1),
    ] {
        let values: Vec<f64> = bins
            .iter()
            .map(|b| {
                100.0 * if value == 0 { b.reads } else { b.bases } as f64 / total.max(1) as f64
            })
            .collect();
        let top = values.iter().copied().fold(0.0, f64::max) * 1.05;
        let bars = bins
            .iter()
            .zip(&values)
            .map(|(b, &v)| Bar {
                x0: b.start,
                x1: b.end,
                y: v,
                tip: format!(
                    "{}–{} bases: {} reads, {} bases ({v:.2}%)",
                    si(b.start),
                    si(b.end),
                    thousands(b.reads),
                    si(b.bases as f64)
                ),
            })
            .collect();
        let chart = Chart::new(
            title,
            x.clone(),
            Axis::linear(
                &format!("% of {}", if value == 0 { "reads" } else { "bases" }),
                0.0,
                top,
                Units::Percent,
            ),
        )
        .layer(Layer::Bars { series, bars });
        figure(page, title, chart, &[]);
    }
}

fn read_quality(page: &mut Page, group: &Group) {
    let bins = &group.data.read_q;
    if bins.is_empty() {
        return;
    }
    // Half-Q bins read better than the tenth-Q bins the data keeps.
    let mut halves: Vec<(f64, u64, u64)> = Vec::new();
    for b in bins {
        let start = (b.q * 2.0).floor() / 2.0;
        match halves.last_mut() {
            Some(last) if (last.0 - start).abs() < 1e-9 => {
                last.1 += b.reads;
                last.2 += b.bases;
            }
            _ => halves.push((start, b.reads, b.bases)),
        }
    }
    let reads: u64 = halves.iter().map(|h| h.1).sum();
    let bases: u64 = halves.iter().map(|h| h.2).sum();
    let max_q = halves.last().map_or(40.0, |h| h.0 + 0.5).max(15.0);
    let x = Axis::linear("read Q", 0.0, max_q.ceil(), Units::Plain);
    let read_pct: Vec<f64> = halves
        .iter()
        .map(|h| 100.0 * h.1 as f64 / reads.max(1) as f64)
        .collect();
    let base_pct: Vec<f64> = halves
        .iter()
        .map(|h| 100.0 * h.2 as f64 / bases.max(1) as f64)
        .collect();
    let top = read_pct
        .iter()
        .chain(&base_pct)
        .copied()
        .fold(0.0, f64::max)
        * 1.05;
    let bars = halves
        .iter()
        .zip(&read_pct)
        .map(|(h, &p)| Bar {
            x0: h.0,
            x1: h.0 + 0.5,
            y: p,
            tip: format!(
                "Q{:.1}–{:.1}: {} reads ({p:.2}%)",
                h.0,
                h.0 + 0.5,
                thousands(h.1)
            ),
        })
        .collect();
    let line: Vec<(f64, f64)> = halves
        .iter()
        .zip(&base_pct)
        .map(|(h, &p)| (h.0 + 0.25, p))
        .collect();
    let chart = Chart::new(
        "Read quality",
        x,
        Axis::linear("% of reads / bases", 0.0, top, Units::Percent),
    )
    .layer(quality_zones_x(group, max_q.ceil()))
    .layer(Layer::Bars { series: 0, bars })
    .layer(Layer::Line {
        name: "% of bases".into(),
        series: 1,
        points: line,
        style: LineStyle::Solid,
    });
    figure(
        page,
        "Reads (bars) and bases (line) by read Q",
        chart,
        &[
            ("reads".into(), 0, LineStyle::Solid),
            ("bases".into(), 1, LineStyle::Solid),
        ],
    );
}

/// Sums `(column, row, count)` cells into blocks of `merge.0` columns by `merge.1` rows,
/// giving each block's first column and row with its total.
fn merge_grid(cells: &[(usize, usize, u64)], merge: (usize, usize)) -> Vec<(usize, usize, u64)> {
    let mut blocks: BTreeMap<(usize, usize), u64> = BTreeMap::new();
    for &(column, row, count) in cells {
        let block = (column / merge.0 * merge.0, row / merge.1 * merge.1);
        *blocks.entry(block).or_default() += count;
    }
    blocks
        .into_iter()
        .map(|((column, row), count)| (column, row, count))
        .collect()
}

/// "Q18", or "Q18–19" for a block of several Q values.
fn q_range(q: usize, span: usize) -> String {
    if span > 1 {
        format!("Q{q}–{}", q + span - 1)
    } else {
        format!("Q{q}")
    }
}

/// Heatmap layers for a grid of read counts: cells colored by log count, and hover
/// rectangles that may span several cells.
///
/// `x_edges` gives a column's lower and upper x for a block starting at a column and
/// spanning some columns; `tip` names a block from its column, row and spans.
fn count_heatmap(
    shades: u32,
    counts: &[(usize, usize, u64)],
    cell_merge: (usize, usize),
    hit_merge: (usize, usize),
    x_edges: &dyn Fn(usize, usize) -> (f64, f64),
    tip: &dyn Fn(usize, usize, (usize, usize)) -> String,
) -> [Layer; 2] {
    let drawn = merge_grid(counts, cell_merge);
    let max = drawn.iter().map(|c| c.2).max().unwrap_or(1) as f64;
    let cells = drawn
        .iter()
        .map(|&(column, row, count)| {
            let (x0, x1) = x_edges(column, cell_merge.0);
            Cell {
                x0,
                x1,
                y0: row as f64 - 0.5,
                y1: (row + cell_merge.1) as f64 - 0.5,
                color: chart::ramp((count as f64).ln_1p() / max.ln_1p(), shades),
            }
        })
        .collect();
    let hits = merge_grid(counts, hit_merge)
        .into_iter()
        .map(|(column, row, count)| {
            let (x0, x1) = x_edges(column, hit_merge.0);
            Hit {
                x0,
                x1,
                y0: row as f64 - 0.5,
                y1: (row + hit_merge.1) as f64 - 0.5,
                tip: format!(
                    "{}: {} reads",
                    tip(column, row, hit_merge),
                    thousands(count)
                ),
            }
        })
        .collect();
    [Layer::Cells(cells), Layer::Hits(hits)]
}

fn length_vs_quality(page: &mut Page, group: &Group) {
    let counts = &group.data.length_q;
    if counts.is_empty() {
        return;
    }
    let first_len = counts.iter().map(|c| c.0).min().unwrap_or(0);
    let last_len = counts.iter().map(|c| c.0).max().unwrap_or(0);
    let max_q = counts.iter().map(|c| c.1).max().unwrap_or(30).max(15) + 1;
    // Hover blocks span several cells; compact pages also draw coarser cells.
    let (cell_merge, hit_merge) = match page.detail {
        Detail::Full => ((1, 1), (2, 2)),
        Detail::Compact => ((2, 2), (6, 4)),
    };
    let x_edges = |bin: usize, span: usize| {
        (
            bins::length_bin_start(bin),
            bins::length_bin_start(bin + span),
        )
    };
    let [cells, hits] = count_heatmap(
        page.detail.shades(),
        counts,
        cell_merge,
        hit_merge,
        &x_edges,
        &|bin, q, (bins, qs)| {
            let (start, end) = x_edges(bin, bins);
            format!("{}–{} bases, {}", si(start), si(end), q_range(q, qs))
        },
    );
    let chart = Chart::new(
        "Read length vs quality",
        Axis::log(
            "read length (bases)",
            bins::length_bin_start(first_len),
            bins::length_bin_start(last_len + 1),
            Units::Si,
        ),
        Axis::linear(
            "read Q",
            0.0,
            (max_q.min(READ_Q_LEVELS)) as f64,
            Units::Plain,
        ),
    )
    .layer(cells)
    .layer(hits)
    .tall(340.0);
    figure(
        page,
        "Reads by length and read Q (color: log read count, dark to bright)",
        chart,
        &[],
    );
}

/// The x value at which a position bin is drawn: its middle, counting from 1.
fn bin_x(report: &Report, bin: usize) -> f64 {
    let (start, end) = report.position_bins[bin];
    match end {
        Some(end) if end - start > 1 => (start + end) as f64 / 2.0 + 0.5,
        _ => start as f64 + 1.0,
    }
}

fn position_axis(report: &Report, last_bin: usize, end: ReadEnd) -> Axis {
    let label = match end {
        ReadEnd::Start => "position from read start (bases)",
        ReadEnd::End => "position from read end (bases)",
    };
    Axis::log(
        label,
        1.0,
        bin_x(report, last_bin).max(10.0) * 1.1,
        Units::Si,
    )
}

fn base_quality(page: &mut Page, report: &Report, group: &Group) {
    for (end, rows) in [
        (ReadEnd::Start, &group.data.quality_from_start),
        (ReadEnd::End, &group.data.quality_from_end),
    ] {
        let Some(last) = rows.last() else { continue };
        let top = rows.iter().map(|r| r.p90).fold(0.0, f64::max).max(20.0) + 2.0;
        let pts = |f: fn(&PositionQuality) -> f64| -> Vec<(f64, f64)> {
            rows.iter().map(|r| (bin_x(report, r.bin), f(r))).collect()
        };
        let chart = Chart::new(
            "Base quality",
            position_axis(report, last.bin, end),
            Axis::linear("base Q", 0.0, top, Units::Plain),
        )
        .layer(quality_zones(group, top))
        .layer(Layer::Band {
            series: 0,
            points: rows
                .iter()
                .map(|r| (bin_x(report, r.bin), r.p10, r.p90))
                .collect(),
            opacity: 0.15,
        })
        .layer(Layer::Band {
            series: 0,
            points: rows
                .iter()
                .map(|r| (bin_x(report, r.bin), r.p25, r.p75))
                .collect(),
            opacity: 0.3,
        })
        .layer(Layer::Line {
            name: "median".into(),
            series: 0,
            points: pts(|r| r.median),
            style: LineStyle::Solid,
        })
        .layer(Layer::Line {
            name: "mean".into(),
            series: 1,
            points: pts(|r| r.mean),
            style: LineStyle::Dashed,
        });
        let caption = match end {
            ReadEnd::Start => "From the read start (5')",
            ReadEnd::End => "From the read end (3')",
        };
        figure(
            page,
            caption,
            chart,
            &[
                ("median".into(), 0, LineStyle::Solid),
                ("mean".into(), 1, LineStyle::Dashed),
            ],
        );
    }
    // Overall base quality distribution.
    let quals = &group.data.base_quals;
    let total: u64 = quals.iter().sum();
    if total > 0 {
        let values: Vec<f64> = quals
            .iter()
            .map(|&c| 100.0 * c as f64 / total as f64)
            .collect();
        let top = values.iter().copied().fold(0.0, f64::max) * 1.05;
        let bars = values
            .iter()
            .enumerate()
            .map(|(q, &p)| Bar {
                x0: q as f64 - 0.5,
                x1: q as f64 + 0.5,
                y: p,
                tip: format!("Q{q}: {p:.2}% of bases"),
            })
            .collect();
        let chart = Chart::new(
            "Base qualities",
            Axis::linear("base Q", -0.5, quals.len() as f64 - 0.5, Units::Plain),
            Axis::linear("% of bases", 0.0, top, Units::Percent),
        )
        .layer(Layer::Bars { series: 2, bars });
        figure(page, "All bases by quality", chart, &[]);
    }
}

fn content_chart(
    report: &Report,
    rows: &[PositionContent],
    end: ReadEnd,
    with_bases: bool,
) -> Option<Chart> {
    let last = rows.last()?;
    let series = |f: fn(&PositionContent) -> f64| -> Vec<(f64, f64)> {
        rows.iter().map(|r| (bin_x(report, r.bin), f(r))).collect()
    };
    let mut chart = if with_bases {
        let top = rows
            .iter()
            .map(|r| r.a.max(r.c).max(r.g).max(r.t))
            .fold(0.0, f64::max)
            .max(40.0)
            + 5.0;
        Chart::new(
            "Base content",
            position_axis(report, last.bin, end),
            Axis::linear("% of bases", 0.0, top.min(100.0), Units::Percent),
        )
    } else {
        let top = rows.iter().map(|r| r.n).fold(0.0, f64::max).max(1.0) * 1.2;
        Chart::new(
            "N content",
            position_axis(report, last.bin, end),
            Axis::linear("% N", 0.0, top, Units::Percent),
        )
    };
    if with_bases {
        for (name, i, f) in [
            (
                "A",
                2usize,
                (|r: &PositionContent| r.a) as fn(&PositionContent) -> f64,
            ),
            ("C", 0, |r| r.c),
            ("G", 4, |r| r.g),
            ("T", 1, |r| r.t),
        ] {
            chart = chart.layer(Layer::Line {
                name: name.into(),
                series: i,
                points: series(f),
                style: LineStyle::Solid,
            });
        }
    } else {
        chart = chart.layer(Layer::Line {
            name: "N".into(),
            series: 11,
            points: series(|r| r.n),
            style: LineStyle::Solid,
        });
    }
    Some(chart)
}

fn base_content(page: &mut Page, report: &Report, group: &Group) {
    let legend = [
        ("A".to_owned(), 2, LineStyle::Solid),
        ("C".to_owned(), 0, LineStyle::Solid),
        ("G".to_owned(), 4, LineStyle::Solid),
        ("T".to_owned(), 1, LineStyle::Solid),
    ];
    for (end, rows, caption) in [
        (
            ReadEnd::Start,
            &group.data.content_from_start,
            "From the read start (5')",
        ),
        (
            ReadEnd::End,
            &group.data.content_from_end,
            "From the read end (3')",
        ),
    ] {
        if let Some(chart) = content_chart(report, rows, end, true) {
            figure(page, caption, chart, &legend);
        }
    }
}

fn n_content(page: &mut Page, report: &Report, group: &Group) {
    // Modern basecallers (dorado, ccs) never write N: a line along the zero axis would look
    // like an empty chart, so say so instead.
    let rows = group
        .data
        .content_from_start
        .iter()
        .chain(&group.data.content_from_end);
    if rows.clone().next().is_some() && rows.clone().all(|c| c.n == 0.0) {
        page.push_str(
            r#"<p class="note">No read has an N (or other ambiguous) base at any position.</p>"#,
        );
        return;
    }
    for (end, rows, caption) in [
        (
            ReadEnd::Start,
            &group.data.content_from_start,
            "From the read start (5')",
        ),
        (
            ReadEnd::End,
            &group.data.content_from_end,
            "From the read end (3')",
        ),
    ] {
        if let Some(chart) = content_chart(report, rows, end, false) {
            figure(page, caption, chart, &[]);
        }
    }
}

fn gc_content(page: &mut Page, group: &Group) {
    let total: u64 = group.data.gc.iter().sum();
    if total == 0 {
        return;
    }
    let observed: Vec<(f64, f64)> = group
        .data
        .gc
        .iter()
        .enumerate()
        .map(|(i, &c)| (i as f64, 100.0 * c as f64 / total as f64))
        .collect();
    let expected: Vec<(f64, f64)> = group
        .data
        .gc_expected
        .iter()
        .enumerate()
        .map(|(i, &e)| (i as f64, 100.0 * e / total as f64))
        .collect();
    let top = observed
        .iter()
        .chain(&expected)
        .map(|p| p.1)
        .fold(0.0, f64::max)
        * 1.08;
    let chart = Chart::new(
        "GC content",
        Axis::linear("GC of read (%)", 0.0, 100.0, Units::Plain),
        Axis::linear("% of reads", 0.0, top, Units::Percent),
    )
    .layer(Layer::Line {
        name: "observed".into(),
        series: 0,
        points: observed,
        style: LineStyle::Solid,
    })
    .layer(Layer::Line {
        name: "normal fit".into(),
        series: 1,
        points: expected,
        style: LineStyle::Dashed,
    });
    figure(
        page,
        "Reads by GC percent",
        chart,
        &[
            ("observed".into(), 0, LineStyle::Solid),
            ("normal fit".into(), 1, LineStyle::Dashed),
        ],
    );
}

fn duplication(page: &mut Page, group: &Group) {
    let Some(d) = &group.data.duplication else {
        return;
    };
    let names: Vec<String> = d.levels.iter().map(|(n, _)| n.clone()).collect();
    let top = d.levels.iter().map(|l| l.1).fold(0.0, f64::max).max(1.0) * 1.08;
    let bars = d
        .levels
        .iter()
        .enumerate()
        .map(|(i, (name, p))| Bar {
            x0: i as f64 - 0.4,
            x1: i as f64 + 0.4,
            y: *p,
            tip: format!("{name} copies: {p:.2}% of sampled reads"),
        })
        .collect();
    let chart = Chart::new(
        "Duplication levels",
        Axis::categories("copies of the read (by both 50-base ends)", names),
        Axis::linear("% of sampled reads", 0.0, top, Units::Percent),
    )
    .layer(Layer::Bars { series: 0, bars });
    figure(
        page,
        &format!(
            "{} sampled reads; {:.2}% would remain after deduplication",
            thousands(d.sampled as u64),
            d.distinct_pct
        ),
        chart,
        &[],
    );
}

fn overrepresented(page: &mut Page, group: &Group) {
    if group.data.overrepresented.is_empty() {
        page.push_str(r#"<p class="note">No read start or end above 0.1% of sampled reads.</p>"#);
        return;
    }
    page.push_str(r#"<div class="table-wrap"><table><tr><th class="l">End</th><th class="l">Sequence (50 bases)</th><th>Reads</th><th>% sampled</th><th class="l">Possible source</th></tr>"#);
    for o in &group.data.overrepresented {
        let end = match o.end {
            ReadEnd::Start => "start",
            ReadEnd::End => "end",
        };
        let _ = write!(
            page,
            r#"<tr><td class="l">{end}</td><td class="l seq">{}</td><td>{}</td><td>{:.2}%</td><td class="l">{}</td></tr>"#,
            o.sequence,
            thousands(o.count),
            o.pct,
            escape(&o.source)
        );
    }
    page.push_str("</table></div>");
}

fn adapters(page: &mut Page, report: &Report, group: &Group) {
    let data = &group.data;
    page.push_str(r#"<div class="table-wrap"><table><tr><th class="l">Adapter or primer</th><th>At start</th><th>At end</th><th>Inside read</th></tr>"#);
    if data.adapters.is_empty() {
        page.push_str(
            r#"<tr><td class="l" colspan="4">No known adapter or primer found.</td></tr>"#,
        );
    }
    for a in &data.adapters {
        let _ = write!(
            page,
            r#"<tr><td class="l">{}</td><td>{:.3}%</td><td>{:.3}%</td><td>{:.3}%</td></tr>"#,
            escape(&a.name),
            a.start_pct,
            a.end_pct,
            a.internal_pct
        );
    }
    let _ = write!(
        page,
        r#"<tr><td class="l">poly(T) of 15+ in the first 100 bases</td><td>{}</td><td></td><td></td></tr><tr><td class="l">poly(A) of 15+ in the last 100 bases</td><td></td><td>{}</td><td></td></tr></table></div>"#,
        opt(data.poly_t_head_pct, 3, "%"),
        opt(data.poly_a_tail_pct, 3, "%")
    );
    if data.adapters.is_empty() {
        return;
    }
    for end in [ReadEnd::Start, ReadEnd::End] {
        let mut chart_top: f64 = 0.0;
        let mut last_bin = 0;
        let mut layers = Vec::new();
        let mut legend = Vec::new();
        for (i, a) in data.adapters.iter().enumerate() {
            let values = match end {
                ReadEnd::Start => &a.cumulative_from_start,
                ReadEnd::End => &a.cumulative_from_end,
            };
            let used = values
                .iter()
                .rposition(|&v| v > 0.0)
                .map_or(0, |p| (p + 2).min(values.len() - 1));
            last_bin = last_bin
                .max(used)
                .max(group.data.content_from_start.last().map_or(0, |c| c.bin));
            chart_top = chart_top.max(values.iter().copied().fold(0.0, f64::max));
            layers.push((a.name.clone(), i, values.clone()));
            legend.push((a.name.clone(), i, LineStyle::Solid));
        }
        let mut chart = Chart::new(
            "Adapter content",
            position_axis(report, last_bin, end),
            Axis::linear(
                "% of reads (cumulative)",
                0.0,
                (chart_top * 1.1).max(0.5),
                Units::Percent,
            ),
        );
        for (name, i, values) in layers {
            let points = values
                .iter()
                .enumerate()
                .take(last_bin + 1)
                .map(|(b, &v)| (bin_x(report, b), v))
                .collect();
            chart = chart.layer(Layer::Line {
                name,
                series: i,
                points,
                style: LineStyle::Solid,
            });
        }
        let caption = match end {
            ReadEnd::Start => "Reads with the adapter at or before each position, from the start",
            ReadEnd::End => "Reads with the adapter at or before each position, from the end",
        };
        figure(page, caption, chart, &legend);
    }
}

fn yield_over_time(page: &mut Page, group: &Group) {
    let time = &group.data.time;
    let Some(last) = time.last() else { return };
    let hours = (last.hours + 0.01).max(0.1);
    let cumulative: Vec<(f64, f64)> = time
        .iter()
        .map(|t| (t.hours, t.cumulative_bases as f64))
        .collect();
    let chart = Chart::new(
        "Cumulative yield",
        Axis::linear("hours", 0.0, hours, Units::Plain),
        Axis::linear("bases", 0.0, last.cumulative_bases as f64 * 1.05, Units::Si),
    )
    .layer(Layer::Line {
        name: "bases".into(),
        series: 0,
        points: cumulative,
        style: LineStyle::Solid,
    });
    figure(page, "Cumulative bases", chart, &[]);

    let width = time.get(1).map_or(1.0 / 60.0, |t| t.hours - time[0].hours);
    let per_hour = |bases: u64| bases as f64 / width;
    let top = time.iter().map(|t| per_hour(t.bases)).fold(0.0, f64::max) * 1.05;
    let bars = time
        .iter()
        .map(|t| Bar {
            x0: t.hours,
            x1: t.hours + width,
            y: per_hour(t.bases),
            tip: format!(
                "{:.2} h: {} reads, {} bases",
                t.hours,
                thousands(t.reads),
                si(t.bases as f64)
            ),
        })
        .collect();
    let chart = Chart::new(
        "Yield rate",
        Axis::linear("hours", 0.0, hours + width, Units::Plain),
        Axis::linear("bases per hour", 0.0, top, Units::Si),
    )
    .layer(Layer::Bars { series: 2, bars });
    figure(page, "Bases per hour", chart, &[]);

    let q: Vec<(f64, f64)> = time
        .iter()
        .filter_map(|t| t.mean_q.map(|q| (t.hours + width / 2.0, q)))
        .collect();
    if !q.is_empty() {
        let top = q.iter().map(|p| p.1).fold(0.0, f64::max).max(15.0) * 1.1;
        let chart = Chart::new(
            "Read Q over time",
            Axis::linear("hours", 0.0, hours + width, Units::Plain),
            Axis::linear("mean read Q", 0.0, top, Units::Plain),
        )
        .layer(quality_zones(group, top))
        .layer(Layer::Line {
            name: "mean read Q".into(),
            series: 0,
            points: q,
            style: LineStyle::Solid,
        });
        figure(page, "Mean read Q", chart, &[]);
    }
    let pass: Vec<(f64, f64)> = time
        .iter()
        .filter_map(|t| t.pass_pct.map(|p| (t.hours + width / 2.0, p)))
        .collect();
    if pass.iter().any(|&(_, p)| p < 100.0) {
        let chart = Chart::new(
            "Pass rate over time",
            Axis::linear("hours", 0.0, hours + width, Units::Plain),
            Axis::linear("% of reads passing", 0.0, 100.0, Units::Plain),
        )
        .layer(Layer::Line {
            name: "pass rate".into(),
            series: 1,
            points: pass,
            style: LineStyle::Solid,
        });
        figure(
            page,
            "Reads passing, as a percent of each slice's reads",
            chart,
            &[],
        );
    }
    let lengths: Vec<(f64, f64)> = time
        .iter()
        .filter_map(|t| t.mean_length.map(|l| (t.hours + width / 2.0, l)))
        .collect();
    let top = lengths.iter().map(|p| p.1).fold(0.0, f64::max) * 1.08;
    let chart = Chart::new(
        "Read length over time",
        Axis::linear("hours", 0.0, hours + width, Units::Plain),
        Axis::linear("mean read length", 0.0, top.max(1.0), Units::Si),
    )
    .layer(Layer::Line {
        name: "mean length".into(),
        series: 3,
        points: lengths,
        style: LineStyle::Solid,
    });
    figure(
        page,
        "Mean read length (shearing or library changes show as drift)",
        chart,
        &[],
    );
}

fn channels(page: &mut Page, group: &Group) {
    let channels = &group.data.channels;
    let Some(max_channel) = channels.iter().map(|c| c.channel).max() else {
        return;
    };
    // Rows of 32 channels for MinION-sized flow cells, wider for PromethION.
    let columns: u32 = match max_channel {
        0..=128 => 16,
        129..=512 => 32,
        _ => 100,
    };
    let rows = max_channel.div_ceil(columns);
    // Channels per hover block and per drawn cell; both divide the row length.
    let (cell_block, hit_block): (u32, u32) = match (page.detail, max_channel) {
        (Detail::Full, 0..=512) => (1, 1),
        (Detail::Full, _) => (1, 10),
        (Detail::Compact, 0..=128) => (1, 1),
        (Detail::Compact, 129..=512) => (2, 8),
        (Detail::Compact, _) => (10, 20),
    };
    // Reads and bases of each block of channels, keyed by row and block within the row.
    let blocks = |size: u32| {
        let mut blocks: BTreeMap<(u32, u32), (u64, u64)> = BTreeMap::new();
        for c in channels {
            let index = c.channel.saturating_sub(1);
            let block = blocks
                .entry((index / columns, index % columns / size))
                .or_default();
            block.0 += c.reads;
            block.1 += c.bases;
        }
        blocks
    };
    let edges = |row: u32, block: u32, size: u32| {
        let x0 = f64::from(block * size);
        let y0 = f64::from(rows - row) - 1.0;
        (x0, x0 + f64::from(size), y0, y0 + 1.0)
    };
    let drawn = blocks(cell_block);
    let max_reads = drawn.values().map(|b| b.0).max().unwrap_or(1) as f64;
    let cells = drawn
        .iter()
        .map(|(&(row, block), &(reads, _))| {
            let (x0, x1, y0, y1) = edges(row, block, cell_block);
            Cell {
                x0,
                x1: x1 - 0.08,
                y0,
                y1: y1 - 0.08,
                color: chart::ramp(reads as f64 / max_reads, page.detail.shades()),
            }
        })
        .collect();
    let hits = blocks(hit_block)
        .iter()
        .map(|(&(row, block), &(reads, bases))| {
            let (x0, x1, y0, y1) = edges(row, block, hit_block);
            let first = row * columns + block * hit_block + 1;
            let last = (first + hit_block - 1).min(max_channel);
            let name = if last > first {
                format!("channels {first}–{last}")
            } else {
                format!("channel {first}")
            };
            Hit {
                x0,
                x1,
                y0,
                y1,
                tip: format!(
                    "{name}: {} reads, {} bases",
                    thousands(reads),
                    si(bases as f64)
                ),
            }
        })
        .collect();
    let height = (f64::from(rows) * 14.0).clamp(160.0, 420.0) + 60.0;
    let chart = Chart::new(
        "Channels",
        Axis::linear(
            &format!("channel number, {columns} per row"),
            0.0,
            f64::from(columns),
            Units::Plain,
        ),
        Axis::linear("row", 0.0, f64::from(rows), Units::Plain),
    )
    .layer(Layer::Cells(cells))
    .layer(Layer::Hits(hits))
    .tall(height);
    let blocked = if cell_block > 1 {
        format!(", drawn in blocks of {cell_block}")
    } else {
        String::new()
    };
    figure(
        page,
        &format!(
            "{} channels with reads{blocked} (color: reads, dark to bright)",
            channels.len()
        ),
        chart,
        &[],
    );

    let mut per_channel: Vec<u64> = channels.iter().map(|c| c.reads).collect();
    per_channel.sort_unstable();
    let top_reads = per_channel.last().copied().unwrap_or(1) as f64;
    let width = (top_reads / 40.0).max(1.0).ceil();
    let mut hist = vec![0u64; (top_reads / width) as usize + 1];
    for &r in &per_channel {
        hist[(r as f64 / width) as usize] += 1;
    }
    let ymax = hist.iter().copied().max().unwrap_or(1) as f64 * 1.08;
    let bars = hist
        .iter()
        .enumerate()
        .map(|(i, &n)| Bar {
            x0: i as f64 * width,
            x1: (i + 1) as f64 * width,
            y: n as f64,
            tip: format!(
                "{}–{} reads: {n} channels",
                si(i as f64 * width),
                si((i + 1) as f64 * width)
            ),
        })
        .collect();
    let chart = Chart::new(
        "Reads per channel",
        Axis::linear(
            "reads in the channel",
            0.0,
            hist.len() as f64 * width,
            Units::Si,
        ),
        Axis::linear("channels", 0.0, ymax, Units::Si),
    )
    .layer(Layer::Bars { series: 5, bars });
    figure(page, "Spread of reads per channel", chart, &[]);
}

fn speed(page: &mut Page, group: &Group) {
    let speed = &group.data.speed;
    let total: u64 = speed.iter().sum();
    if total == 0 {
        return;
    }
    let top_bin = speed.iter().rposition(|&c| c > 0).unwrap_or(0) + 1;
    let values: Vec<f64> = speed
        .iter()
        .map(|&c| 100.0 * c as f64 / total as f64)
        .collect();
    let ymax = values.iter().copied().fold(0.0, f64::max) * 1.08;
    let bars = values
        .iter()
        .enumerate()
        .take(top_bin)
        .map(|(i, &p)| {
            let lo = i as f64 * SPEED_BIN_WIDTH;
            Bar {
                x0: lo,
                x1: lo + SPEED_BIN_WIDTH,
                y: p,
                tip: format!("{lo:.0}–{:.0} b/s: {p:.2}% of reads", lo + SPEED_BIN_WIDTH),
            }
        })
        .collect();
    let chart = Chart::new(
        "Translocation speed",
        Axis::linear(
            "bases per second",
            0.0,
            top_bin as f64 * SPEED_BIN_WIDTH,
            Units::Plain,
        ),
        Axis::linear("% of reads", 0.0, ymax, Units::Percent),
    )
    .layer(Layer::Bars { series: 3, bars });
    figure(page, "Read length over read duration", chart, &[]);
}

fn end_reasons(page: &mut Page, group: &Group) {
    let total: u64 = group.data.end_reasons.iter().map(|r| r.1).sum();
    page.push_str(r#"<div class="table-wrap"><table><tr><th class="l">End reason</th><th>Reads</th><th>%</th></tr>"#);
    for (reason, count) in &group.data.end_reasons {
        let _ = write!(
            page,
            r#"<tr><td class="l">{}</td><td>{}</td><td>{:.2}%</td></tr>"#,
            escape(reason),
            thousands(*count),
            100.0 * *count as f64 / total.max(1) as f64
        );
    }
    page.push_str("</table></div>");
}

fn histogram_chart(title: &str, x_label: &str, values: &[u64], series: usize) -> Option<Chart> {
    let total: u64 = values.iter().sum();
    if total == 0 {
        return None;
    }
    let pcts: Vec<f64> = values
        .iter()
        .map(|&c| 100.0 * c as f64 / total as f64)
        .collect();
    let ymax = pcts.iter().copied().fold(0.0, f64::max) * 1.08;
    let first = values.iter().position(|&c| c > 0).unwrap_or(0);
    let bars = pcts
        .iter()
        .enumerate()
        .skip(first)
        .map(|(i, &p)| Bar {
            x0: i as f64 - 0.45,
            x1: i as f64 + 0.45,
            y: p,
            tip: format!("{i}: {} ({p:.2}%)", thousands(values[i])),
        })
        .collect();
    Some(
        Chart::new(
            title,
            Axis::linear(
                x_label,
                first as f64 - 0.5,
                values.len() as f64 - 0.5,
                Units::Plain,
            ),
            Axis::linear("% of reads", 0.0, ymax, Units::Percent),
        )
        .layer(Layer::Bars { series, bars }),
    )
}

fn pacbio(page: &mut Page, group: &Group) {
    let Some(p) = &group.data.pacbio else { return };
    let _ = write!(
        page,
        r#"<p class="note">{:.1}% of reads carry kinetics tags.</p>"#,
        p.kinetics_pct
    );
    for (title, label, values, series) in [
        (
            "CCS passes",
            "passes (np; 100 means 100 or more)",
            &p.passes,
            0,
        ),
        (
            "Effective coverage",
            "effective coverage (ec)",
            &p.effective_coverage,
            2,
        ),
        (
            "Barcode score",
            "lima barcode score (bq)",
            &p.barcode_score,
            3,
        ),
    ] {
        if let Some(chart) = histogram_chart(title, label, values, series) {
            figure(page, title, chart, &[]);
        }
    }
    if !p.passes_q.is_empty() {
        let max_passes = p.passes_q.iter().map(|c| c.0).max().unwrap_or(10) + 1;
        let max_q = p.passes_q.iter().map(|c| c.1).max().unwrap_or(30) + 1;
        let hit_merge = match page.detail {
            Detail::Full => (1, 1),
            Detail::Compact => (2, 2),
        };
        let [cells, hits] = count_heatmap(
            page.detail.shades(),
            &p.passes_q,
            (1, 1),
            hit_merge,
            &|passes, span| (passes as f64 - 0.5, (passes + span) as f64 - 0.5),
            &|passes, q, (span, qs)| {
                if span > 1 {
                    let last = passes + span - 1;
                    format!("{passes}–{last} passes, {}", q_range(q, qs))
                } else {
                    format!("{passes} passes, {}", q_range(q, qs))
                }
            },
        );
        let chart = Chart::new(
            "Passes vs read Q",
            Axis::linear(
                "passes (60 means 60 or more)",
                -0.5,
                max_passes as f64 - 0.5,
                Units::Plain,
            ),
            Axis::linear("read Q", 0.0, max_q as f64, Units::Plain),
        )
        .layer(cells)
        .layer(hits);
        figure(
            page,
            "Reads by passes and read Q (color: log read count)",
            chart,
            &[],
        );
    }
    if p.zmws > 0 {
        let total: u64 = p.reads_per_zmw.iter().sum();
        page.push_str(
            r#"<div class="table-wrap"><table><tr><th class="l">Reads per ZMW</th><th>ZMWs</th><th>%</th></tr>"#,
        );
        for (i, &n) in p.reads_per_zmw.iter().enumerate() {
            let label = if i + 1 == p.reads_per_zmw.len() {
                format!("{} or more", i + 1)
            } else {
                (i + 1).to_string()
            };
            let _ = write!(
                page,
                r#"<tr><td class="l">{label}</td><td>{}</td><td>{:.2}%</td></tr>"#,
                thousands(n),
                100.0 * n as f64 / total.max(1) as f64
            );
        }
        page.push_str(
            r#"</table></div><p class="explain">Two reads from one ZMW are by-strand HiFi reads (one per strand) or the same molecule read twice.</p>"#,
        );
    }
}

fn modifications(page: &mut Page, group: &Group) {
    let mods = &group.data.mods;
    page.push_str(r#"<div class="table-wrap"><table><tr><th class="l">Modification</th><th>Calls</th><th>Calls per read</th><th>Called modified</th><th>CpG calls</th><th>CpG called modified</th></tr>"#);
    for m in mods {
        let _ = write!(
            page,
            r#"<tr><td class="l">{}</td><td>{}</td><td>{:.0}</td><td>{:.2}%</td><td>{}</td><td>{}</td></tr>"#,
            escape(&m.label),
            thousands(m.calls),
            m.calls_per_read,
            m.modified_pct,
            if m.cpg_calls > 0 {
                thousands(m.cpg_calls)
            } else {
                "–".to_owned()
            },
            opt(m.cpg_modified_pct, 2, "%")
        );
    }
    page.push_str("</table></div>");
    let mut chart_top: f64 = 0.0;
    let mut lines = Vec::new();
    let mut legend = Vec::new();
    for (i, m) in mods.iter().enumerate().take(PALETTE) {
        let total: u64 = m.probs.iter().sum();
        let points: Vec<(f64, f64)> = m
            .probs
            .iter()
            .enumerate()
            .map(|(b, &c)| {
                (
                    (b as f64 * 10.0 + 5.0).min(255.0) / 256.0,
                    100.0 * c as f64 / total.max(1) as f64,
                )
            })
            .collect();
        chart_top = chart_top.max(points.iter().map(|p| p.1).fold(0.0, f64::max));
        lines.push(Layer::Line {
            name: m.label.clone(),
            series: i,
            points,
            style: LineStyle::Solid,
        });
        legend.push((m.label.clone(), i, LineStyle::Solid));
    }
    if lines.is_empty() {
        return;
    }
    let mut chart = Chart::new(
        "Modification probabilities",
        Axis::linear(
            "probability of the modification (ML)",
            0.0,
            1.0,
            Units::Plain,
        ),
        Axis::linear("% of calls", 0.0, chart_top * 1.08, Units::Percent),
    );
    for line in lines {
        chart = chart.layer(line);
    }
    figure(page, "Calls by modification probability", chart, &legend);
}

fn alignment(page: &mut Page, group: &Group) {
    let Some(m) = &group.data.mapping else { return };
    let _ = write!(
        page,
        r#"<div class="table-wrap"><table><tr><th class="l">Primary mapped</th><th>Unmapped</th><th>Secondary</th><th>Supplementary</th><th>Split reads</th><th>Median identity</th><th>Median gap-compressed</th></tr><tr><td class="l">{} ({:.2}%)</td><td>{}</td><td>{}</td><td>{}</td><td>{:.2}%</td><td>{}</td><td>{}</td></tr></table></div>"#,
        thousands(m.mapped),
        m.mapped_pct,
        thousands(m.unmapped),
        thousands(m.secondary),
        thousands(m.supplementary),
        m.split_pct,
        opt(m.median_identity_pct, 2, "%"),
        opt(m.median_gap_compressed_pct, 2, "%"),
    );
    let total: u64 = m.identity.iter().sum();
    if total > 0 {
        let to_points = |values: &[u64]| -> Vec<(f64, f64)> {
            values
                .iter()
                .enumerate()
                .map(|(i, &c)| {
                    (
                        m.identity_floor + i as f64 / 10.0,
                        100.0 * c as f64 / total as f64,
                    )
                })
                .collect()
        };
        let blast = to_points(&m.identity);
        let gap = to_points(&m.gap_compressed_identity);
        let top = blast.iter().chain(&gap).map(|p| p.1).fold(0.0, f64::max) * 1.08;
        let chart = Chart::new(
            "Identity",
            Axis::linear("identity (%)", m.identity_floor, 100.0, Units::Plain),
            Axis::linear("% of mapped reads", 0.0, top, Units::Percent),
        )
        .layer(Layer::Line {
            name: "BLAST identity".into(),
            series: 0,
            points: blast,
            style: LineStyle::Thin,
        })
        .layer(Layer::Line {
            name: "gap-compressed".into(),
            series: 1,
            points: gap,
            style: LineStyle::Thin,
        });
        figure(
            page,
            "Mapped reads by identity",
            chart,
            &[
                ("BLAST identity".into(), 0, LineStyle::Solid),
                ("gap-compressed".into(), 1, LineStyle::Solid),
            ],
        );
    }
    if let Some(chart) =
        histogram_chart("MAPQ", "mapping quality (60 means 60 or more)", &m.mapq, 2)
    {
        figure(page, "Primary alignments by MAPQ", chart, &[]);
    }
    if let Some(chart) = histogram_chart(
        "Aligned fraction",
        "percent of the read in the primary alignment",
        &m.aligned_fraction,
        3,
    ) {
        figure(page, "Fraction of each read aligned", chart, &[]);
    }
}

/// The overview: run tiles, status matrix, per-group table and comparisons.
fn overview(page: &mut Page, report: &Report) {
    let all = &report.groups[0];
    let s = &all.summary;
    page.push_str(r#"<section class="group" id="overview"><h2>Overview</h2>"#);
    let barcodes = report
        .groups
        .iter()
        .filter(|g| g.kind == GroupKind::Barcode)
        .count();
    page.push_str(r#"<div class="tiles">"#);
    tile(page, &thousands(s.reads), "reads");
    tile(page, &si(s.bases as f64), "bases");
    tile(page, &si(s.n50 as f64), "read N50");
    tile(page, &opt(s.median_read_q, 1, ""), "median read Q");
    if report.barcoded {
        tile(page, &barcodes.to_string(), "barcodes");
        if let Some(u) = report
            .groups
            .iter()
            .find(|g| g.kind == GroupKind::Unclassified)
        {
            tile(
                page,
                &format!("{:.1}%", u.summary.share_pct),
                "unclassified",
            );
        }
    }
    page.push_str("</div>");

    let rows = rows(report);
    status_matrix(page, &rows);
    summary_table(page, &rows);
    if report.barcoded && report.groups.len() > 2 {
        let groups: Vec<&Group> = report.groups.iter().skip(1).collect();
        let labels: Vec<String> = groups.iter().map(|g| display_name(g)).collect();
        comparisons(page, report, &groups, &labels, "Barcodes compared");
    }
    extras::run_info(page, report);
    page.push_str("</section>");
}

/// One line of the overview tables: a group, its label and where it links.
struct Row<'a> {
    label: String,
    link: Option<String>,
    group: &'a Group,
}

fn rows(report: &Report) -> Vec<Row<'_>> {
    report
        .groups
        .iter()
        .enumerate()
        .map(|(index, group)| Row {
            label: display_name(group),
            link: Some(group_id(index, group)),
            group,
        })
        .collect()
}

fn row_label(row: &Row<'_>) -> String {
    match &row.link {
        Some(link) => format!(r##"<a href="#{link}">{}</a>"##, escape(&row.label)),
        None => escape(&row.label),
    }
}

fn status_matrix(page: &mut Page, rows: &[Row<'_>]) {
    let shown: Vec<ModuleId> = ModuleId::ALL
        .iter()
        .copied()
        .filter(|&id| {
            rows.iter().any(|r| {
                r.group
                    .module(id)
                    .is_some_and(|m| m.status != Status::Absent)
            })
        })
        .collect();
    page.push_str(r#"<div class="panel"><h3>Checks</h3>"#);
    let many = rows.len() > MANY_GROUPS;
    if many {
        check_tally(page, rows, &shown);
    }
    let _ = write!(
        page,
        r#"<div class="table-wrap"><table class="matrix{}"><tr><th class="l">Group</th>"#,
        if many { " dense" } else { "" }
    );
    for id in &shown {
        let _ = write!(page, r#"<th class="rot"><div>{}</div></th>"#, id.title());
    }
    page.push_str("</tr>");
    for row in rows {
        let _ = write!(page, r#"<tr><td class="l">{}</td>"#, row_label(row));
        for &id in &shown {
            match row.group.module(id) {
                Some(m) if m.status != Status::Absent => {
                    let href = row.link.as_ref().map_or_else(String::new, |link| {
                        format!(r##" href="#{link}-{}""##, id.key())
                    });
                    let _ = write!(
                        page,
                        r#"<td><a class="dot {}"{href} title="{}: {}"></a></td>"#,
                        status_class(m.status),
                        m.status.label(),
                        escape(&m.note)
                    );
                }
                _ => page.push_str("<td></td>"),
            }
        }
        page.push_str("</tr>");
    }
    page.push_str(r#"</table></div><div class="legend"><span><span class="dot pass"></span>pass</span><span><span class="dot warn"></span>warn</span><span><span class="dot fail"></span>fail</span><span><span class="dot info"></span>shown without a verdict</span></div></div>"#);
}

/// Groups a check's tally names before it says how many more there are.
const TALLY_NAMES: usize = 8;

/// Per check, how many groups pass, warn and fail, and which ones warn or fail: the
/// short answer when the matrix has too many rows to take in.
fn check_tally(page: &mut Page, rows: &[Row<'_>], shown: &[ModuleId]) {
    page.push_str(r#"<div class="table-wrap"><table class="tally"><tr><th class="l">Check</th><th>Pass</th><th>Warn</th><th>Fail</th><th class="l">Share of groups</th><th class="l">Groups that warn or fail</th></tr>"#);
    for &id in shown {
        let graded = |status: Status| {
            rows.iter()
                .filter(move |r| r.group.module(id).is_some_and(|m| m.status == status))
        };
        let counts = [Status::Pass, Status::Warn, Status::Fail].map(|s| graded(s).count());
        let total: usize = counts.iter().sum();
        if total == 0 {
            continue;
        }
        let _ = write!(
            page,
            r#"<tr><td class="l">{}</td><td>{}</td><td>{}</td><td>{}</td><td class="l"><div class="stack mini">"#,
            id.title(),
            counts[0],
            counts[1],
            counts[2]
        );
        for (count, class) in counts.iter().zip(["pass", "warn", "fail"]) {
            if *count > 0 {
                let _ = write!(
                    page,
                    r#"<div class="{class}" style="width:{:.2}%"></div>"#,
                    100.0 * *count as f64 / total as f64
                );
            }
        }
        page.push_str(r#"</div></td><td class="l names">"#);
        let flagged: Vec<(&Row<'_>, Status)> = graded(Status::Fail)
            .map(|r| (r, Status::Fail))
            .chain(graded(Status::Warn).map(|r| (r, Status::Warn)))
            .collect();
        for (row, status) in flagged.iter().take(TALLY_NAMES) {
            let name = escape(&row.label);
            let class = status_class(*status);
            match &row.link {
                Some(link) => {
                    let _ = write!(
                        page,
                        r##"<a class="flag {class}" href="#{link}-{}">{name}</a> "##,
                        id.key()
                    );
                }
                None => {
                    let _ = write!(page, r#"<span class="flag {class}">{name}</span> "#);
                }
            }
        }
        if flagged.len() > TALLY_NAMES {
            let _ = write!(page, "+{} more", flagged.len() - TALLY_NAMES);
        }
        page.push_str("</td></tr>");
    }
    page.push_str("</table></div>");
}

fn summary_table(page: &mut Page, rows: &[Row<'_>]) {
    page.push_str(r#"<div class="panel"><h3>Summary</h3><div class="table-wrap"><table><tr><th class="l">Group</th><th>Reads</th><th>% reads</th><th>Bases</th><th>Mean length</th><th>N50</th><th>Median Q</th><th>≥Q20 reads</th><th>GC</th>"#);
    let any = |f: &dyn Fn(&Group) -> bool| rows.iter().any(|r| f(r.group));
    let mapped = any(&|g| g.summary.mapped_pct.is_some());
    let methylated = any(&|g| g.summary.methylated_pct.is_some());
    let coverage = any(&|g| g.summary.coverage.is_some());
    let contamination = any(&|g| !g.data.contaminants.is_empty());
    if mapped {
        page.push_str("<th>Mapped</th><th>Identity</th>");
    }
    if methylated {
        page.push_str("<th>CpG 5mC</th>");
    }
    if coverage {
        page.push_str("<th>Coverage</th>");
    }
    if contamination {
        page.push_str("<th>Contaminants</th>");
    }
    page.push_str("</tr>");
    for row in rows {
        let group = row.group;
        let s = &group.summary;
        let _ = write!(
            page,
            r#"<tr><td class="l">{}</td><td>{}</td><td>{:.2}%</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td>"#,
            row_label(row),
            thousands(s.reads),
            s.share_pct,
            si(s.bases as f64),
            si(s.mean_length),
            si(s.n50 as f64),
            opt(s.median_read_q, 1, ""),
            opt(s.reads_q20_pct, 1, "%"),
            opt(s.gc_pct, 1, "%"),
        );
        if mapped {
            let _ = write!(
                page,
                "<td>{}</td><td>{}</td>",
                opt(s.mapped_pct, 2, "%"),
                opt(s.median_identity_pct, 2, "%")
            );
        }
        if methylated {
            let _ = write!(page, "<td>{}</td>", opt(s.methylated_pct, 1, "%"));
        }
        if coverage {
            let _ = write!(page, "<td>{}</td>", opt(s.coverage, 1, "x"));
        }
        if contamination {
            let contaminant =
                super::checks::foreign_contaminant_pct(&group.data.contaminants).unwrap_or(0.0);
            let _ = write!(page, "<td>{contaminant:.2}%</td>");
        }
        page.push_str("</tr>");
    }
    page.push_str("</table></div></div>");
}

/// Line style and color of the `i`th barcode in overlays.
fn barcode_style(i: usize) -> (usize, LineStyle) {
    let style = if (i / PALETTE) % 2 == 1 {
        LineStyle::Dashed
    } else {
        LineStyle::Solid
    };
    (i % PALETTE, style)
}

/// How a barcode's curve reads at an x where it has no point.
#[derive(Debug, Clone, Copy)]
enum Gap {
    /// Histograms: no bin means no reads.
    Zero,
    /// Per-position values: the barcode has no data there.
    Missing,
    /// Cumulative curves: the last value stands.
    Hold,
}

/// One overlay chart of the group comparison: a curve per group on shared axes.
struct Overlay<'a> {
    title: &'a str,
    caption: &'a str,
    x: Axis,
    y_label: &'a str,
    units: Units,
    gap: Gap,
    curves: Vec<Vec<(f64, f64)>>,
}

/// Barcodes with fewer reads are too noisy to single out as outliers.
const OUTLIER_MIN_READS: u64 = 200;

/// Most outlier barcodes drawn over a band.
const MAX_OUTLIERS: usize = 5;

/// A curve is an outlier when it leaves the 10-90% band by this share of the chart height.
const OUTLIER_EXCURSION: f64 = 0.06;

/// The value at `q` (0..=1) of sorted values, interpolating between neighbors.
fn quantile(sorted: &[f64], q: f64) -> f64 {
    let at = q * (sorted.len() - 1) as f64;
    let (low, high) = (at.floor() as usize, at.ceil() as usize);
    sorted[low] + (sorted[high] - sorted[low]) * (at - low as f64)
}

/// A curve's values at each of `xs` (sorted), filling gaps as `gap` says.
fn resample(curve: &[(f64, f64)], xs: &[f64], gap: Gap) -> Vec<Option<f64>> {
    let mut next = 0;
    xs.iter()
        .map(|&x| {
            while next < curve.len() && curve[next].0 < x {
                next += 1;
            }
            match curve.get(next) {
                // Groups share their binning, so the same bin has the same x in each.
                Some(&(at, y)) if at == x => Some(y),
                _ => match gap {
                    Gap::Zero => Some(0.0),
                    Gap::Missing => None,
                    Gap::Hold => Some(next.checked_sub(1).map_or(0.0, |before| curve[before].1)),
                },
            }
        })
        .collect()
}

impl Overlay<'_> {
    /// One line per group, for reports with few groups.
    fn lines(self, page: &mut Page, labels: &[String], legend: &[(String, usize, LineStyle)]) {
        let top = self
            .curves
            .iter()
            .flatten()
            .map(|p| p.1)
            .fold(0.0, f64::max)
            .max(1e-9)
            * 1.08;
        let mut chart = Chart::new(
            self.title,
            self.x,
            Axis::linear(self.y_label, 0.0, top, self.units),
        );
        for (i, points) in self.curves.into_iter().enumerate() {
            let (series, style) = barcode_style(i);
            chart = chart.layer(Layer::Line {
                name: labels[i].clone(),
                series,
                points,
                style,
            });
        }
        figure(page, self.caption, chart, legend);
    }

    /// The median and spread of the barcodes as bands, with only the outlier barcodes
    /// and the unassigned reads drawn as lines: readable however many barcodes there are.
    fn band(self, page: &mut Page, groups: &[&Group], labels: &[String]) {
        let mut xs: Vec<f64> = self.curves.iter().flatten().map(|p| p.0).collect();
        xs.sort_by(f64::total_cmp);
        xs.dedup();
        let sampled: Vec<Vec<Option<f64>>> = self
            .curves
            .iter()
            .map(|curve| resample(curve, &xs, self.gap))
            .collect();
        let is_member =
            |group: &Group| !matches!(group.kind, GroupKind::Unclassified | GroupKind::NoBarcode);
        let members: Vec<usize> = (0..groups.len())
            .filter(|&i| is_member(groups[i]) && !self.curves[i].is_empty())
            .collect();

        // Quantiles of the member curves at each x where at least a quarter have data.
        let spread: Vec<Option<[f64; 5]>> = (0..xs.len())
            .map(|at| {
                let mut values: Vec<f64> = members.iter().filter_map(|&i| sampled[i][at]).collect();
                values.sort_by(f64::total_cmp);
                (values.len() >= (members.len() / 4).max(2))
                    .then(|| [0.1, 0.25, 0.5, 0.75, 0.9].map(|q| quantile(&values, q)))
            })
            .collect();
        let band_top = spread.iter().flatten().map(|s| s[4]).fold(0.0, f64::max);

        // How far each member strays outside the 10-90% band, at its worst.
        let excursion = |i: usize| {
            sampled[i]
                .iter()
                .zip(&spread)
                .filter_map(|(value, spread)| {
                    let (value, spread) = ((*value)?, (*spread)?);
                    Some((spread[0] - value).max(value - spread[4]))
                })
                .fold(0.0, f64::max)
        };
        let mut outliers: Vec<(usize, f64)> = members
            .iter()
            .filter(|&&i| groups[i].summary.reads >= OUTLIER_MIN_READS)
            .map(|&i| (i, excursion(i)))
            .filter(|&(_, excursion)| excursion > OUTLIER_EXCURSION * band_top.max(1e-9))
            .collect();
        outliers.sort_by(|a, b| b.1.total_cmp(&a.1));
        outliers.truncate(MAX_OUTLIERS);
        outliers.sort_by_key(|&(i, _)| i);
        let others: Vec<usize> = (0..groups.len())
            .filter(|&i| !is_member(groups[i]) && !self.curves[i].is_empty())
            .collect();

        let drawn_top = outliers
            .iter()
            .map(|&(i, _)| i)
            .chain(others.iter().copied())
            .flat_map(|i| self.curves[i].iter().map(|p| p.1))
            .fold(band_top, f64::max);
        let band = |low: usize, high: usize, opacity: f64| Layer::Band {
            series: 0,
            points: xs
                .iter()
                .zip(&spread)
                .filter_map(|(&x, spread)| spread.map(|s| (x, s[low], s[high])))
                .collect(),
            opacity,
        };
        let mut chart = Chart::new(
            self.title,
            self.x,
            Axis::linear(self.y_label, 0.0, drawn_top.max(1e-9) * 1.08, self.units),
        )
        .layer(band(0, 4, 0.16))
        .layer(band(1, 3, 0.3))
        .layer(Layer::Line {
            name: "median".into(),
            series: 0,
            points: xs
                .iter()
                .zip(&spread)
                .filter_map(|(&x, spread)| spread.map(|s| (x, s[2])))
                .collect(),
            style: LineStyle::Solid,
        });
        let mut legend: Vec<(String, usize, LineStyle)> = Vec::new();
        for (rank, &(i, _)) in outliers.iter().enumerate() {
            legend.push((labels[i].clone(), rank + 1, LineStyle::Thin));
        }
        for &i in &others {
            legend.push((labels[i].clone(), PALETTE - 1, LineStyle::Dashed));
        }
        let lines = outliers.iter().map(|&(i, _)| i).chain(others);
        for (i, (name, series, style)) in lines.zip(&legend) {
            chart = chart.layer(Layer::Line {
                name: name.clone(),
                series: *series,
                points: self.curves[i].clone(),
                style: *style,
            });
        }

        let mut keys = format!(
            r#"<div class="legend"><span class="s0"><i></i>median of {} groups</span><span class="s0"><i class="box mid"></i>25–75%</span><span class="s0"><i class="box outer"></i>10–90%</span>"#,
            members.len()
        );
        for (name, series, style) in &legend {
            keys.push_str(&legend_item(name, *series, *style));
        }
        if outliers.is_empty() {
            keys.push_str("<span>no outlier groups</span>");
        }
        keys.push_str("</div>");
        page.figures.push(Figure {
            caption: self.caption.to_owned(),
            chart,
            legend: keys,
            span: Span::Column,
        });
    }
}

fn comparisons(
    page: &mut Page,
    report: &Report,
    groups: &[&Group],
    labels: &[String],
    heading: &str,
) {
    let many = groups.len() > MANY_GROUPS;
    let _ = write!(page, r#"<div class="panel"><h3>{}</h3>"#, escape(heading));
    if many {
        let _ = write!(
            page,
            r#"<p class="explain">With {} groups, the curves show the median group and the spread of all of them. Only groups that leave the 10–90% band by a wide margin are drawn by name (at most {MAX_OUTLIERS}, and none with fewer than {OUTLIER_MIN_READS} reads), with the unassigned reads dashed. Every group has its own section with the full charts.</p>"#,
            groups.len()
        );
    }

    let mut bar_chart =
        |caption: &str, label: &str, units: Units, value: &dyn Fn(&Group) -> Option<f64>| {
            let values: Vec<Option<f64>> = groups.iter().map(|g| value(g)).collect();
            let top = values
                .iter()
                .flatten()
                .copied()
                .fold(0.0, f64::max)
                .max(1e-9)
                * 1.08;
            let text = |v: f64| match units {
                Units::Si => si(v),
                Units::Percent => format!("{v:.2}%"),
                Units::Plain => format!("{v:.2}"),
            };
            let bars = values
                .iter()
                .enumerate()
                .filter_map(|(i, v)| {
                    v.map(|v| Bar {
                        x0: i as f64 - 0.4,
                        x1: i as f64 + 0.4,
                        y: v,
                        tip: format!("{}: {}", labels[i], text(v)),
                    })
                })
                .collect();
            let mut chart = Chart::new(
                caption,
                Axis::categories("", labels.to_vec()),
                Axis::linear(label, 0.0, top, units),
            )
            .layer(Layer::Bars { series: 0, bars });
            if many {
                let mut present: Vec<f64> = values.iter().flatten().copied().collect();
                present.sort_by(f64::total_cmp);
                if !present.is_empty() {
                    let median = quantile(&present, 0.5);
                    chart = chart.layer(Layer::Level {
                        name: format!("median {}", text(median)),
                        series: 1,
                        y: median,
                    });
                }
                page.figures.push(Figure {
                    caption: caption.to_owned(),
                    chart,
                    legend: String::new(),
                    span: Span::Row,
                });
            } else {
                figure(page, caption, chart, &[]);
            }
        };
    bar_chart("Reads per group", "reads", Units::Si, &|g| {
        Some(g.summary.reads as f64)
    });
    bar_chart("Bases per group", "bases", Units::Si, &|g| {
        Some(g.summary.bases as f64)
    });
    bar_chart("Read N50 per group", "read N50", Units::Si, &|g| {
        Some(g.summary.n50 as f64)
    });
    bar_chart(
        "Median read Q per group",
        "median read Q",
        Units::Plain,
        &|g| g.summary.median_read_q,
    );

    let legend: Vec<(String, usize, LineStyle)> = labels
        .iter()
        .enumerate()
        .map(|(i, n)| {
            let (series, style) = barcode_style(i);
            (n.clone(), series, style)
        })
        .collect();
    let mut draw = |overlay: Overlay<'_>| {
        if many {
            overlay.band(page, groups, labels);
        } else {
            overlay.lines(page, labels, &legend);
        }
    };
    let curves = |curve: &dyn Fn(&Group) -> Vec<(f64, f64)>| -> Vec<Vec<(f64, f64)>> {
        groups.iter().map(|g| curve(g)).collect()
    };

    let min_len = groups
        .iter()
        .filter_map(|g| g.data.length_bins.first())
        .map(|b| b.start)
        .fold(f64::INFINITY, f64::min);
    let max_len = groups
        .iter()
        .filter_map(|g| g.data.length_bins.last())
        .map(|b| b.end)
        .fold(0.0, f64::max);
    if min_len.is_finite() {
        draw(Overlay {
            title: "Bases by length",
            caption: "Bases by read length",
            x: Axis::log("read length (bases)", min_len, max_len, Units::Si),
            y_label: "% of bases",
            units: Units::Percent,
            gap: Gap::Zero,
            curves: curves(&|g| {
                let total = g.summary.bases.max(1) as f64;
                g.data
                    .length_bins
                    .iter()
                    .map(|b| ((b.start * b.end).sqrt(), 100.0 * b.bases as f64 / total))
                    .collect()
            }),
        });
    }
    let max_q = groups
        .iter()
        .filter_map(|g| g.data.read_q.last())
        .map(|b| b.q)
        .fold(15.0, f64::max);
    if groups.iter().any(|g| !g.data.read_q.is_empty()) {
        draw(Overlay {
            title: "Reads by read Q",
            caption: "Reads by read Q",
            x: Axis::linear("read Q", 0.0, max_q.ceil(), Units::Plain),
            y_label: "% of reads",
            units: Units::Percent,
            gap: Gap::Zero,
            curves: curves(&|g| {
                let total = g.data.read_q.iter().map(|b| b.reads).sum::<u64>().max(1) as f64;
                // Whole-Q bins keep overlays readable.
                let mut by_q: Vec<(f64, f64)> = Vec::new();
                for b in &g.data.read_q {
                    let q = b.q.floor() + 0.5;
                    match by_q.last_mut() {
                        Some(last) if (last.0 - q).abs() < 1e-9 => {
                            last.1 += 100.0 * b.reads as f64 / total;
                        }
                        _ => by_q.push((q, 100.0 * b.reads as f64 / total)),
                    }
                }
                by_q
            }),
        });
    }
    let last_bin = groups
        .iter()
        .filter_map(|g| g.data.quality_from_start.last())
        .map(|p| p.bin)
        .max();
    if let Some(last_bin) = last_bin {
        draw(Overlay {
            title: "Median base quality",
            caption: "Median base quality from the read start",
            x: position_axis(report, last_bin, ReadEnd::Start),
            y_label: "median base Q",
            units: Units::Plain,
            gap: Gap::Missing,
            curves: curves(&|g| {
                g.data
                    .quality_from_start
                    .iter()
                    .map(|p| (bin_x(report, p.bin), p.median))
                    .collect()
            }),
        });
    }
    if groups.iter().any(|g| g.summary.gc_pct.is_some()) {
        draw(Overlay {
            title: "GC content",
            caption: "Reads by GC percent",
            x: Axis::linear("GC of read (%)", 0.0, 100.0, Units::Plain),
            y_label: "% of reads",
            units: Units::Percent,
            gap: Gap::Zero,
            curves: curves(&|g| {
                let total = g.data.gc.iter().sum::<u64>().max(1) as f64;
                g.data
                    .gc
                    .iter()
                    .enumerate()
                    .map(|(i, &c)| (i as f64, 100.0 * c as f64 / total))
                    .collect()
            }),
        });
    }
    let max_hours = groups
        .iter()
        .filter_map(|g| g.data.time.last())
        .map(|t| t.hours)
        .fold(0.0, f64::max);
    if max_hours > 0.0 {
        draw(Overlay {
            title: "Cumulative yield",
            caption: "Cumulative bases over the run",
            x: Axis::linear("hours", 0.0, max_hours, Units::Plain),
            y_label: "bases",
            units: Units::Si,
            gap: Gap::Hold,
            curves: curves(&|g| {
                g.data
                    .time
                    .iter()
                    .map(|t| (t.hours, t.cumulative_bases as f64))
                    .collect()
            }),
        });
    }
    page.push_str("</div>");
}

/// Opens the wrapper of a reference table, which scrolls in a fixed height once the
/// table has many rows.
fn long_table_wrap(rows: usize) -> &'static str {
    if rows > MANY_GROUPS {
        r#"<div class="table-wrap long">"#
    } else {
        r#"<div class="table-wrap">"#
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn formats_thousands() {
        assert_eq!(thousands(0), "0");
        assert_eq!(thousands(999), "999");
        assert_eq!(thousands(1_234_567), "1,234,567");
    }
}
