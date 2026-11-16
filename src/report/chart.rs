//! Static SVG charts: lines, bands, bars and heatmap cells on linear or log axes.
//!
//! Colors come from CSS classes (`s0`..`s11`, `zone-*`) so the page's light and dark
//! themes restyle them; hover text uses SVG `<title>` elements, so no script is needed.
//!
//! Lines and bars carry no hover elements of their own. Each chart gets one row of
//! invisible full-height strips instead, each naming the values of every series at that x,
//! which keeps the markup small however many points or series a chart has.

use std::collections::BTreeMap;
use std::fmt::Write as _;

/// Chart size in SVG units; the page scales it to its container.
const WIDTH: f64 = 760.0;
/// Width of a chart that spans two columns of the page's chart grid: two charts and the
/// gap between them, so its text comes out the same size as theirs.
const WIDE_WIDTH: f64 = 2.0 * WIDTH + 26.0;
const HEIGHT: f64 = 300.0;
const LEFT: f64 = 66.0;
const RIGHT: f64 = 18.0;
const TOP: f64 = 14.0;
const BOTTOM: f64 = 46.0;

/// Bars at least this wide, in SVG units, are drawn with a gap of [`BAR_GAP`] after them.
const NARROW_BAR: f64 = 3.0;
const BAR_GAP: f64 = 0.6;

/// Narrowest hover strip on a categorical axis, where every category keeps its own strip.
const CATEGORY_STRIP: f64 = 2.0;

/// Half-width of a hover strip that has no neighbor to share the plot with.
const LONE_STRIP: f64 = 8.0;

/// How much hover and path detail a chart carries.
///
/// Reports with many groups render their charts `Compact` so the page stays small.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Detail {
    #[default]
    Full,
    Compact,
}

impl Detail {
    /// Narrowest hover strip, in SVG units; closer points share a strip.
    fn min_strip(self) -> f64 {
        match self {
            Self::Full => 6.0,
            Self::Compact => 13.0,
        }
    }

    /// Line points are dropped while the drawn line stays within this distance, in SVG
    /// units, of the full one.
    fn tolerance(self) -> f64 {
        match self {
            Self::Full => 0.3,
            Self::Compact => 0.6,
        }
    }

    /// Number of distinct heatmap colors. Fewer shades give longer runs of equal cells,
    /// which draw as one stroke.
    pub fn shades(self) -> u32 {
        match self {
            Self::Full => 48,
            Self::Compact => 16,
        }
    }
}

/// Number of series colors defined in the stylesheet.
pub const PALETTE: usize = 12;

/// Axis mapping.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Scale {
    Linear,
    Log10,
}

/// How tick labels are written.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Units {
    /// 1.5k, 20M.
    Si,
    /// Plain numbers with up to one decimal.
    Plain,
    Percent,
}

/// One chart axis.
#[derive(Debug, Clone)]
pub struct Axis {
    pub label: String,
    pub scale: Scale,
    pub units: Units,
    pub min: f64,
    pub max: f64,
    /// Category names at 0, 1, 2, ...; empty for numeric axes.
    pub categories: Vec<String>,
}

impl Axis {
    pub fn linear(label: &str, min: f64, max: f64, units: Units) -> Self {
        Self {
            label: label.to_owned(),
            scale: Scale::Linear,
            units,
            min,
            max: if max > min { max } else { min + 1.0 },
            categories: Vec::new(),
        }
    }

    /// A categorical axis: one slot per name, centered on 0, 1, 2, ...
    pub fn categories(label: &str, names: Vec<String>) -> Self {
        Self {
            label: label.to_owned(),
            scale: Scale::Linear,
            units: Units::Plain,
            min: -0.5,
            max: names.len().max(1) as f64 - 0.5,
            categories: names,
        }
    }

    pub fn log(label: &str, min: f64, max: f64, units: Units) -> Self {
        let min = min.max(1e-9);
        Self {
            label: label.to_owned(),
            scale: Scale::Log10,
            units,
            min,
            max: if max > min { max } else { min * 10.0 },
            categories: Vec::new(),
        }
    }

    /// Position of `value` as a fraction of the axis.
    fn fraction(&self, value: f64) -> f64 {
        match self.scale {
            Scale::Linear => (value - self.min) / (self.max - self.min),
            Scale::Log10 => {
                (value.max(self.min).log10() - self.min.log10())
                    / (self.max.log10() - self.min.log10())
            }
        }
    }

    fn tick_label(&self, tick: f64) -> String {
        if self.categories.is_empty() {
            format_tick(tick, self.units)
        } else {
            self.categories
                .get(tick as usize)
                .cloned()
                .unwrap_or_default()
        }
    }

    fn ticks(&self) -> Vec<f64> {
        if !self.categories.is_empty() {
            // Label at most about 40 categories.
            let step = self.categories.len().div_ceil(40);
            return (0..self.categories.len())
                .step_by(step)
                .map(|i| i as f64)
                .collect();
        }
        match self.scale {
            Scale::Linear => linear_ticks(self.min, self.max),
            Scale::Log10 => log_ticks(self.min, self.max),
        }
    }
}

/// Something drawn on a chart.
#[derive(Debug, Clone)]
pub enum Layer {
    /// Horizontal background bands (quality zones), as y ranges with a CSS class.
    Zones(Vec<(f64, f64, &'static str)>),
    /// Vertical background bands, as x ranges with a CSS class: quality zones of a chart
    /// with read Q along x.
    XZones(Vec<(f64, f64, &'static str)>),
    /// A filled band between two curves sharing x values.
    Band {
        series: usize,
        points: Vec<(f64, f64, f64)>,
        opacity: f64,
    },
    Line {
        name: String,
        series: usize,
        points: Vec<(f64, f64)>,
        style: LineStyle,
    },
    Bars {
        series: usize,
        bars: Vec<Bar>,
    },
    /// A dashed horizontal reference line across the plot, labeled at its right end.
    Level {
        name: String,
        series: usize,
        y: f64,
    },
    /// Heatmap cells, each with its own fill color.
    Cells(Vec<Cell>),
    /// Invisible rectangles with hover text, laid over heatmap cells. They may be coarser
    /// than the cells, so a detailed heatmap needs few titles.
    Hits(Vec<Hit>),
}

/// Line dash style; dashing distinguishes series once the palette repeats.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LineStyle {
    Solid,
    Dashed,
    Thin,
}

/// One bar, from `x0` to `x1` on the x axis and from 0 (or the axis floor) to `y`.
#[derive(Debug, Clone)]
pub struct Bar {
    pub x0: f64,
    pub x1: f64,
    pub y: f64,
    pub tip: String,
}

/// One heatmap cell.
#[derive(Debug, Clone)]
pub struct Cell {
    pub x0: f64,
    pub x1: f64,
    pub y0: f64,
    pub y1: f64,
    pub color: String,
}

/// One hover rectangle over a heatmap.
#[derive(Debug, Clone)]
pub struct Hit {
    pub x0: f64,
    pub x1: f64,
    pub y0: f64,
    pub y1: f64,
    pub tip: String,
}

/// A chart: axes plus layers drawn in order.
#[derive(Debug, Clone)]
pub struct Chart {
    pub title: String,
    pub x: Axis,
    pub y: Axis,
    pub layers: Vec<Layer>,
    pub height: f64,
    width: f64,
    /// Space below the plot for tick labels and the axis label.
    bottom: f64,
}

impl Chart {
    pub fn new(title: &str, x: Axis, y: Axis) -> Self {
        // Rotated category labels need room below the axis.
        let bottom = if x.categories.is_empty() {
            BOTTOM
        } else {
            let longest = x
                .categories
                .iter()
                .map(|c| c.chars().count())
                .max()
                .unwrap_or(0);
            (BOTTOM + longest as f64 * 4.0).min(150.0)
        };
        let height = HEIGHT + bottom - BOTTOM;
        Self {
            title: title.to_owned(),
            x,
            y,
            layers: Vec::new(),
            height,
            width: WIDTH,
            bottom,
        }
    }

    pub fn layer(mut self, layer: Layer) -> Self {
        self.layers.push(layer);
        self
    }

    /// Widens the chart to span two columns of the page's chart grid.
    pub fn wide(mut self) -> Self {
        self.width = WIDE_WIDTH;
        self
    }

    pub fn tall(mut self, height: f64) -> Self {
        self.height = height + self.bottom - BOTTOM;
        self
    }

    fn px(&self, x: f64) -> f64 {
        LEFT + self.x.fraction(x).clamp(0.0, 1.0) * (self.width - LEFT - RIGHT)
    }

    fn py(&self, y: f64) -> f64 {
        let plot = self.height - TOP - self.bottom;
        TOP + (1.0 - self.y.fraction(y).clamp(0.0, 1.0)) * plot
    }

    /// Renders the chart as an inline SVG element.
    pub fn render(&self, detail: Detail) -> String {
        let mut svg = String::with_capacity(16 * 1024);
        let h = self.height;
        let _ = write!(
            svg,
            r#"<svg class="chart" viewBox="0 0 {} {h}" role="img" aria-label="{}">"#,
            self.width,
            escape(&self.title)
        );
        let (x0, x1) = (LEFT, self.width - RIGHT);
        let (y0, y1) = (TOP, h - self.bottom);

        for layer in &self.layers {
            match layer {
                Layer::Zones(zones) => {
                    for &(lo, hi, class) in zones {
                        let (top, bottom) = (self.py(hi), self.py(lo));
                        let _ = write!(
                            svg,
                            r#"<rect class="{class}" x="{x0:.1}" y="{top:.1}" width="{:.1}" height="{:.1}"/>"#,
                            x1 - x0,
                            (bottom - top).max(0.0)
                        );
                    }
                }
                Layer::XZones(zones) => {
                    for &(lo, hi, class) in zones {
                        let (left, right) = (self.px(lo), self.px(hi));
                        let _ = write!(
                            svg,
                            r#"<rect class="{class}" x="{left:.1}" y="{y0:.1}" width="{:.1}" height="{:.1}"/>"#,
                            (right - left).max(0.0),
                            y1 - y0
                        );
                    }
                }
                _ => {}
            }
        }

        // Grid lines as one path, and tick labels grouped by their anchor.
        let mut grid = String::new();
        svg.push_str(r#"<g text-anchor="end">"#);
        for tick in self.y.ticks() {
            let y = self.py(tick);
            let _ = write!(grid, "M{},{}H{}", num(x0), num(y), num(x1));
            let _ = write!(
                svg,
                r#"<text x="{}" y="{}">{}</text>"#,
                num(x0 - 6.0),
                num(y + 4.0),
                escape(&self.y.tick_label(tick))
            );
        }
        let numeric = self.x.categories.is_empty();
        if numeric {
            svg.push_str(r#"</g><g text-anchor="middle">"#);
        }
        for tick in self.x.ticks() {
            let x = num(self.px(tick));
            let label = escape(&self.x.tick_label(tick));
            if numeric {
                let _ = write!(grid, "M{x},{}V{}", num(y0), num(y1));
                let _ = write!(
                    svg,
                    r#"<text x="{x}" y="{}">{label}</text>"#,
                    num(y1 + 16.0)
                );
            } else {
                let _ = write!(
                    svg,
                    r#"<text transform="translate({x} {}) rotate(-35)">{label}</text>"#,
                    num(y1 + 12.0)
                );
            }
        }
        let _ = write!(svg, r#"</g><path class="grid" d="{grid}"/>"#);

        for layer in &self.layers {
            match layer {
                Layer::Zones(_) | Layer::XZones(_) => {}
                Layer::Band {
                    series,
                    points,
                    opacity,
                } => self.band(&mut svg, *series, points, *opacity, detail),
                Layer::Line {
                    series,
                    points,
                    style,
                    ..
                } => self.line(&mut svg, *series, points, *style, detail),
                Layer::Bars { series, bars } => self.bars(&mut svg, *series, bars),
                Layer::Level { name, series, y } => self.level(&mut svg, name, *series, *y),
                Layer::Cells(cells) => self.cells(&mut svg, cells),
                Layer::Hits(hits) => self.hits(&mut svg, hits),
            }
        }

        let _ = write!(svg, r#"<path class="axis" d="M{x0},{y0}V{y1}H{x1}"/>"#);
        let _ = write!(
            svg,
            r#"<text class="label" x="{:.1}" y="{:.1}" text-anchor="middle">{}</text>"#,
            (x0 + x1) / 2.0,
            h - 8.0,
            escape(&self.x.label)
        );
        let _ = write!(
            svg,
            r#"<text class="label" transform="translate(14 {:.1}) rotate(-90)" text-anchor="middle">{}</text>"#,
            (y0 + y1) / 2.0,
            escape(&self.y.label)
        );
        self.hover(&mut svg, detail);
        svg.push_str("</svg>");
        svg
    }

    fn band(
        &self,
        svg: &mut String,
        series: usize,
        points: &[(f64, f64, f64)],
        opacity: f64,
        detail: Detail,
    ) {
        if points.is_empty() {
            return;
        }
        let edge = |y: fn(&(f64, f64, f64)) -> f64| -> Vec<(f64, f64)> {
            let pixels: Vec<(f64, f64)> = points
                .iter()
                .map(|p| (self.px(p.0), self.py(y(p))))
                .collect();
            simplify(&pixels, detail.tolerance())
        };
        let mut outline = edge(|p| p.1);
        outline.extend(edge(|p| p.2).into_iter().rev());
        let _ = write!(
            svg,
            r#"<path class="band s{}" d="{}z" fill-opacity="{opacity}"/>"#,
            series % PALETTE,
            polyline(&outline)
        );
    }

    fn line(
        &self,
        svg: &mut String,
        series: usize,
        points: &[(f64, f64)],
        style: LineStyle,
        detail: Detail,
    ) {
        if points.is_empty() {
            return;
        }
        let pixels: Vec<(f64, f64)> = points
            .iter()
            .map(|&(x, y)| (self.px(x), self.py(y)))
            .collect();
        let style_class = match style {
            LineStyle::Solid => "",
            LineStyle::Dashed => " dashed",
            LineStyle::Thin => " thin",
        };
        let _ = write!(
            svg,
            r#"<path class="line s{}{style_class}" d="{}"/>"#,
            series % PALETTE,
            polyline(&simplify(&pixels, detail.tolerance()))
        );
    }

    /// All bars of a layer as one path; their hover text lives in the hover strips.
    ///
    /// Bars that touch are drawn as one stepped outline, which is what a histogram of
    /// many narrow bins comes to.
    fn bars(&self, svg: &mut String, series: usize, bars: &[Bar]) {
        let base = num(self.py(self.y.min));
        let mut path = String::new();
        // Right edge of the bar drawn last, while its outline is still open.
        let mut open: Option<f64> = None;
        for bar in bars {
            let (left, right) = (self.px(bar.x0), self.px(bar.x1));
            let narrow = right - left < NARROW_BAR;
            // Wide bars keep a gap between them; narrow ones would lose their width to it.
            let right = (right - if narrow { 0.0 } else { BAR_GAP }).max(left + 0.8);
            let top = num(self.py(bar.y));
            if top == base {
                // Nothing to draw: an empty bin ends the outline before it.
                if open.take().is_some() {
                    let _ = write!(path, "V{base}z");
                }
                continue;
            }
            match open {
                // A narrow bar joins the outline of the bar it touches or nearly touches.
                Some(end) if narrow && (left - end).abs() < 1.0 => {
                    let _ = write!(path, "V{top}H{}", num(right));
                }
                _ => {
                    if open.is_some() {
                        let _ = write!(path, "V{base}z");
                    }
                    let _ = write!(path, "M{},{base}V{top}H{}", num(left), num(right));
                }
            }
            open = Some(right);
        }
        if open.is_some() {
            let _ = write!(path, "V{base}z");
        }
        let _ = write!(
            svg,
            r#"<path class="bars s{}" d="{path}"/>"#,
            series % PALETTE
        );
    }

    fn level(&self, svg: &mut String, name: &str, series: usize, y: f64) {
        let (x0, x1, y) = (LEFT, self.width - RIGHT, self.py(y));
        let _ = write!(
            svg,
            r#"<line class="level s{}" x1="{x0}" x2="{x1}" y1="{y:.1}" y2="{y:.1}"/><text class="halo" x="{:.1}" y="{:.1}" text-anchor="end">{}</text>"#,
            series % PALETTE,
            x1 - 4.0,
            y - 4.0,
            escape(name)
        );
    }

    /// Cells as one path per color and cell height: each run of touching cells in a row
    /// is a horizontal stroke as thick as the cells are tall, which takes far less markup
    /// than rectangles.
    fn cells(&self, svg: &mut String, cells: &[Cell]) {
        // Strokes by color and height, as (middle y, left, right) in tenths of a unit.
        let mut strokes: BTreeMap<(&str, i64), Vec<Stroke>> = BTreeMap::new();
        let tenths = |value: f64| (value * 10.0).round() as i64;
        for cell in cells {
            let (left, top) = (self.px(cell.x0), self.py(cell.y1));
            let right = self.px(cell.x1).max(left + 0.5);
            let height = (self.py(cell.y0) - top).max(0.5);
            strokes
                .entry((&cell.color, tenths(height)))
                .or_default()
                .push((tenths(top + height / 2.0), tenths(left), tenths(right)));
        }
        svg.push_str(r#"<g class="cells">"#);
        for ((color, height), mut row) in strokes {
            row.sort_unstable();
            let mut path = String::new();
            let mut run: Option<Stroke> = None;
            for (y, left, right) in row {
                match &mut run {
                    Some((run_y, _, run_right)) if *run_y == y && left - *run_right <= 1 => {
                        *run_right = right;
                    }
                    _ => {
                        push_stroke(&mut path, run);
                        run = Some((y, left, right));
                    }
                }
            }
            push_stroke(&mut path, run);
            let _ = write!(
                svg,
                r#"<path stroke="{color}" stroke-width="{}" d="{path}"/>"#,
                num(height as f64 / 10.0)
            );
        }
        svg.push_str("</g>");
    }

    fn hits(&self, svg: &mut String, hits: &[Hit]) {
        svg.push_str(r#"<g class="hits">"#);
        for hit in hits {
            let mut path = String::new();
            let (left, top) = (self.px(hit.x0), self.py(hit.y1));
            rect_path(&mut path, left, top, self.px(hit.x1), self.py(hit.y0));
            let _ = write!(
                svg,
                r#"<path d="{path}"><title>{}</title></path>"#,
                escape(&hit.tip)
            );
        }
        svg.push_str("</g>");
    }

    /// Invisible full-height strips whose titles give every line's and bar's value there.
    ///
    /// Points closer together than the detail's strip width share a strip, which then shows
    /// each line's point nearest its middle and each bar layer's tallest bar.
    fn hover(&self, svg: &mut String, detail: Detail) {
        let heatmap = |layer: &Layer| matches!(layer, Layer::Cells(_) | Layer::Hits(_));
        if self.layers.iter().any(heatmap) {
            // Heatmaps bring their own hover rectangles; strips would cover them.
            return;
        }
        let mut xs: Vec<f64> = Vec::new();
        for layer in &self.layers {
            match layer {
                Layer::Line { points, .. } => xs.extend(points.iter().map(|p| self.px(p.0))),
                Layer::Bars { bars, .. } => {
                    xs.extend(bars.iter().map(|b| self.px((b.x0 + b.x1) / 2.0)));
                }
                Layer::Zones(_)
                | Layer::XZones(_)
                | Layer::Band { .. }
                | Layer::Level { .. }
                | Layer::Cells(_)
                | Layer::Hits(_) => {}
            }
        }
        xs.sort_by(f64::total_cmp);
        let min_strip = if self.x.categories.is_empty() {
            detail.min_strip()
        } else {
            CATEGORY_STRIP
        };
        let mut anchors: Vec<f64> = Vec::new();
        for x in xs {
            if anchors.last().is_none_or(|&last| x - last >= min_strip) {
                anchors.push(x);
            }
        }
        let (Some(&first), Some(&last)) = (anchors.first(), anchors.last()) else {
            return;
        };
        // Strip i spans edges[i]..edges[i + 1]; inner edges lie midway between anchors.
        let outer = anchors
            .get(1)
            .map_or(LONE_STRIP, |next| (next - first) / 2.0);
        let mut edges = Vec::with_capacity(anchors.len() + 1);
        edges.push((first - outer).max(LEFT));
        edges.extend(anchors.windows(2).map(|pair| (pair[0] + pair[1]) / 2.0));
        edges.push((last + outer).min(self.width - RIGHT));
        let strip_of = |px: f64| edges[1..anchors.len()].partition_point(|&edge| edge <= px);

        let mut tips: Vec<Vec<Tip<'_>>> = vec![Vec::new(); anchors.len()];
        for layer in &self.layers {
            match layer {
                Layer::Line { name, points, .. } => {
                    let mut nearest: Vec<Option<(f64, (f64, f64))>> = vec![None; anchors.len()];
                    for &point in points {
                        let px = self.px(point.0);
                        let strip = strip_of(px);
                        let distance = (px - anchors[strip]).abs();
                        if nearest[strip].is_none_or(|(best, _)| distance < best) {
                            nearest[strip] = Some((distance, point));
                        }
                    }
                    for (strip, found) in nearest.into_iter().enumerate() {
                        if let Some((_, (x, y))) = found {
                            tips[strip].push(Tip::Point { name, x, y });
                        }
                    }
                }
                Layer::Bars { bars, .. } => {
                    let mut tallest: Vec<Option<&Bar>> = vec![None; anchors.len()];
                    for bar in bars {
                        let strip = strip_of(self.px((bar.x0 + bar.x1) / 2.0));
                        if tallest[strip].is_none_or(|best| bar.y > best.y) {
                            tallest[strip] = Some(bar);
                        }
                    }
                    for (strip, found) in tallest.into_iter().enumerate() {
                        if let Some(bar) = found {
                            tips[strip].push(Tip::Bar(&bar.tip));
                        }
                    }
                }
                Layer::Zones(_)
                | Layer::XZones(_)
                | Layer::Band { .. }
                | Layer::Level { .. }
                | Layer::Cells(_)
                | Layer::Hits(_) => {}
            }
        }

        let height = num(self.height - self.bottom - TOP);
        svg.push_str(r#"<g class="hover">"#);
        for (strip, tips) in tips.iter().enumerate() {
            if tips.is_empty() {
                continue;
            }
            let (left, right) = (edges[strip], edges[strip + 1]);
            let width = num(right - left);
            let _ = write!(
                svg,
                r#"<path d="M{},{TOP}h{width}v{height}h-{width}z"><title>{}</title></path>"#,
                num(left),
                escape(&self.tip_text(tips))
            );
        }
        svg.push_str("</g>");
    }

    /// The hover text of one strip, one line per series.
    fn tip_text(&self, tips: &[Tip<'_>]) -> String {
        let x_text = |x: f64| {
            if self.x.categories.is_empty() {
                format_value(x, self.x.units)
            } else {
                self.x.tick_label(x.round())
            }
        };
        let mut point_xs = tips.iter().filter_map(|tip| match tip {
            Tip::Point { x, .. } => Some(x_text(*x)),
            Tip::Bar(_) => None,
        });
        // Lines that share an x get it once, as a heading.
        let shared_x = point_xs
            .next()
            .filter(|first| tips.len() > 1 && point_xs.all(|x| x == *first));
        let mut text = String::new();
        for tip in tips {
            if let Tip::Bar(tip) = tip {
                text.push_str(tip);
                text.push('\n');
            }
        }
        if let Some(x) = &shared_x {
            let has_bar = tips.iter().any(|tip| matches!(tip, Tip::Bar(_)));
            if !has_bar {
                if self.x.label.is_empty() {
                    let _ = writeln!(text, "{x}");
                } else {
                    let _ = writeln!(text, "{}: {x}", self.x.label);
                }
            }
        }
        for tip in tips {
            if let Tip::Point { name, x, y } = tip {
                let y = format_value(*y, self.y.units);
                if shared_x.is_some() {
                    let _ = writeln!(text, "{name}: {y}");
                } else {
                    let _ = writeln!(text, "{name}: {} → {y}", x_text(*x));
                }
            }
        }
        text.truncate(text.trim_end().len());
        text
    }
}

/// One series' contribution to a hover strip.
#[derive(Debug, Clone, Copy)]
enum Tip<'a> {
    Point { name: &'a str, x: f64, y: f64 },
    Bar(&'a str),
}

/// Appends a closed rectangle, at least half a unit wide and tall, to an SVG path.
fn rect_path(path: &mut String, left: f64, top: f64, right: f64, bottom: f64) {
    let width = num((right - left).max(0.5));
    let _ = write!(
        path,
        "M{},{}h{width}v{}h-{width}z",
        num(left),
        num(top),
        num((bottom - top).max(0.5))
    );
}

/// A horizontal stroke as (middle y, left, right), in tenths of a unit.
type Stroke = (i64, i64, i64);

/// Appends a horizontal stroke to an SVG path.
fn push_stroke(path: &mut String, run: Option<Stroke>) {
    if let Some((y, left, right)) = run {
        let _ = write!(
            path,
            "M{},{}h{}",
            num(left as f64 / 10.0),
            num(y as f64 / 10.0),
            num((right - left) as f64 / 10.0)
        );
    }
}

/// The points of a polyline that matter: those left out lie within `tolerance` of the
/// line through the points kept (Ramer-Douglas-Peucker).
fn simplify(points: &[(f64, f64)], tolerance: f64) -> Vec<(f64, f64)> {
    if points.len() < 3 {
        return points.to_vec();
    }
    let mut keep = vec![false; points.len()];
    keep[0] = true;
    keep[points.len() - 1] = true;
    let mut spans = vec![(0, points.len() - 1)];
    while let Some((first, last)) = spans.pop() {
        let (a, b) = (points[first], points[last]);
        let (dx, dy) = (b.0 - a.0, b.1 - a.1);
        let length = dx.hypot(dy);
        let distance = |p: (f64, f64)| {
            if length < 1e-9 {
                (p.0 - a.0).hypot(p.1 - a.1)
            } else {
                (dy * (p.0 - a.0) - dx * (p.1 - a.1)).abs() / length
            }
        };
        let farthest = (first + 1..last)
            .map(|i| (i, distance(points[i])))
            .max_by(|x, y| x.1.total_cmp(&y.1));
        if let Some((split, distance)) = farthest {
            if distance > tolerance {
                keep[split] = true;
                spans.push((first, split));
                spans.push((split, last));
            }
        }
    }
    points
        .iter()
        .zip(keep)
        .filter_map(|(&point, kept)| kept.then_some(point))
        .collect()
}

/// SVG path data for a polyline: an absolute start, then relative steps in tenths.
fn polyline(points: &[(f64, f64)]) -> String {
    let tenths = |value: f64| (value * 10.0).round() as i64;
    let mut path = String::with_capacity(points.len() * 8);
    let mut last: Option<(i64, i64)> = None;
    let mut stepping = false;
    for &(x, y) in points {
        let (x, y) = (tenths(x), tenths(y));
        match last {
            None => {
                let _ = write!(path, "M{},{}", num(x as f64 / 10.0), num(y as f64 / 10.0));
            }
            Some((last_x, last_y)) => {
                let (dx, dy) = (x - last_x, y - last_y);
                // One `l` starts the run of relative steps; the rest follow implicitly.
                path.push(if stepping { ' ' } else { 'l' });
                stepping = true;
                let _ = write!(path, "{},{}", num(dx as f64 / 10.0), num(dy as f64 / 10.0));
            }
        }
        last = Some((x, y));
    }
    path
}

/// A coordinate with one decimal, without a trailing `.0`.
fn num(value: f64) -> String {
    let text = format!("{value:.1}");
    match text.strip_suffix(".0") {
        Some(whole) => whole.to_owned(),
        None => text,
    }
}

/// Round-number ticks spanning `min..=max`.
fn linear_ticks(min: f64, max: f64) -> Vec<f64> {
    let span = max - min;
    let raw = span / 6.0;
    let magnitude = 10f64.powf(raw.log10().floor());
    let step = [1.0, 2.0, 2.5, 5.0, 10.0]
        .iter()
        .map(|m| m * magnitude)
        .find(|&s| s >= raw)
        .unwrap_or(10.0 * magnitude);
    let mut ticks = Vec::new();
    let mut tick = (min / step).ceil() * step;
    while tick <= max + step * 1e-9 {
        ticks.push(if tick.abs() < step * 1e-9 { 0.0 } else { tick });
        tick += step;
    }
    ticks
}

/// 1, 2 and 5 times powers of ten within `min..=max`, thinned to whole decades when wide.
fn log_ticks(min: f64, max: f64) -> Vec<f64> {
    let decades = (max.log10() - min.log10()).max(0.0);
    let multipliers: &[f64] = if decades > 4.0 {
        &[1.0]
    } else {
        &[1.0, 2.0, 5.0]
    };
    let mut ticks = Vec::new();
    let mut power = 10f64.powf(min.log10().floor());
    while power <= max * 1.0001 {
        for m in multipliers {
            let tick = power * m;
            if tick >= min * 0.9999 && tick <= max * 1.0001 {
                ticks.push(tick);
            }
        }
        power *= 10.0;
    }
    ticks
}

fn format_tick(value: f64, units: Units) -> String {
    match units {
        Units::Si => si(value),
        Units::Plain => trim_number(value),
        Units::Percent => format!("{}%", trim_number(value)),
    }
}

fn format_value(value: f64, units: Units) -> String {
    match units {
        Units::Si => si(value),
        Units::Plain => format!("{value:.2}")
            .trim_end_matches('0')
            .trim_end_matches('.')
            .to_owned(),
        Units::Percent => format!("{value:.2}%"),
    }
}

fn trim_number(value: f64) -> String {
    if (value - value.round()).abs() < 1e-9 {
        format!("{value:.0}")
    } else {
        format!("{value:.1}")
    }
}

/// Short SI form: 950, 1.5k, 20M, 3.2G, 1.1T.
pub fn si(value: f64) -> String {
    let abs = value.abs();
    let (scaled, suffix) = if abs >= 1e12 {
        (value / 1e12, "T")
    } else if abs >= 1e9 {
        (value / 1e9, "G")
    } else if abs >= 1e6 {
        (value / 1e6, "M")
    } else if abs >= 1e3 {
        (value / 1e3, "k")
    } else {
        (value, "")
    };
    let digits = if scaled.abs() >= 100.0 || (scaled - scaled.round()).abs() < 0.05 {
        0
    } else {
        1
    };
    format!("{scaled:.digits$}{suffix}")
}

/// Escapes text for SVG and HTML.
pub fn escape(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    for c in text.chars() {
        match c {
            '&' => out.push_str("&amp;"),
            '<' => out.push_str("&lt;"),
            '>' => out.push_str("&gt;"),
            '"' => out.push_str("&quot;"),
            '\'' => out.push_str("&#39;"),
            _ => out.push(c),
        }
    }
    out
}

/// A color from a perceptually ordered ramp (dark blue to yellow, like viridis), for
/// `fraction` in 0..=1, rounded to one of `shades` evenly spaced colors (at least 2) so
/// that heatmap cells sharing a shade can be drawn together.
pub fn ramp(fraction: f64, shades: u32) -> String {
    const STOPS: [(f64, [f64; 3]); 5] = [
        (0.0, [68.0, 1.0, 84.0]),
        (0.25, [59.0, 82.0, 139.0]),
        (0.5, [33.0, 145.0, 140.0]),
        (0.75, [94.0, 201.0, 98.0]),
        (1.0, [253.0, 231.0, 37.0]),
    ];
    let last_shade = f64::from(shades.max(2) - 1);
    let f = (fraction.clamp(0.0, 1.0) * last_shade).round() / last_shade;
    let i = STOPS
        .iter()
        .rposition(|&(at, _)| at <= f)
        .unwrap_or(0)
        .min(STOPS.len() - 2);
    let (a, ca) = STOPS[i];
    let (b, cb) = STOPS[i + 1];
    let t = (f - a) / (b - a);
    let mix = |k: usize| (ca[k] + (cb[k] - ca[k]) * t).round() as u8;
    format!("#{:02x}{:02x}{:02x}", mix(0), mix(1), mix(2))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn formats_si() {
        assert_eq!(si(950.0), "950");
        assert_eq!(si(1500.0), "1.5k");
        assert_eq!(si(20_000_000.0), "20M");
        assert_eq!(si(153_260_000_000.0), "153G");
    }

    #[test]
    fn ticks_are_round_and_in_range() {
        assert_eq!(linear_ticks(0.0, 40.0), vec![0.0, 10.0, 20.0, 30.0, 40.0]);
        assert_eq!(log_ticks(100.0, 1000.0), vec![100.0, 200.0, 500.0, 1000.0]);
        assert_eq!(
            log_ticks(10.0, 1e7),
            vec![10.0, 100.0, 1e3, 1e4, 1e5, 1e6, 1e7]
        );
    }

    #[test]
    fn draws_x_zones_as_vertical_bands() {
        // Read Q along x from 0 to 40: fail below 15, warn to 20, pass above.
        let chart = Chart::new(
            "q",
            Axis::linear("read Q", 0.0, 40.0, Units::Plain),
            Axis::linear("% of reads", 0.0, 8.0, Units::Percent),
        )
        .layer(Layer::XZones(vec![
            (0.0, 15.0, "zone-bad"),
            (15.0, 20.0, "zone-warn"),
            (20.0, 40.0, "zone-good"),
        ]));
        let svg = chart.render(Detail::Full);
        let width = |class: &str| -> f64 {
            let at = svg
                .find(&format!(r#"class="{class}""#))
                .expect("band drawn");
            let rest = &svg[at..];
            let w = rest.find("width=\"").expect("width") + 7;
            rest[w..w + rest[w..].find('"').expect("quote")]
                .parse()
                .expect("number")
        };
        // Each band is as wide as its share of the x axis, not the whole plot.
        let (bad, warn, good) = (width("zone-bad"), width("zone-warn"), width("zone-good"));
        assert!((bad / good - 15.0 / 20.0).abs() < 0.01, "{bad} {good}");
        assert!((warn / good - 5.0 / 20.0).abs() < 0.01, "{warn} {good}");
    }

    #[test]
    fn renders_escaped_svg() {
        let chart = Chart::new(
            "a<b",
            Axis::linear("x", 0.0, 1.0, Units::Plain),
            Axis::linear("y", 0.0, 1.0, Units::Plain),
        )
        .layer(Layer::Line {
            name: "s&t".into(),
            series: 0,
            points: vec![(0.0, 0.0), (1.0, 1.0)],
            style: LineStyle::Solid,
        });
        let svg = chart.render(Detail::Full);
        assert!(svg.starts_with("<svg") && svg.ends_with("</svg>"));
        assert!(svg.contains("a&lt;b") && svg.contains("s&amp;t"));
    }

    #[test]
    fn simplify_keeps_corners_and_drops_straight_runs() {
        let flat: Vec<(f64, f64)> = (0..100).map(|i| (f64::from(i), 5.0)).collect();
        assert_eq!(simplify(&flat, 0.3), vec![(0.0, 5.0), (99.0, 5.0)]);
        let spike = [(0.0, 0.0), (1.0, 0.0), (2.0, 9.0), (3.0, 0.0), (4.0, 0.0)];
        assert_eq!(simplify(&spike, 0.3), spike.to_vec());
    }

    #[test]
    fn polyline_steps_are_relative() {
        let path = polyline(&[(10.0, 20.0), (12.5, 19.0), (15.0, 19.0)]);
        assert_eq!(path, "M10,20l2.5,-1 2.5,0");
    }

    #[test]
    fn touching_cells_of_one_color_share_a_stroke() {
        let cell = |x: f64, color: &str| Cell {
            x0: x,
            x1: x + 1.0,
            y0: 0.0,
            y1: 1.0,
            color: color.to_owned(),
        };
        let svg = Chart::new(
            "cells",
            Axis::linear("x", 0.0, 4.0, Units::Plain),
            Axis::linear("y", 0.0, 1.0, Units::Plain),
        )
        .layer(Layer::Cells(vec![
            cell(0.0, "#111111"),
            cell(1.0, "#111111"),
            cell(2.0, "#222222"),
            cell(3.0, "#111111"),
        ]))
        .render(Detail::Full);
        assert_eq!(svg.matches("<path stroke=").count(), 2);
        assert_eq!(svg.matches(",134h").count(), 3, "{svg}");
    }

    fn bars_chart(count: usize) -> Chart {
        let bars = (0..count)
            .map(|i| Bar {
                x0: i as f64,
                x1: i as f64 + 1.0,
                y: (i % 7) as f64,
                tip: format!("bar {i}"),
            })
            .collect();
        Chart::new(
            "bars",
            Axis::linear("x", 0.0, count as f64, Units::Plain),
            Axis::linear("y", 0.0, 7.0, Units::Plain),
        )
        .layer(Layer::Bars { series: 0, bars })
    }

    #[test]
    fn wide_bars_each_get_a_hover_strip() {
        let svg = bars_chart(20).render(Detail::Full);
        assert_eq!(svg.matches("<title>").count(), 20);
        assert!(svg.contains("<title>bar 19</title>"));
    }

    #[test]
    fn narrow_bars_share_strips_and_show_the_tallest() {
        let full = bars_chart(2000).render(Detail::Full);
        let compact = bars_chart(2000).render(Detail::Compact);
        let strips = |svg: &str| svg.matches("<title>").count();
        assert!(strips(&full) <= 113, "full strips: {}", strips(&full));
        assert!(strips(&compact) < strips(&full));
        // Bar 6 is the first of height 6, the tallest in the first strip.
        assert!(full.contains("<title>bar 6</title>"));
    }

    #[test]
    fn lines_sharing_x_get_one_heading() {
        let line = |name: &str, y: f64| Layer::Line {
            name: name.into(),
            series: 0,
            points: vec![(0.0, y), (10.0, y)],
            style: LineStyle::Solid,
        };
        let svg = Chart::new(
            "lines",
            Axis::linear("pos", 0.0, 10.0, Units::Plain),
            Axis::linear("y", 0.0, 5.0, Units::Plain),
        )
        .layer(line("a", 1.0))
        .layer(line("b", 2.0))
        .render(Detail::Full);
        assert!(svg.contains("<title>pos: 10\na: 1\nb: 2</title>"), "{svg}");
    }

    #[test]
    fn ramp_spans_its_stops() {
        assert_eq!(ramp(0.0, 48), "#440154");
        assert_eq!(ramp(1.0, 48), "#fde725");
        assert_eq!(ramp(0.4, 2), "#440154");
    }
}
