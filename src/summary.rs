//! MinKNOW and dorado sequencing summaries: one tab-separated row per read, no sequence.

use anyhow::{bail, Context, Result};

use crate::barcode::LabelCache;
use crate::read::{DuplexStatus, Filter, Origin, Platform, Read};
use crate::source::Source;

/// Positions of the columns lrqc reads.
#[derive(Debug, Clone, Default)]
pub struct Columns {
    read_id: usize,
    length: usize,
    qscore: Option<usize>,
    passes_filtering: Option<usize>,
    channel: Option<usize>,
    start_time: Option<usize>,
    duration: Option<usize>,
    barcode: Option<usize>,
    end_reason: Option<usize>,
    run_id: Option<usize>,
    mux: Option<usize>,
    /// dorado duplex summaries: the template parent of a duplex read, `-` otherwise.
    duplex_parent: Option<usize>,
    /// The read a split read was cut from; the read itself when whole.
    parent_read_id: Option<usize>,
}

impl Columns {
    /// Finds the columns in a summary's header line.
    ///
    /// # Errors
    ///
    /// Fails when the header lacks `read_id` or `sequence_length_template`.
    pub fn from_header(line: &str) -> Result<Self> {
        let names: Vec<&str> = line.trim_end().split('\t').collect();
        let find = |name: &str| names.iter().position(|n| *n == name);
        let Some(read_id) = find("read_id") else {
            bail!("sequencing summary header has no read_id column");
        };
        Ok(Self {
            read_id,
            length: find("sequence_length_template")
                .context("sequencing summary header has no sequence_length_template column")?,
            qscore: find("mean_qscore_template"),
            passes_filtering: find("passes_filtering"),
            channel: find("channel"),
            start_time: find("start_time"),
            duration: find("duration"),
            barcode: find("barcode_arrangement"),
            end_reason: find("end_reason"),
            run_id: find("run_id"),
            mux: find("mux"),
            duplex_parent: find("duplex_parent_template"),
            parent_read_id: find("parent_read_id"),
        })
    }
}

/// Decodes one summary row.
///
/// # Errors
///
/// Fails when the row has fewer columns than the header or an unreadable length.
pub fn decode(line: &[u8], source: &Source, cache: &mut LabelCache, read: &mut Read) -> Result<()> {
    read.clear();
    let line = line.strip_suffix(b"\n").unwrap_or(line);
    let line = line.strip_suffix(b"\r").unwrap_or(line);
    let fields: Vec<&[u8]> = line.split(|&b| b == b'\t').collect();
    let columns = &source.columns;
    let field = |i: usize| fields.get(i).copied();
    let opt = |i: Option<usize>| i.and_then(field);

    let id = field(columns.read_id).context("summary row is shorter than its header")?;
    read.name.extend_from_slice(id);
    read.name_hash = crate::read::hash_name(id);
    read.len = field(columns.length)
        .and_then(parse)
        .context("summary row has no readable sequence_length_template")?;
    read.platform = Platform::Ont;
    read.barcode = source.barcode.clone();
    read.reported_q = opt(columns.qscore).and_then(parse);
    read.filter = match opt(columns.passes_filtering) {
        Some(b"TRUE" | b"True" | b"true" | b"1") => Some(Filter::Pass),
        Some(b"FALSE" | b"False" | b"false" | b"0") => Some(Filter::Fail),
        _ => source.filter,
    };
    read.channel = opt(columns.channel).and_then(parse);
    read.start_secs = opt(columns.start_time).and_then(parse);
    read.duration_secs = opt(columns.duration).and_then(parse);
    read.end_reason = opt(columns.end_reason)
        .filter(|r| !r.is_empty())
        .map(|r| cache.intern_verbatim(r));
    read.mux = opt(columns.mux).and_then(parse);
    read.run = opt(columns.run_id)
        .filter(|r| !r.is_empty())
        .map(|r| cache.intern_verbatim(&r[..r.len().min(8)]));
    if let Some(label) = opt(columns.barcode).filter(|b| !b.is_empty()) {
        read.barcode = Some(cache.intern(label));
    }
    // A duplex read is named `template;complement` and names its template parent. Its
    // simplex parents are rows of their own, which the summary does not mark, so they
    // count as simplex: the duplex rate (duplex bases against all simplex bases) holds.
    if columns.duplex_parent.is_some() {
        let duplex = opt(columns.duplex_parent).is_some_and(|p| !p.is_empty() && p != b"-")
            || id.contains(&b';');
        read.duplex = Some(if duplex {
            DuplexStatus::Duplex
        } else {
            DuplexStatus::Simplex
        });
    }
    if opt(columns.parent_read_id).is_some_and(|p| !p.is_empty() && p != b"-" && p != id) {
        read.origin = Origin::SplitFromParent;
    }
    read.sanitize();
    Ok(())
}

fn parse<T: std::str::FromStr>(bytes: &[u8]) -> Option<T> {
    std::str::from_utf8(bytes).ok()?.trim().parse().ok()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn decodes_minknow_rows() {
        let header = "filename\tread_id\trun_id\tchannel\tmux\tstart_time\tduration\tpasses_filtering\tsequence_length_template\tmean_qscore_template\tbarcode_arrangement\tend_reason";
        let source = Source {
            columns: Columns::from_header(header).expect("valid header"),
            ..Source::default()
        };
        let row = b"f.pod5\tabc\trun\t42\t1\t3600.5\t12.0\tTRUE\t4500\t17.25\tbarcode09\tsignal_positive\n";
        let mut read = Read::default();
        decode(row, &source, &mut LabelCache::default(), &mut read).expect("valid row");
        assert_eq!(read.len, 4500);
        assert_eq!(read.channel, Some(42));
        assert_eq!(read.filter, Some(Filter::Pass));
        assert_eq!(read.barcode.as_deref(), Some("barcode09"));
        assert_eq!(read.end_reason.as_deref(), Some("signal_positive"));
        assert_eq!(read.reported_q, Some(17.25));
    }

    #[test]
    fn reads_duplex_and_split_rows() {
        let header = "read_id\tparent_read_id\tduplex_parent_template\tsequence_length_template\tmean_qscore_template";
        let source = Source {
            columns: Columns::from_header(header).expect("valid header"),
            ..Source::default()
        };
        let mut read = Read::default();
        let mut cache = LabelCache::default();
        decode(b"a;b\ta;b\ta\t9000\t30.1\n", &source, &mut cache, &mut read).expect("valid");
        assert_eq!(read.duplex, Some(DuplexStatus::Duplex));
        decode(b"a\ta\t-\t9100\t20.5\n", &source, &mut cache, &mut read).expect("valid");
        assert_eq!(read.duplex, Some(DuplexStatus::Simplex));
        assert_eq!(read.origin, Origin::Whole);
        decode(b"c\tp\t-\t4000\tnan\n", &source, &mut cache, &mut read).expect("valid");
        assert_eq!(read.origin, Origin::SplitFromParent);
        assert_eq!(read.reported_q, None, "NaN is no quality");
    }

    #[test]
    fn rejects_headers_without_lengths() {
        assert!(Columns::from_header("read_id\tchannel").is_err());
    }
}
