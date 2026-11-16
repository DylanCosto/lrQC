//! The per-read table: one tab-separated row per read, BGZF-compressed, in input order.
//!
//! Workers format and compress the rows of each batch; a writer thread puts the compressed
//! chunks back in input order, so the table is the same at any thread count.

use std::collections::{BTreeMap, HashMap};
use std::fs::{File, OpenOptions};
use std::io::{BufWriter, Write as _};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use crossbeam_channel::Receiver;

use crate::bgzf;
use crate::read::{Alignment, DuplexStatus, Filter, Origin, Platform, Read};
use crate::stats::ReadOutcome;
use crate::time;

/// The table's columns.
pub const COLUMNS: [&str; 27] = [
    "read_id",
    "group",
    "run",
    "platform",
    "length",
    "read_q",
    "filter",
    "gc_pct",
    "channel",
    "mux",
    "start_time",
    "duration",
    "duplex",
    "split",
    "passes",
    "accuracy",
    "effective_coverage",
    "zmw",
    "barcode_score",
    "poly_a",
    "screen",
    "adapters",
    "barcodes_inside",
    "kinnex_adapters",
    "mod_calls",
    "mapq",
    "identity",
];

/// Where the table goes, and whether earlier passes (watch mode) are kept.
#[derive(Debug, Clone)]
pub struct Output {
    pub path: PathBuf,
    pub mode: Mode,
}

/// How the table file is opened.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Mode {
    Create,
    /// Add to the table of an earlier pass; gzip members concatenate.
    Append,
}

/// What a read's start time counts from.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Clock {
    /// Unix time (BAM and FASTQ tags): written as ISO 8601.
    Unix,
    /// Seconds from the start of the run (sequencing summaries).
    RunRelative,
}

/// Names the row formatter needs.
#[derive(Debug, Clone, Copy)]
pub struct Names<'run> {
    pub adapters: &'run [String],
    pub screen: &'run [String],
}

/// Appends one read's row.
pub fn write_row(
    out: &mut String,
    read: &Read,
    group: &str,
    outcome: &ReadOutcome,
    names: Names<'_>,
    clock: Clock,
) {
    const NA: &str = "NA";
    let opt = |value: Option<String>| value.unwrap_or_else(|| NA.to_owned());
    let platform = match read.platform {
        Platform::Ont => "ont",
        Platform::PacBio => "pacbio",
        Platform::Unknown => NA,
    };
    let filter = match outcome.filter {
        Some(Filter::Pass) => "pass",
        Some(Filter::Fail) => "fail",
        None => NA,
    };
    let duplex = match read.duplex {
        Some(DuplexStatus::Simplex) => "simplex",
        Some(DuplexStatus::Duplex) => "duplex",
        Some(DuplexStatus::SimplexWithDuplexOffspring) => "simplex_with_duplex_offspring",
        None => NA,
    };
    let start = read.start_secs.map(|secs| match clock {
        Clock::Unix => time::format_iso8601(secs),
        Clock::RunRelative => format!("{secs:.3}"),
    });
    let (mapq, identity) = match read.alignment {
        Alignment::Mapped(m) => (
            Some(m.mapq.to_string()),
            outcome.identity.map(|i| format!("{:.4}", i * 100.0)),
        ),
        Alignment::Unmapped => (Some("unmapped".to_owned()), None),
        Alignment::NotAligned => (None, None),
    };
    let list = |items: Vec<&str>| {
        if items.is_empty() {
            String::new()
        } else {
            items.join(",")
        }
    };
    let mut adapters: Vec<&str> = outcome
        .adapters
        .iter()
        .map(|&i| names.adapters[i].as_str())
        .collect();
    adapters.dedup();
    let mut barcodes: Vec<&str> = outcome.barcodes.iter().map(String::as_str).collect();
    barcodes.dedup();

    let fields = [
        clean(&String::from_utf8_lossy(&read.name)),
        clean(group),
        opt(read.run.as_deref().map(clean)),
        platform.to_owned(),
        read.len.to_string(),
        opt(outcome.read_q.map(|q| format!("{q:.2}"))),
        filter.to_owned(),
        opt(outcome.gc_pct.map(|g| g.to_string())),
        opt(read.channel.map(|c| c.to_string())),
        opt(read.mux.map(|m| m.to_string())),
        opt(start),
        opt(read.duration_secs.map(|d| format!("{d:.3}"))),
        duplex.to_owned(),
        u8::from(read.origin == Origin::SplitFromParent).to_string(),
        opt(read.passes.map(|p| p.to_string())),
        opt(read.accuracy.map(|a| format!("{a:.6}"))),
        opt(read.effective_coverage.map(|e| format!("{e:.2}"))),
        opt(read.zmw.map(|z| z.to_string())),
        opt(read.barcode_score.map(|b| b.to_string())),
        opt(read.poly_a.map(|p| p.to_string())),
        // NA when not screened (no sequence, --no-screen); empty when nothing matched.
        match (outcome.screened, outcome.screen) {
            (false, _) => "NA".to_owned(),
            (true, None) => String::new(),
            (true, Some(i)) => names.screen[i].clone(),
        },
        clean(&list(adapters)),
        list(barcodes),
        opt(outcome.kinnex_adapters.map(|k| k.to_string())),
        read.mod_probs.len().to_string(),
        opt(mapq),
        opt(identity),
    ];
    for (i, field) in fields.iter().enumerate() {
        if i > 0 {
            out.push('\t');
        }
        out.push_str(field);
    }
    out.push('\n');
}

/// Read names and labels go in a TSV cell: tabs and newlines become spaces.
fn clean(text: &str) -> String {
    text.replace(['\t', '\n', '\r'], " ")
}

/// A message to the writer.
#[derive(Debug)]
pub enum Message {
    /// Compressed rows of one batch: file index, batch index within the file, BGZF bytes.
    Chunk {
        file: usize,
        batch: u64,
        data: Vec<u8>,
    },
    /// A file is fully read, in this many batches.
    FileDone { file: usize, batches: u64 },
}

/// Receives compressed chunks in any order and writes them in input order.
pub fn write(output: &Output, receiver: &Receiver<Message>) -> Result<()> {
    let exists = output.path.metadata().is_ok_and(|m| m.len() > 0);
    let file = match output.mode {
        Mode::Append => OpenOptions::new()
            .create(true)
            .append(true)
            .open(&output.path),
        Mode::Create => File::create(&output.path),
    }
    .with_context(|| format!("opening {}", output.path.display()))?;
    let mut out = BufWriter::with_capacity(1 << 20, file);
    let mut writer = Ordered::default();
    if output.mode == Mode::Create || !exists {
        let mut header = COLUMNS.join("\t");
        header.push('\n');
        let mut block = Vec::new();
        bgzf::compress(header.as_bytes(), &mut block);
        out.write_all(&block)?;
    }
    for message in receiver {
        writer.accept(message);
        for chunk in writer.ready() {
            out.write_all(&chunk)
                .with_context(|| format!("writing {}", output.path.display()))?;
        }
    }
    for chunk in writer.rest() {
        out.write_all(&chunk)?;
    }
    out.write_all(&bgzf::EOF_BLOCK)?;
    out.flush()
        .with_context(|| format!("writing {}", output.path.display()))
}

/// Reorders chunks: file 0's batches in order, then file 1's, and so on.
#[derive(Debug, Default)]
struct Ordered {
    pending: BTreeMap<(usize, u64), Vec<u8>>,
    finished: HashMap<usize, u64>,
    next_file: usize,
    next_batch: u64,
}

impl Ordered {
    fn accept(&mut self, message: Message) {
        match message {
            Message::Chunk { file, batch, data } => {
                self.pending.insert((file, batch), data);
            }
            Message::FileDone { file, batches } => {
                self.finished.insert(file, batches);
            }
        }
    }

    /// Chunks that can be written now.
    fn ready(&mut self) -> Vec<Vec<u8>> {
        let mut out = Vec::new();
        loop {
            if let Some(data) = self.pending.remove(&(self.next_file, self.next_batch)) {
                out.push(data);
                self.next_batch += 1;
            } else if self.finished.get(&self.next_file) == Some(&self.next_batch) {
                self.next_file += 1;
                self.next_batch = 0;
            } else {
                return out;
            }
        }
    }

    /// Everything left, in order: after an error or an early stop not every file finishes.
    fn rest(&mut self) -> Vec<Vec<u8>> {
        std::mem::take(&mut self.pending).into_values().collect()
    }
}

/// The table's path in an output directory.
pub fn path_in(dir: &Path) -> PathBuf {
    dir.join("reads.tsv.gz")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn writes_chunks_in_input_order() {
        let mut ordered = Ordered::default();
        let chunk = |file, batch| Message::Chunk {
            file,
            batch,
            data: vec![
                u8::try_from(file * 10).expect("small") + u8::try_from(batch).expect("small"),
            ],
        };
        ordered.accept(chunk(1, 0));
        ordered.accept(chunk(0, 1));
        assert!(ordered.ready().is_empty(), "file 0 batch 0 is missing");
        ordered.accept(chunk(0, 0));
        assert_eq!(ordered.ready(), vec![vec![0], vec![1]]);
        ordered.accept(Message::FileDone {
            file: 0,
            batches: 2,
        });
        assert_eq!(ordered.ready(), vec![vec![10]]);
        ordered.accept(chunk(1, 1));
        assert_eq!(ordered.ready(), vec![vec![11]]);
    }

    #[test]
    fn formats_rows() {
        let read = Read {
            name: b"r1".to_vec(),
            len: 1000,
            platform: Platform::Ont,
            start_secs: time::parse_iso8601(b"2024-08-20T10:39:34.457Z"),
            channel: Some(42),
            duplex: Some(DuplexStatus::Duplex),
            ..Read::default()
        };
        let outcome = ReadOutcome {
            counted: true,
            read_q: Some(18.25),
            filter: Some(Filter::Pass),
            gc_pct: Some(41),
            screen: Some(1),
            screened: true,
            adapters: vec![0, 0],
            ..ReadOutcome::default()
        };
        let adapters = ["LSK".to_owned()];
        let screen = ["lambda".to_owned(), "ecoli".to_owned()];
        let mut row = String::new();
        write_row(
            &mut row,
            &read,
            "barcode01",
            &outcome,
            Names {
                adapters: &adapters,
                screen: &screen,
            },
            Clock::Unix,
        );
        let fields: Vec<&str> = row.trim_end().split('\t').collect();
        assert_eq!(fields.len(), COLUMNS.len());
        let field = |name: &str| fields[COLUMNS.iter().position(|c| *c == name).expect("column")];
        assert_eq!(field("read_id"), "r1");
        assert_eq!(field("read_q"), "18.25");
        assert_eq!(field("start_time"), "2024-08-20T10:39:34.457Z");
        assert_eq!(field("duplex"), "duplex");
        assert_eq!(field("screen"), "ecoli");
        assert_eq!(field("adapters"), "LSK");
        assert_eq!(field("identity"), "NA");
    }
}
