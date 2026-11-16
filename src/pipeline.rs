//! Reading the inputs on reader threads and counting their reads on worker threads.

use std::collections::{BTreeMap, HashMap};
use std::io::{BufRead, BufReader};
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};
use std::sync::Arc;
use std::thread;

use anyhow::{anyhow, bail, Context, Result};
use crossbeam_channel::{Receiver, Sender};
use serde::{Deserialize, Serialize};

use crate::adapters::Adapters;
use crate::barcode::{LabelCache, NO_BARCODE};
use crate::barcodes::Barcodes;
use crate::input::{self, Format, InputFile};
use crate::read::{Platform, Read};
use crate::screen::Screen;
use crate::source::Source;
use crate::stats::{Scratch, Settings, Stats};
use crate::{bam, bgzf, fastx, per_read, summary};

/// Bytes of records per batch handed to a worker.
const BATCH_BYTES: usize = 4 << 20;

/// Records per batch at most, for inputs of small records (summaries).
const BATCH_RECORDS: usize = 8192;

/// What the pipeline needs to know besides the files.
#[derive(Debug)]
pub struct Options {
    pub threads: usize,
    /// Platform forced on every read, overriding headers and tags.
    pub platform: Option<Platform>,
    pub barcode_names: Arc<Vec<String>>,
    pub max_reads: Option<u64>,
    pub min_q: Option<f32>,
    pub adapters: Adapters,
    pub barcodes: Barcodes,
    /// The contaminant screen; `None` turns it off.
    pub screen: Option<Screen>,
    pub kinnex: Option<Adapters>,
    pub group_by: GroupBy,
    /// Where to write the per-read table, if anywhere.
    pub per_read: Option<per_read::Output>,
}

/// How reads are grouped.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GroupBy {
    Barcode,
    /// Each run (flow cell or movie) and barcode apart: `PAY87954 · barcode01`.
    RunAndBarcode,
}

/// What reading the inputs produced.
#[derive(Debug)]
pub struct Collected {
    /// Statistics by group name.
    pub groups: BTreeMap<String, Stats>,
    pub files: Vec<FileInfo>,
}

impl Collected {
    /// Adds another pass's reads, for watch mode.
    pub fn merge(&mut self, other: Collected) {
        for (name, stats) in other.groups {
            match self.groups.get_mut(&name) {
                Some(existing) => existing.merge(stats),
                None => {
                    self.groups.insert(name, stats);
                }
            }
        }
        self.files.extend(other.files);
        self.files.sort_by(|a, b| a.path.cmp(&b.path));
    }
}

/// One input file, as the report describes it.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FileInfo {
    pub path: String,
    pub format: Format,
    pub records: u64,
    pub path_barcode: Option<String>,
    /// `@RG` and `@PG` header lines of BAM inputs.
    pub header_lines: Vec<String>,
}

/// Records of one file, ready for decoding.
#[derive(Debug)]
struct Batch {
    source: Arc<Source>,
    format: Format,
    path: Arc<str>,
    data: Vec<u8>,
    /// End offset of each record in `data`.
    ends: Vec<usize>,
    /// Position in the input, for writing the per-read table in order.
    file: usize,
    index: u64,
}

impl Batch {
    fn new(source: &Arc<Source>, format: Format, path: &Arc<str>, file: usize, index: u64) -> Self {
        Self {
            source: Arc::clone(source),
            format,
            path: Arc::clone(path),
            data: Vec::with_capacity(BATCH_BYTES + (1 << 16)),
            ends: Vec::new(),
            file,
            index,
        }
    }

    fn is_full(&self) -> bool {
        self.data.len() >= BATCH_BYTES || self.ends.len() >= BATCH_RECORDS
    }
}

/// Shared state of the reader threads.
#[derive(Debug)]
struct Readers<'run> {
    files: &'run [InputFile],
    options: &'run Options,
    next_file: AtomicUsize,
    records: AtomicU64,
}

/// Reads every file and returns per-group statistics.
///
/// # Errors
///
/// Fails on the first unreadable or malformed input.
pub fn run(files: &[InputFile], options: &Options) -> Result<Collected> {
    let settings = Settings {
        min_q: options.min_q,
        adapters: &options.adapters,
        barcodes: &options.barcodes,
        screen: options.screen.as_ref(),
        kinnex: options.kinnex.as_ref(),
    };
    let readers = Readers {
        files,
        options,
        next_file: AtomicUsize::new(0),
        records: AtomicU64::new(0),
    };
    // Decompression runs on the rayon pool, so a few readers keep the workers busy.
    let reader_count = files.len().clamp(1, (options.threads / 4).max(1));
    let (sender, receiver) = crossbeam_channel::bounded::<Batch>(options.threads * 2);
    // Out-of-order chunks wait in the writer, so the channel needs no bound.
    let (rows, row_receiver) = crossbeam_channel::unbounded::<per_read::Message>();
    let rows = options.per_read.as_ref().map(|_| rows);
    let screen_names: Vec<String> = options
        .screen
        .as_ref()
        .map(|s| s.references().iter().map(|r| r.key.clone()).collect())
        .unwrap_or_default();
    let names = per_read::Names {
        adapters: options.adapters.names(),
        screen: &screen_names,
    };

    let (file_infos, worker_groups) = thread::scope(|scope| {
        let writer = options.per_read.as_ref().map(|output| {
            let row_receiver = row_receiver.clone();
            scope.spawn(move || per_read::write(output, &row_receiver))
        });
        drop(row_receiver);
        let reader_handles: Vec<_> = (0..reader_count)
            .map(|_| {
                let sender = sender.clone();
                let rows = rows.clone();
                let readers = &readers;
                scope.spawn(move || readers.read_all(&sender, rows.as_ref()))
            })
            .collect();
        drop(sender);
        let worker_handles: Vec<_> = (0..options.threads)
            .map(|_| {
                let receiver = receiver.clone();
                let rows = rows.clone();
                let group_by = options.group_by;
                scope.spawn(move || work(&receiver, &settings, group_by, rows, names))
            })
            .collect();
        drop(rows);
        let mut file_infos = Vec::new();
        let mut first_error = None;
        for handle in reader_handles {
            match handle
                .join()
                .map_err(|_| anyhow!("a reader thread panicked"))
            {
                Ok(Ok(infos)) => file_infos.extend(infos),
                Ok(Err(e)) | Err(e) => {
                    first_error.get_or_insert(e);
                }
            }
        }
        let mut groups = Vec::new();
        for handle in worker_handles {
            match handle
                .join()
                .map_err(|_| anyhow!("a worker thread panicked"))
            {
                Ok(Ok(g)) => groups.push(g),
                Ok(Err(e)) | Err(e) => {
                    first_error.get_or_insert(e);
                }
            }
        }
        if let Some(writer) = writer {
            match writer
                .join()
                .map_err(|_| anyhow!("the per-read writer panicked"))
            {
                Ok(Ok(())) => {}
                Ok(Err(e)) | Err(e) => {
                    first_error.get_or_insert(e);
                }
            }
        }
        match first_error {
            Some(e) => Err(e),
            None => Ok((file_infos, groups)),
        }
    })?;

    let mut groups: BTreeMap<String, Stats> = BTreeMap::new();
    for worker in worker_groups {
        for (name, stats) in worker {
            match groups.get_mut(&*name) {
                Some(existing) => existing.merge(stats),
                None => {
                    groups.insert(name.to_string(), stats);
                }
            }
        }
    }
    let mut files = file_infos;
    files.sort_by(|a, b| a.path.cmp(&b.path));
    Ok(Collected { groups, files })
}

impl Readers<'_> {
    fn read_all(
        &self,
        sender: &Sender<Batch>,
        rows: Option<&Sender<per_read::Message>>,
    ) -> Result<Vec<FileInfo>> {
        let mut infos = Vec::new();
        loop {
            let index = self.next_file.fetch_add(1, Ordering::Relaxed);
            let Some(file) = self.files.get(index) else {
                break;
            };
            if self.limit_reached() {
                break;
            }
            let info = self
                .read_file(index, file, sender, rows)
                .with_context(|| format!("reading {}", file.path.display()))?;
            infos.push(info);
        }
        Ok(infos)
    }

    fn limit_reached(&self) -> bool {
        self.options
            .max_reads
            .is_some_and(|max| self.records.load(Ordering::Relaxed) >= max)
    }

    fn read_file(
        &self,
        file_index: usize,
        file: &InputFile,
        sender: &Sender<Batch>,
        rows: Option<&Sender<per_read::Message>>,
    ) -> Result<FileInfo> {
        let (mut reader, format) = match file.format {
            Some(format) => (
                BufReader::with_capacity(1 << 20, bgzf::open(&file.path)?),
                format,
            ),
            None => {
                let mut reader = BufReader::with_capacity(1 << 20, bgzf::open_stdin()?);
                let format = sniff(&mut reader)?;
                (reader, format)
            }
        };
        let mut source = Source {
            platform: self.options.platform.unwrap_or_default(),
            barcode: file.path_barcode.as_deref().map(Arc::from),
            filter: file.path_filter,
            barcode_names: Arc::clone(&self.options.barcode_names),
            run: input::run_from_path(&file.path).map(Arc::from),
            ..Source::default()
        };
        let mut info = FileInfo {
            path: file.path.display().to_string(),
            format,
            records: 0,
            path_barcode: file.path_barcode.clone(),
            header_lines: Vec::new(),
        };
        let mut fasta = fastx::FastaSplitter::default();
        match format {
            Format::Bam => {
                let header = bam::read_header(&mut reader, &self.options.barcode_names)?;
                if self.options.platform.is_none() {
                    source.platform = bam::header_platform(&header);
                }
                source.aligned = header.references > 0;
                source.read_groups = bam::read_group_map(&header);
                info.header_lines = header
                    .read_groups
                    .iter()
                    .map(|g| {
                        let fields: Vec<String> =
                            g.fields.iter().map(|(k, v)| format!("{k}:{v}")).collect();
                        format!("@RG\t{}", fields.join("\t"))
                    })
                    .chain(header.programs.iter().cloned())
                    .collect();
            }
            Format::Summary => {
                let mut line = String::new();
                reader.read_line(&mut line)?;
                source.columns = summary::Columns::from_header(&line)?;
            }
            Format::Fastq | Format::Fasta => {}
        }

        let source = Arc::new(source);
        let path: Arc<str> = info.path.as_str().into();
        let mut batches = 0u64;
        let mut batch = Batch::new(&source, format, &path, file_index, batches);
        loop {
            let more = match format {
                Format::Bam => bam::read_record(&mut reader, &mut batch.data)?,
                Format::Fastq => fastx::read_fastq_record(&mut reader, &mut batch.data)?,
                Format::Fasta => fasta.next_record(&mut reader, &mut batch.data)?,
                Format::Summary => {
                    let n = reader.read_until(b'\n', &mut batch.data)?;
                    let line_start = batch.data.len() - n;
                    if n > 0 && batch.data[line_start..].iter().all(u8::is_ascii_whitespace) {
                        batch.data.truncate(line_start);
                        continue;
                    }
                    n > 0
                }
            };
            if !more {
                break;
            }
            batch.ends.push(batch.data.len());
            info.records += 1;
            if batch.is_full() {
                let records = batch.ends.len() as u64;
                batches += 1;
                let next = Batch::new(&source, format, &path, file_index, batches);
                let full = std::mem::replace(&mut batch, next);
                sender.send(full).context("workers stopped early")?;
                if self.records.fetch_add(records, Ordering::Relaxed) + records
                    >= self.options.max_reads.unwrap_or(u64::MAX)
                {
                    break;
                }
            }
        }
        if !batch.ends.is_empty() {
            self.records
                .fetch_add(batch.ends.len() as u64, Ordering::Relaxed);
            sender.send(batch).context("workers stopped early")?;
            batches += 1;
        }
        if let Some(rows) = rows {
            let done = per_read::Message::FileDone {
                file: file_index,
                batches,
            };
            rows.send(done)
                .context("the per-read writer stopped early")?;
        }
        Ok(info)
    }
}

type Groups = HashMap<Arc<str>, Stats>;

/// Format of a stream from its first decompressed bytes.
fn sniff(reader: &mut impl BufRead) -> Result<Format> {
    let head = reader.fill_buf().context("reading standard input")?;
    Ok(match head.first() {
        None => bail!("standard input is empty"),
        _ if head.starts_with(b"BAM\x01") => Format::Bam,
        Some(b'@') => Format::Fastq,
        Some(b'>') => Format::Fasta,
        _ if head
            .split(|&b| b == b'\n')
            .next()
            .is_some_and(|line| line.split(|&b| b == b'\t').any(|field| field == b"read_id")) =>
        {
            Format::Summary
        }
        _ => bail!("standard input is not BAM, FASTQ, FASTA or a sequencing summary"),
    })
}

/// A read's run (when runs are split) and barcode.
type GroupKey = (Option<Arc<str>>, Option<Arc<str>>);

/// The name of a read's group.
#[derive(Debug, Default)]
struct GroupNames {
    combined: HashMap<GroupKey, Arc<str>>,
}

impl GroupNames {
    fn name(&mut self, read: &Read, group_by: GroupBy) -> Arc<str> {
        let barcode = read.barcode.clone();
        let run = match group_by {
            GroupBy::Barcode => None,
            GroupBy::RunAndBarcode => Some(read.run.clone()),
        };
        let key = (run.clone().flatten(), barcode.clone());
        if let Some(name) = self.combined.get(&key) {
            return Arc::clone(name);
        }
        let barcode = barcode.as_deref().unwrap_or(NO_BARCODE).to_owned();
        let name: Arc<str> = match run {
            None => barcode.into(),
            Some(run) => format!("{} · {barcode}", run.as_deref().unwrap_or("unknown run")).into(),
        };
        self.combined.insert(key, Arc::clone(&name));
        name
    }
}

/// Decodes and counts batches until the readers finish. After an error it keeps draining
/// the channel, so readers never block on a full channel.
fn work(
    receiver: &Receiver<Batch>,
    settings: &Settings<'_>,
    group_by: GroupBy,
    rows: Option<Sender<per_read::Message>>,
    names: per_read::Names<'_>,
) -> Result<Groups> {
    let mut table = String::new();
    let mut groups = Groups::new();
    let mut cache = LabelCache::default();
    let mut group_names = GroupNames::default();
    let mut read = Read::default();
    let mut scratch = Scratch::default();
    let adapters = settings.adapters.names().len();
    let mut error = None;
    for batch in receiver {
        if error.is_some() {
            continue;
        }
        table.clear();
        let clock = match batch.format {
            Format::Summary => per_read::Clock::RunRelative,
            Format::Bam | Format::Fastq | Format::Fasta => per_read::Clock::Unix,
        };
        let mut start = 0;
        for &end in &batch.ends {
            let record = &batch.data[start..end];
            start = end;
            let decoded = match batch.format {
                Format::Bam => bam::decode(record, &batch.source, &mut cache, &mut read),
                Format::Fastq | Format::Fasta => {
                    fastx::decode(record, &batch.source, &mut cache, &mut read)
                }
                Format::Summary => summary::decode(record, &batch.source, &mut cache, &mut read),
            };
            if let Err(e) = decoded {
                error = Some(e.context(format!("decoding a record of {}", batch.path)));
                break;
            }
            let key = group_names.name(&read, group_by);
            groups
                .entry(Arc::clone(&key))
                .or_insert_with(|| Stats::new(adapters))
                .add(&read, settings, &mut scratch);
            if rows.is_some() && scratch.outcome().counted {
                per_read::write_row(&mut table, &read, &key, scratch.outcome(), names, clock);
            }
        }
        if let Some(rows) = &rows {
            let mut data = Vec::new();
            bgzf::compress(table.as_bytes(), &mut data);
            let chunk = per_read::Message::Chunk {
                file: batch.file,
                batch: batch.index,
                data,
            };
            if rows.send(chunk).is_err() {
                error.get_or_insert_with(|| anyhow!("the per-read writer stopped early"));
            }
        }
    }
    match error {
        Some(e) => Err(e),
        None => Ok(groups),
    }
}
