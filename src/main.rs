//! lrqc: per-barcode quality control reports for Oxford Nanopore and PacBio long reads.

mod adapters;
mod bam;
mod barcode;
mod barcodes;
mod bgzf;
mod bins;
mod fastx;
mod input;
mod per_read;
mod pipeline;
mod read;
mod report;
mod sampler;
mod screen;
mod source;
mod stats;
mod summary;
mod time;

use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::fs;
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::sync::Arc;
use std::thread;
use std::time::{Duration, Instant, SystemTime};

use anyhow::{bail, Context, Result};
use clap::{Args, Parser, Subcommand, ValueEnum};

use crate::adapters::Adapters;
use crate::barcodes::Barcodes;
use crate::input::Format;
use crate::pipeline::{Collected, GroupBy, Options};
use crate::read::Platform;
use crate::report::checks::{ModuleId, Status};
use crate::report::{Report, ReportOptions, Thresholds};
use crate::screen::Screen;

/// Per-barcode QC reports for Oxford Nanopore and PacBio long reads.
///
/// Reads BAM (aligned or not), FASTQ and FASTA files, gzipped or not, MinKNOW output
/// directories, sequencing summaries and standard input (`-`), and writes report.html,
/// summary.tsv, qc.json and MultiQC files.
#[derive(Debug, Parser)]
#[command(version, about, args_conflicts_with_subcommands = true)]
struct Cli {
    #[command(subcommand)]
    command: Option<Command>,

    #[command(flatten)]
    run: RunArgs,
}

#[derive(Debug, Subcommand)]
enum Command {
    /// Compare the qc.json files of earlier runs in one report.
    Compare(CompareArgs),
    /// Build assets/human.sketch.gz from a human genome FASTA (for developers).
    #[command(hide = true)]
    BuildHumanSketch {
        /// The genome, such as T2T-CHM13v2.0's `GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz`.
        fasta: PathBuf,
        /// Where to write the sketch.
        output: PathBuf,
    },
}

#[derive(Debug, Args)]
struct RunArgs {
    /// Read files, directories to search for them, or `-` for standard input.
    #[arg(required_unless_present = "print_thresholds")]
    inputs: Vec<PathBuf>,

    /// Output directory.
    #[arg(short, long, default_value = "lrqc")]
    output: PathBuf,

    /// Worker threads.
    #[arg(short, long, default_value_t = default_threads())]
    threads: usize,

    /// Platform of the reads; auto reads it from headers and tags.
    #[arg(long, value_enum, default_value_t = PlatformArg::Auto)]
    platform: PlatformArg,

    /// Read Q at which reads pass when the input does not say, and the read and base
    /// quality warn level [default: 10 ONT, 20 PacBio].
    #[arg(long)]
    min_qscore: Option<f32>,

    /// Haploid genome size (3.1g, 4.6m, 48502), to report coverage.
    #[arg(long)]
    genome_size: Option<String>,

    /// MinKNOW or SMRT Link sample sheet (CSV), to name barcodes by sample.
    #[arg(long)]
    sample_sheet: Option<PathBuf>,

    /// FASTA of the barcodes given to lima: names `bc` tag indices, and its sequences find
    /// PacBio barcodes inside reads [default: PacBio's SMRTbell adapter indexes, for the
    /// barcodes the file names or sample sheet name].
    #[arg(long)]
    barcode_names: Option<PathBuf>,

    /// FASTA of extra adapter or primer sequences to search for.
    #[arg(long)]
    adapters: Option<PathBuf>,

    /// FASTA of Kinnex/MAS segmentation adapters (skera's), to count array segments.
    #[arg(long)]
    kinnex_adapters: Option<PathBuf>,

    /// TOML file overriding warn and fail levels (see --print-thresholds).
    #[arg(long)]
    thresholds: Option<PathBuf>,

    /// Print the default thresholds as TOML and exit.
    #[arg(long)]
    print_thresholds: bool,

    /// Exit with status 2 when any check warns, or fails.
    #[arg(long, value_enum)]
    fail_on: Option<FailOn>,

    /// Report each run (flow cell or movie) and barcode apart.
    #[arg(long)]
    split_runs: bool,

    /// Skip the screen for controls and contaminants.
    #[arg(long)]
    no_screen: bool,

    /// Also write reads.tsv.gz: one row per read with everything measured about it.
    #[arg(long)]
    per_read: bool,

    /// Watch the inputs during a run: every this many seconds, read new files and rewrite
    /// the report, until MinKNOW writes its final summary.
    #[arg(long, value_name = "SECONDS")]
    watch: Option<u64>,

    /// Report title [default: the first input's name].
    #[arg(long)]
    title: Option<String>,

    /// Stop after about this many reads, for a quick look.
    #[arg(long)]
    max_reads: Option<u64>,
}

#[derive(Debug, Args)]
struct CompareArgs {
    /// qc.json files, or directories holding one.
    #[arg(required = true)]
    reports: Vec<PathBuf>,

    /// Output directory.
    #[arg(short, long, default_value = "lrqc-compare")]
    output: PathBuf,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
enum PlatformArg {
    Auto,
    Ont,
    Pacbio,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
enum FailOn {
    Warn,
    Fail,
}

/// Kinnex adapters are 16 bases on HiFi reads, so they are matched exactly.
const KINNEX_MIN_SHARE: f64 = 1.0;

/// Exit status when a check reaches the `--fail-on` level; 1 is reserved for errors.
const QC_FAILED: u8 = 2;

fn default_threads() -> usize {
    thread::available_parallelism().map_or(4, usize::from)
}

fn main() -> Result<ExitCode> {
    // A usage error exits 1, like any other failure: 2 means a check reached --fail-on.
    let cli = match Cli::try_parse() {
        Ok(cli) => cli,
        Err(e) if e.use_stderr() => {
            let _ = e.print();
            return Ok(ExitCode::FAILURE);
        }
        Err(e) => {
            let _ = e.print();
            return Ok(ExitCode::SUCCESS);
        }
    };
    match cli.command {
        Some(Command::Compare(args)) => compare(&args),
        Some(Command::BuildHumanSketch { fasta, output }) => {
            screen::write_human_sketch(&fasta, &output)?;
            Ok(ExitCode::SUCCESS)
        }
        None => run(&cli.run),
    }
}

/// Barcode names the run's file paths and sample sheet give, one per barcode of a pair:
/// `m84011.bc2001--bc2002.bam` names `bc2001` and `bc2002`.
fn run_barcode_names<'a>(
    files: impl IntoIterator<Item = &'a input::InputFile>,
    samples: &HashMap<String, String>,
) -> BTreeSet<String> {
    files
        .into_iter()
        .filter_map(|f| f.path_barcode.clone())
        .chain(samples.keys().cloned())
        .flat_map(|label| label.split("--").map(str::to_owned).collect::<Vec<_>>())
        .collect()
}

/// The barcode index, with PacBio's built-in barcodes for `names`.
fn built_in_barcodes(names: &BTreeSet<String>) -> Barcodes {
    Barcodes::new(&barcodes::pacbio_built_in(names.iter().map(String::as_str)))
}

fn run(args: &RunArgs) -> Result<ExitCode> {
    if args.print_thresholds {
        print!("{}", Thresholds::default().to_toml());
        return Ok(ExitCode::SUCCESS);
    }
    let started = Instant::now();
    let threads = args.threads.max(1);
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .context("starting the thread pool")?;

    let barcode_fasta = match &args.barcode_names {
        Some(path) => adapters::read_fasta(path)?,
        None => Vec::new(),
    };
    let samples = match &args.sample_sheet {
        Some(path) => report::samples::load(path)?,
        None => HashMap::new(),
    };
    // Watching finds its files, and the barcodes they name, pass by pass.
    let files = match args.watch {
        Some(_) => Vec::new(),
        None => input::discover(&args.inputs)?,
    };
    let barcodes = if barcode_fasta.is_empty() {
        built_in_barcodes(&run_barcode_names(&files, &samples))
    } else {
        Barcodes::new(&barcode_fasta)
    };
    let mut options = Options {
        threads,
        platform: match args.platform {
            PlatformArg::Auto => None,
            PlatformArg::Ont => Some(Platform::Ont),
            PlatformArg::Pacbio => Some(Platform::PacBio),
        },
        barcode_names: Arc::new(barcode_fasta.iter().map(|(name, _)| name.clone()).collect()),
        max_reads: args.max_reads,
        min_q: args.min_qscore,
        adapters: Adapters::load(args.adapters.as_deref())?,
        barcodes,
        screen: if args.no_screen {
            None
        } else {
            Some(Screen::load()?)
        },
        kinnex: match &args.kinnex_adapters {
            Some(path) => Some(Adapters::from_sequences(
                adapters::read_fasta(path)?,
                KINNEX_MIN_SHARE,
            )),
            None => None,
        },
        group_by: if args.split_runs {
            GroupBy::RunAndBarcode
        } else {
            GroupBy::Barcode
        },
        per_read: per_read_output(args)?,
    };
    let report_options = ReportOptions {
        title: args
            .title
            .clone()
            .unwrap_or_else(|| default_title(&args.inputs)),
        command: std::env::args().collect::<Vec<_>>().join(" "),
        thresholds: match &args.thresholds {
            Some(path) => Thresholds::load(path)?,
            None => Thresholds::default(),
        },
        genome_size: args
            .genome_size
            .as_deref()
            .map(report::parse_genome_size)
            .transpose()?,
        samples,
        refresh_secs: args.watch,
    };

    let report = match args.watch {
        Some(interval) => watch(args, &mut options, &report_options, interval)?,
        None => {
            eprintln!(
                "lrqc: reading {} input(s) with {threads} threads",
                files.len()
            );
            let collected = pipeline::run(&files, &options)?;
            let report = report::build(&collected, &options, &report_options);
            write_outputs(&args.output, &report)?;
            report
        }
    };

    let all = &report.groups[0].summary;
    eprintln!(
        "lrqc: {} reads, {} bases in {} group(s); wrote {} in {:.1}s",
        all.reads,
        all.bases,
        report.groups.len().saturating_sub(1).max(1),
        args.output.join("report.html").display(),
        started.elapsed().as_secs_f64()
    );
    Ok(exit_code(&report, args.fail_on))
}

/// Where `--per-read` writes. Watch mode appends each pass to a table started afresh.
fn per_read_output(args: &RunArgs) -> Result<Option<per_read::Output>> {
    if !args.per_read {
        return Ok(None);
    }
    fs::create_dir_all(&args.output)
        .with_context(|| format!("creating {}", args.output.display()))?;
    let path = per_read::path_in(&args.output);
    let mode = if args.watch.is_some() {
        if path.exists() {
            fs::remove_file(&path).with_context(|| format!("removing {}", path.display()))?;
        }
        per_read::Mode::Append
    } else {
        per_read::Mode::Create
    };
    Ok(Some(per_read::Output { path, mode }))
}

/// Exit status for `--fail-on`, naming the worst check on standard error.
fn exit_code(report: &Report, fail_on: Option<FailOn>) -> ExitCode {
    let Some(level) = fail_on else {
        return ExitCode::SUCCESS;
    };
    let threshold = match level {
        FailOn::Warn => Status::Warn,
        FailOn::Fail => Status::Fail,
    };
    match report::checks::worst(&report.groups) {
        Some((status, what)) if status >= threshold => {
            eprintln!("lrqc: {}: {what}", status.label());
            ExitCode::from(QC_FAILED)
        }
        _ => ExitCode::SUCCESS,
    }
}

/// Reads files as they appear and rewrites the report after each pass, until MinKNOW's
/// final summary appears and every file is read.
fn watch(
    args: &RunArgs,
    options: &mut Options,
    report_options: &ReportOptions,
    interval: u64,
) -> Result<Report> {
    if args.inputs.iter().any(|p| p.as_os_str() == "-") {
        bail!("--watch reads files as they appear; it cannot watch standard input");
    }
    let interval = Duration::from_secs(interval.max(1));
    let mut done: HashSet<PathBuf> = HashSet::new();
    let mut read_files: Vec<input::InputFile> = Vec::new();
    // PacBio barcodes named so far, and whether files were read before the latest of them
    // was added to the search (they are read again once the run finishes).
    let mut barcode_names = run_barcode_names([], &report_options.samples);
    let mut stale = false;
    let mut collected = Collected {
        groups: BTreeMap::new(),
        files: Vec::new(),
    };
    let mut unreadable: Vec<PathBuf> = Vec::new();
    loop {
        let mut files = input::discover_lenient(&args.inputs)?;
        // MinKNOW writes the summary as the run goes; the read files are what is watched.
        if let Some(named) = files
            .iter()
            .find(|f| f.format == Some(Format::Summary) && args.inputs.contains(&f.path))
        {
            bail!(
                "--watch reads the run's read files; {} is a sequencing summary, complete only at the end of the run",
                named.path.display()
            );
        }
        files.retain(|f| f.format != Some(Format::Summary));
        let finished = input::run_finished(&args.inputs);
        // Files read earlier that the run is now read without, as FASTQ once MinKNOW's BAM
        // of the same reads appears: everything is read again without them.
        let current: HashSet<&Path> = files.iter().map(|f| f.path.as_path()).collect();
        let before = read_files.len();
        read_files.retain(|f| current.contains(f.path.as_path()) || !f.path.exists());
        let superseded = read_files.len() != before;
        // A file still being written changes; read it once it has been still for a pass,
        // or at once when the run has finished.
        let ready: Vec<_> = files
            .into_iter()
            .filter(|f| !done.contains(&f.path))
            .filter(|f| finished || settled(&f.path, interval))
            .collect();
        if !ready.is_empty() {
            // New files may name PacBio barcodes to search for.
            if args.barcode_names.is_none() {
                let before = barcode_names.len();
                barcode_names.extend(run_barcode_names(&ready, &HashMap::new()));
                if barcode_names.len() != before {
                    options.barcodes = built_in_barcodes(&barcode_names);
                    stale |= !read_files.is_empty();
                }
            }
            eprintln!("lrqc: reading {} new file(s)", ready.len());
            let (new, read, failed) = read_pass(ready, options)?;
            collected.merge(new);
            done.extend(read.iter().map(|f| f.path.clone()));
            read_files.extend(read);
            // While the run goes on, a file that cannot be read yet (still being copied,
            // say) is read again on the next pass; the rest of the pass is counted.
            for (file, e) in &failed {
                if finished {
                    eprintln!("lrqc: {e:#}");
                } else {
                    eprintln!(
                        "lrqc: {e:#}; reading {} again on the next pass",
                        file.path.display()
                    );
                }
            }
            if finished {
                unreadable = failed.into_iter().map(|(f, _)| f.path).collect();
            }
        }
        let reread = if superseded {
            Some("files the run is now read without (FASTQ of reads also in BAM)")
        } else if finished && stale {
            Some("barcodes that appeared during the run")
        } else {
            None
        };
        if let Some(why) = reread {
            eprintln!(
                "lrqc: reading all {} file(s) again for {why}",
                read_files.len()
            );
            // The per-read table is written afresh, then appended to again.
            if let Some(output) = options.per_read.as_mut() {
                output.mode = per_read::Mode::Create;
            }
            collected = pipeline::run(&read_files, &*options)?;
            if let Some(output) = options.per_read.as_mut() {
                output.mode = per_read::Mode::Append;
            }
            stale = false;
        }
        let mut report = report::build(&collected, &*options, report_options);
        if finished {
            report.refresh_secs = None;
        }
        write_outputs(&args.output, &report)?;
        if finished {
            eprintln!("lrqc: the run has finished");
            if !unreadable.is_empty() {
                bail!(
                    "{} file(s) could not be read and are not in the report: {}",
                    unreadable.len(),
                    unreadable
                        .iter()
                        .map(|p| p.display().to_string())
                        .collect::<Vec<_>>()
                        .join(", ")
                );
            }
            return Ok(report);
        }
        thread::sleep(interval);
    }
}

/// One watch pass: what was read, the files read, and the files that failed with why.
type Pass = (
    Collected,
    Vec<input::InputFile>,
    Vec<(input::InputFile, anyhow::Error)>,
);

/// Reads one watch pass. When the batch fails, its files are read one at a time, so one
/// unreadable file does not hold back the rest; the per-read table loses the rows of each
/// failed attempt. Returns what was read, the files read, and the files that failed.
fn read_pass(ready: Vec<input::InputFile>, options: &Options) -> Result<Pass> {
    let table_len = || {
        options
            .per_read
            .as_ref()
            .map(|o| fs::metadata(&o.path).map_or(0, |m| m.len()))
    };
    let restore = |len: Option<u64>| -> Result<()> {
        if let (Some(output), Some(len)) = (options.per_read.as_ref(), len) {
            fs::OpenOptions::new()
                .write(true)
                .open(&output.path)
                .and_then(|f| f.set_len(len))
                .with_context(|| format!("restoring {}", output.path.display()))?;
        }
        Ok(())
    };
    let len = table_len();
    if let Ok(collected) = pipeline::run(&ready, options) {
        return Ok((collected, ready, Vec::new()));
    }
    restore(len)?;
    let mut collected = Collected {
        groups: BTreeMap::new(),
        files: Vec::new(),
    };
    let (mut read, mut failed) = (Vec::new(), Vec::new());
    for file in ready {
        let len = table_len();
        match pipeline::run(std::slice::from_ref(&file), options) {
            Ok(new) => {
                collected.merge(new);
                read.push(file);
            }
            Err(e) => {
                restore(len)?;
                failed.push((file, e));
            }
        }
    }
    Ok((collected, read, failed))
}

fn settled(path: &Path, interval: Duration) -> bool {
    fs::metadata(path)
        .and_then(|m| m.modified())
        .ok()
        .and_then(|modified| SystemTime::now().duration_since(modified).ok())
        .is_some_and(|age| age >= interval)
}

/// Writes every output, each through a temporary file so readers never see half of one.
fn write_outputs(dir: &Path, report: &Report) -> Result<()> {
    fs::create_dir_all(dir).with_context(|| format!("creating {}", dir.display()))?;
    let outputs = [
        ("report.html", report::html::render(report)),
        ("summary.tsv", report::tsv::render(report)),
        (
            "qc.json",
            serde_json::to_string(report).context("serializing the report")?,
        ),
        (
            "lrqc_general_stats_mqc.tsv",
            report::multiqc::general_stats(report),
        ),
        ("lrqc_checks_mqc.tsv", report::multiqc::checks(report)),
    ];
    for (name, contents) in outputs {
        let path = dir.join(name);
        let partial = dir.join(format!(".{name}.partial"));
        fs::write(&partial, contents).with_context(|| format!("writing {}", partial.display()))?;
        fs::rename(&partial, &path).with_context(|| format!("writing {}", path.display()))?;
    }
    Ok(())
}

fn compare(args: &CompareArgs) -> Result<ExitCode> {
    let mut reports = Vec::new();
    for path in &args.reports {
        let file = if path.is_dir() {
            path.join("qc.json")
        } else {
            path.clone()
        };
        let text =
            fs::read_to_string(&file).with_context(|| format!("reading {}", file.display()))?;
        let report = read_report(&text)
            .with_context(|| format!("{} is not an lrQC qc.json", file.display()))?;
        reports.push(report);
    }
    fs::create_dir_all(&args.output)
        .with_context(|| format!("creating {}", args.output.display()))?;
    let path = args.output.join("compare.html");
    fs::write(&path, report::html::render_comparison(&reports))
        .with_context(|| format!("writing {}", path.display()))?;
    eprintln!(
        "lrqc: compared {} reports; wrote {}",
        reports.len(),
        path.display()
    );
    Ok(ExitCode::SUCCESS)
}

/// Reads a qc.json, dropping modules this version no longer has (such as `pass_rate`), so
/// reports from earlier versions can still be compared.
fn read_report(text: &str) -> Result<Report> {
    let mut value: serde_json::Value = serde_json::from_str(text)?;
    // Reports from earlier versions compare where the format allows; a later version's may
    // hold what this one cannot read.
    if let Some(version) = value.get("version").and_then(|v| v.as_str()) {
        let parse = |v: &str| -> Vec<u64> { v.split('.').filter_map(|p| p.parse().ok()).collect() };
        let current = env!("CARGO_PKG_VERSION");
        if parse(version) > parse(current) {
            bail!("it was written by lrQC {version}, newer than this lrQC {current}; compare it with that version");
        } else if version != current {
            eprintln!("lrqc: comparing a report from lrqc {version} with lrqc {current}; checks added since are missing from it");
        }
    }
    let known: HashSet<&str> = ModuleId::ALL.iter().map(|id| id.key()).collect();
    if let Some(groups) = value.get_mut("groups").and_then(|g| g.as_array_mut()) {
        for group in groups {
            if let Some(modules) = group.get_mut("modules").and_then(|m| m.as_array_mut()) {
                modules.retain(|m| {
                    m.get("id")
                        .and_then(|id| id.as_str())
                        .is_some_and(|id| known.contains(id))
                });
            }
        }
    }
    Ok(serde_json::from_value(value)?)
}

fn default_title(inputs: &[PathBuf]) -> String {
    let Some(first) = inputs.first() else {
        return "lrqc".to_owned();
    };
    let name = if first.as_os_str() == "-" {
        "standard input".to_owned()
    } else {
        first.file_name().map_or_else(
            || first.display().to_string(),
            |name| name.to_string_lossy().into_owned(),
        )
    };
    match inputs.len() {
        1 => name,
        n => format!("{name} and {} more", n - 1),
    }
}
