//! Finding input files, recognizing their formats, and what their paths say about the reads.

use std::collections::HashSet;
use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, Ordering};

use anyhow::{bail, Context, Result};
use serde::{Deserialize, Serialize};

use crate::barcode;
use crate::read::Filter;

/// Formats lrqc reads.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize)]
pub enum Format {
    Bam,
    Fastq,
    Fasta,
    /// MinKNOW or dorado `sequencing_summary*.txt`.
    Summary,
}

/// One input file and what its path says about its reads.
#[derive(Debug, Clone, Serialize)]
pub struct InputFile {
    pub path: PathBuf,
    /// `None` for standard input, whose format is read from its first bytes.
    pub format: Option<Format>,
    /// Barcode named by the file or directory name.
    pub path_barcode: Option<String>,
    #[serde(skip)]
    pub path_filter: Option<Filter>,
}

/// Expands files and directories into the read files they hold, sorted for a stable order.
///
/// # Errors
///
/// Fails when a path does not exist, a named file has an unrecognized format, a directory
/// holds no reads, or summaries are mixed with sequence files (their reads would count twice).
pub fn discover(paths: &[PathBuf]) -> Result<Vec<InputFile>> {
    find(paths, Empty::Refuse)
}

/// Like [`discover`], but directories may hold no reads yet (watch mode).
///
/// # Errors
///
/// As [`discover`], except for empty directories.
pub fn discover_lenient(paths: &[PathBuf]) -> Result<Vec<InputFile>> {
    find(paths, Empty::Allow)
}

/// Whether an input directory may hold no reads.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Empty {
    Refuse,
    Allow,
}

fn find(paths: &[PathBuf], empty: Empty) -> Result<Vec<InputFile>> {
    let mut files = Vec::new();
    // Summaries found by walking a directory, as opposed to named on the command line.
    let mut walked_summaries = 0;
    for path in paths {
        if path.as_os_str() == "-" {
            files.push(InputFile {
                path: path.clone(),
                format: None,
                path_barcode: None,
                path_filter: None,
            });
            continue;
        }
        let meta = fs::metadata(path).with_context(|| format!("reading {}", path.display()))?;
        if meta.is_dir() {
            let mut found = Vec::new();
            let mut visited = HashSet::new();
            walk(path, path, &mut found, &mut visited)?;
            if found.is_empty() && empty == Empty::Refuse {
                bail!(
                    "{} holds no BAM, FASTQ, FASTA or sequencing summary files",
                    path.display()
                );
            }
            let reads = found
                .iter()
                .filter(|f| f.format != Some(Format::Summary))
                .count();
            if reads > 0 {
                // A MinKNOW run folder holds its summary beside fastq_pass/bam_pass: the
                // summary describes the same reads, so the read files are used.
                let before = found.len();
                found.retain(|f| f.format != Some(Format::Summary));
                walked_summaries += before - found.len();
            }
            files.extend(found);
        } else {
            let format = detect_format(path).with_context(|| {
                format!(
                    "{}: not a BAM, FASTQ, FASTA or sequencing summary (by name: .bam, .fastq, .fq, .fasta, .fa, .fna, optionally .gz; sequencing_summary*.txt)",
                    path.display()
                )
            })?;
            files.push(InputFile::new(path.clone(), format, path));
        }
    }
    if walked_summaries > 0 && !SUMMARY_NOTED.swap(true, Ordering::Relaxed) {
        eprintln!(
            "lrqc: using the read files, not the {walked_summaries} sequencing summar{} beside them (they describe the same reads)",
            if walked_summaries == 1 { "y" } else { "ies" }
        );
    }
    // Over every input, so `lrqc run/bam_pass run/fastq_pass` counts each read once too.
    prefer_bam_over_fastq(&mut files);
    files.sort_by(|a, b| a.path.cmp(&b.path));
    // The same file reached twice (named twice, through a symlink, or as `a` and `./a`)
    // is read once.
    let mut seen = HashSet::new();
    files.retain(|f| {
        f.format.is_none()
            || seen.insert(fs::canonicalize(&f.path).unwrap_or_else(|_| f.path.clone()))
    });

    let summaries = files
        .iter()
        .filter(|f| f.format == Some(Format::Summary))
        .count();
    if summaries > 0 && summaries < files.len() {
        bail!(
            "found {summaries} sequencing summar{} next to {} read file(s); a summary describes the same reads, so give one or the other",
            if summaries == 1 { "y" } else { "ies" },
            files.len() - summaries
        );
    }
    Ok(files)
}

/// Whether the note about skipped summaries has been printed (watch mode finds them each pass).
static SUMMARY_NOTED: AtomicBool = AtomicBool::new(false);

/// MinKNOW can write the same reads as both FASTQ and BAM (`fastq_pass/` beside
/// `bam_pass/`). Where one run folder holds both, only the BAM files are read: they carry
/// more tags, and reading both would count every read twice.
fn prefer_bam_over_fastq(files: &mut Vec<InputFile>) {
    let bam_runs: HashSet<PathBuf> = files
        .iter()
        .filter_map(|f| minknow_output(&f.path))
        .filter_map(|(run, is_bam)| is_bam.then_some(run))
        .collect();
    if bam_runs.is_empty() {
        return;
    }
    let before = files.len();
    files.retain(
        |f| !matches!(minknow_output(&f.path), Some((run, false)) if bam_runs.contains(&run)),
    );
    let dropped = before - files.len();
    // Said once: --watch discovers the files again on every pass.
    static SAID: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(false);
    if dropped > 0 && !SAID.swap(true, std::sync::atomic::Ordering::Relaxed) {
        eprintln!(
            "lrqc: the run has both BAM and FASTQ output; reading the BAM files and skipping {dropped} FASTQ file(s) of the same reads"
        );
    }
}

/// The run folder holding a file in MinKNOW's `bam_*` or `fastq_*` output folders, and
/// whether that folder holds BAM.
fn minknow_output(path: &Path) -> Option<(PathBuf, bool)> {
    path.ancestors().skip(1).find_map(|dir| {
        let name = dir.file_name()?.to_str()?.to_ascii_lowercase();
        let (kind, filter) = name.split_once('_')?;
        if !matches!(filter, "pass" | "fail") {
            return None;
        }
        let is_bam = match kind {
            "bam" => true,
            "fastq" => false,
            _ => return None,
        };
        Some((dir.parent().unwrap_or(Path::new("")).to_path_buf(), is_bam))
    })
}

/// Whether the run has finished: MinKNOW writes `final_summary_*.txt` last, and PacBio's
/// instruments `<movie>.transferdone` once a movie's files are complete.
pub fn run_finished(paths: &[PathBuf]) -> bool {
    fn holds_summary(dir: &Path) -> bool {
        fs::read_dir(dir).is_ok_and(|entries| {
            entries.flatten().any(|entry| {
                let name = entry.file_name();
                let name = name.to_string_lossy();
                (name.starts_with("final_summary") && name.ends_with(".txt"))
                    || name.ends_with(".transferdone")
            })
        })
    }
    paths.iter().any(|p| {
        // A relative path has no ancestors above the working directory to search.
        let p = &fs::canonicalize(p).unwrap_or_else(|_| p.clone());
        // MinKNOW puts the summary at the top of the run folder: the folder watched, one of
        // its subfolders (watching a folder of runs), or one or two levels above
        // (watching fastq_pass or fastq_pass/barcodeNN). Sibling runs' folders are not
        // searched: their summaries say nothing about this run.
        let below = fs::read_dir(p).is_ok_and(|entries| {
            entries
                .flatten()
                .any(|e| e.path().is_dir() && holds_summary(&e.path()))
        });
        below || p.ancestors().take(3).any(holds_summary)
    })
}

/// Run named by a MinKNOW file name, `<flow cell>_<pass|fail>_..._<run>_<n>.fastq.gz`.
pub fn run_from_path(path: &Path) -> Option<String> {
    let name = path.file_name()?.to_str()?;
    let parts: Vec<&str> = name.split('_').collect();
    (parts.len() >= 4 && matches!(parts[1], "pass" | "fail") && !parts[0].is_empty())
        .then(|| parts[0].to_owned())
}

fn walk(
    root: &Path,
    dir: &Path,
    files: &mut Vec<InputFile>,
    visited: &mut HashSet<PathBuf>,
) -> Result<()> {
    // A directory reached again through a symlink (or a symlink loop) is read once.
    if !visited.insert(fs::canonicalize(dir).unwrap_or_else(|_| dir.to_path_buf())) {
        return Ok(());
    }
    let mut entries: Vec<PathBuf> = fs::read_dir(dir)
        .with_context(|| format!("listing {}", dir.display()))?
        .map(|entry| entry.map(|e| e.path()))
        .collect::<Result<_, _>>()?;
    entries.sort();
    for path in entries {
        // Hidden files include macOS's `._name` AppleDouble files, which copy every read
        // file's name onto shared drives but hold no reads.
        if path
            .file_name()
            .and_then(|n| n.to_str())
            .is_some_and(|n| n.starts_with('.'))
        {
            continue;
        }
        if path.is_dir() {
            walk(root, &path, files, visited)?;
        } else if let Some(format) = detect_format(&path) {
            files.push(InputFile::new(path, format, root));
        }
    }
    Ok(())
}

impl InputFile {
    /// `root` is the file or directory named on the command line that led to `path`.
    fn new(path: PathBuf, format: Format, root: &Path) -> Self {
        let (inner, outer) = path_parts(&path, root);
        Self {
            path_barcode: barcode::from_path_parts(&inner, &outer),
            path_filter: filter_from_path_parts(&inner, &outer),
            path,
            format: Some(format),
        }
    }
}

/// A path's components innermost first, split into those the user named or that lie
/// below what they named (`inner`: the root's own name and everything under it), and the
/// directories above it (`outer`). What a path says about its reads must not depend on
/// whether it was typed absolute or relative, so outer components count only when they
/// are exactly a barcode or MinKNOW folder name.
fn path_parts<'a>(path: &'a Path, root: &Path) -> (Vec<&'a str>, Vec<&'a str>) {
    let all: Vec<&str> = path.iter().filter_map(|c| c.to_str()).collect();
    let split = root.iter().count().max(1).min(all.len()).saturating_sub(1);
    let inner = all[split..].iter().rev().copied().collect();
    let outer = all[..split].iter().rev().copied().collect();
    (inner, outer)
}

/// The format a file name implies.
fn detect_format(path: &Path) -> Option<Format> {
    let name = path.file_name()?.to_str()?.to_ascii_lowercase();
    let name = name.strip_suffix(".gz").unwrap_or(&name);
    if name.ends_with(".bam") {
        Some(Format::Bam)
    } else if [".fastq", ".fq"].iter().any(|ext| name.ends_with(ext)) {
        Some(Format::Fastq)
    } else if [".fasta", ".fa", ".fna"]
        .iter()
        .any(|ext| name.ends_with(ext))
    {
        Some(Format::Fasta)
    } else if name.starts_with("sequencing_summary") && name.ends_with(".txt") {
        Some(Format::Summary)
    } else {
        None
    }
}

/// Pass or fail from MinKNOW's `*_pass`/`*_fail` folders or PacBio's `fail_reads` files.
///
/// Components the user named or that lie below it are matched by word (`HG002.pass.fastq`,
/// `m84011.hifi_reads.bam`); directories above it only when named exactly as MinKNOW
/// names them, so a `qc_fail_review/` folder higher up does not fail every read.
fn filter_from_path_parts(inner: &[&str], outer: &[&str]) -> Option<Filter> {
    let by_word = |component: &&str| {
        let component = component.to_ascii_lowercase();
        let words: Vec<&str> = component
            .split(|c: char| !c.is_ascii_alphanumeric())
            .collect();
        if words.contains(&"fail") {
            Some(Filter::Fail)
        } else if words.contains(&"pass") || words.windows(2).any(|w| w == ["hifi", "reads"]) {
            Some(Filter::Pass)
        } else {
            None
        }
    };
    let minknow_folder = |component: &&str| {
        let component = component.to_ascii_lowercase();
        let (kind, filter) = component.split_once('_')?;
        if !matches!(kind, "fastq" | "bam" | "fasta" | "pod5") {
            return None;
        }
        match filter {
            "pass" => Some(Filter::Pass),
            "fail" => Some(Filter::Fail),
            _ => None,
        }
    };
    inner
        .iter()
        .find_map(by_word)
        .or_else(|| outer.iter().find_map(minknow_folder))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn filter_from_path(path: &Path) -> Option<Filter> {
        let (inner, outer) = path_parts(path, path);
        filter_from_path_parts(&inner, &outer)
    }

    #[test]
    fn reads_paths_the_same_however_they_are_typed() {
        let read = |path: &str, root: &str| {
            let f = InputFile::new(PathBuf::from(path), Format::Fastq, Path::new(root));
            (f.path_barcode, f.path_filter)
        };
        // Folders above what was named count only when exactly a barcode or MinKNOW folder.
        let file = "PAY1_x_abc_0.fastq";
        assert_eq!(
            read(
                &format!("/d/proj_bc01/run1/fastq_pass/{file}"),
                "/d/proj_bc01/run1/fastq_pass"
            ),
            read(&format!("fastq_pass/{file}"), "fastq_pass"),
        );
        assert_eq!(
            read(
                &format!("/d/qc_fail_review/run/{file}"),
                "/d/qc_fail_review/run"
            )
            .1,
            None
        );
        assert_eq!(
            read(
                &format!("/d/unclassified_samples/run/{file}"),
                "/d/unclassified_samples/run"
            )
            .0,
            None
        );
        // Exact barcode and MinKNOW folder names above what was named still count.
        assert_eq!(
            read(
                "/r/fastq_fail/barcode07/reads.fastq",
                "/r/fastq_fail/barcode07/reads.fastq"
            ),
            (Some("barcode07".to_owned()), Some(Filter::Fail))
        );
        // Below what was named, names are searched as before.
        assert_eq!(
            read("/r/run/pass_reads/x.bc2001--bc2001.fastq", "/r/run")
                .0
                .as_deref(),
            Some("bc2001--bc2001")
        );
    }

    #[test]
    fn recognizes_formats() {
        assert_eq!(detect_format(Path::new("a.bam")), Some(Format::Bam));
        assert_eq!(detect_format(Path::new("a.fastq.gz")), Some(Format::Fastq));
        assert_eq!(detect_format(Path::new("a.FQ")), Some(Format::Fastq));
        assert_eq!(detect_format(Path::new("a.fa.gz")), Some(Format::Fasta));
        assert_eq!(
            detect_format(Path::new("sequencing_summary_PAY_abc.txt")),
            Some(Format::Summary)
        );
        assert_eq!(detect_format(Path::new("a.bam.bai")), None);
    }

    #[test]
    fn finds_only_this_runs_final_summary() {
        let root = std::env::temp_dir().join(format!("lrqc-finished-{}", std::process::id()));
        let run = root.join("run_b");
        let pass = run.join("fastq_pass").join("barcode01");
        fs::create_dir_all(&pass).expect("temp dirs");
        fs::create_dir_all(root.join("run_a")).expect("temp dirs");
        fs::write(root.join("run_a").join("final_summary_A.txt"), "").expect("temp file");
        let inputs = [run.join("fastq_pass")];
        assert!(
            !run_finished(&inputs),
            "another run's summary does not count"
        );
        fs::write(run.join("final_summary_B.txt"), "").expect("temp file");
        assert!(run_finished(&inputs));
        assert!(run_finished(&[pass]));
        assert!(
            run_finished(std::slice::from_ref(&root)),
            "watching a folder of runs"
        );
        // PacBio instruments mark a finished movie with `<movie>.transferdone`.
        let movie = root.join("r84011_20240101_000000").join("1_A01");
        fs::create_dir_all(&movie).expect("temp dirs");
        assert!(!run_finished(std::slice::from_ref(&movie)));
        fs::write(movie.join("m84011_240101_000000_s1.transferdone"), "").expect("temp file");
        assert!(run_finished(std::slice::from_ref(&movie)));
        fs::remove_dir_all(&root).expect("cleanup");
    }

    #[test]
    fn reads_runs_from_minknow_names() {
        assert_eq!(
            run_from_path(Path::new("PAY12345_pass_barcode01_abc123_0.fastq.gz")).as_deref(),
            Some("PAY12345")
        );
        assert_eq!(run_from_path(Path::new("reads.fastq.gz")), None);
    }

    #[test]
    fn reads_filter_from_path() {
        assert_eq!(
            filter_from_path(Path::new("/r/fastq_pass/barcode01/x.fastq.gz")),
            Some(Filter::Pass)
        );
        assert_eq!(
            filter_from_path(Path::new("/r/bam_fail/x.bam")),
            Some(Filter::Fail)
        );
        assert_eq!(
            filter_from_path(Path::new("/r/m84011.fail_reads.bam")),
            Some(Filter::Fail)
        );
        assert_eq!(
            filter_from_path(Path::new("/r/m84011.hifi_reads.bam")),
            Some(Filter::Pass)
        );
        assert_eq!(filter_from_path(Path::new("/r/reads.bam")), None);
    }
}
