//! End-to-end runs of the binary on small real-read fixtures in `tests/data`.
//!
//! Expected counts come from samtools and awk over the same files, not from lrqc.

use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

/// One row of `summary.tsv`, by column name.
type Row = HashMap<String, String>;

/// Runs lrqc on `inputs` and returns its summary rows by group, plus the output directory.
fn run(name: &str, inputs: &[&str], extra: &[&str]) -> (HashMap<String, Row>, PathBuf) {
    let data = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data");
    let out = std::env::temp_dir().join(format!("lrqc-test-{name}-{}", std::process::id()));
    let _ = fs::remove_dir_all(&out);
    let status = Command::new(env!("CARGO_BIN_EXE_lrqc"))
        .args(inputs.iter().map(|i| data.join(i)))
        .args(extra)
        .arg("-o")
        .arg(&out)
        .arg("-t")
        .arg("2")
        .status()
        .expect("lrqc runs");
    assert!(status.success(), "lrqc failed on {inputs:?}");
    let tsv = fs::read_to_string(out.join("summary.tsv")).expect("summary.tsv written");
    let mut lines = tsv.lines();
    let header: Vec<&str> = lines.next().expect("header").split('\t').collect();
    let rows = lines
        .map(|line| {
            let row: Row = header
                .iter()
                .zip(line.split('\t'))
                .map(|(k, v)| ((*k).to_owned(), v.to_owned()))
                .collect();
            (row["group"].clone(), row)
        })
        .collect();
    (rows, out)
}

fn assert_yield(rows: &HashMap<String, Row>, group: &str, reads: &str, bases: &str) {
    let row = rows
        .get(group)
        .unwrap_or_else(|| panic!("no group {group} in {:?}", rows.keys()));
    assert_eq!(
        (row["reads"].as_str(), row["bases"].as_str()),
        (reads, bases),
        "group {group}"
    );
}

#[test]
fn splits_dorado_bam_by_bc_tag() {
    let (rows, out) = run("dorado", &["ont_dorado.bam"], &[]);
    assert_yield(&rows, "All reads", "60", "1053839");
    assert_yield(&rows, "barcode01", "20", "508656");
    assert_yield(&rows, "barcode02", "16", "312544");
    assert_yield(&rows, "barcode03", "12", "21181");
    assert_yield(&rows, "barcode04", "4", "48543");
    assert_yield(&rows, "unclassified", "8", "162915");
    assert_eq!(rows["All reads"]["platform"], "Oxford Nanopore");
    // barcode02's qualities were lowered by 9 when the fixture was made.
    let q = |g: &str| rows[g]["median_read_q"].parse::<f64>().expect("a number");
    assert!(q("barcode02") < q("barcode01") - 5.0);

    let html = fs::read_to_string(out.join("report.html")).expect("report written");
    assert!(html.starts_with("<!doctype html>"));
    assert!(!html.contains("<script"), "the report needs no script");
    assert!(html.contains(r#"id="g1-barcode01""#));
    let json: serde_json::Value =
        serde_json::from_str(&fs::read_to_string(out.join("qc.json")).expect("qc.json written"))
            .expect("valid JSON");
    assert_eq!(json["groups"].as_array().map(Vec::len), Some(6));
}

#[test]
fn names_lima_files_by_path() {
    let (rows, _) = run("lima", &["lima"], &[]);
    assert_yield(&rows, "All reads", "50", "915312");
    assert_yield(&rows, "bc2001--bc2001", "25", "481355");
    assert_yield(&rows, "bc2002--bc2002", "25", "433957");
    assert_eq!(rows["All reads"]["platform"], "PacBio");
    assert_eq!(
        rows["bc2001--bc2001"]["pass_pct"], "100.00",
        "hifi_reads files pass"
    );
}

#[test]
fn reads_minknow_folders() {
    let (rows, _) = run("minknow", &["minknow"], &[]);
    assert_yield(&rows, "barcode01", "20", "422247");
    assert_yield(&rows, "barcode02", "10", "242002");
    assert_eq!(
        rows["barcode01"]["pass_pct"], "100.00",
        "fastq_pass reads pass"
    );
    assert_eq!(
        rows["barcode02"]["pass_pct"], "0.00",
        "fastq_fail reads fail"
    );
}

#[test]
fn reads_sequencing_summaries() {
    let (rows, _) = run("summary", &["sequencing_summary_test.txt"], &[]);
    assert_yield(&rows, "All reads", "300", "5721364");
    assert_yield(&rows, "barcode01", "106", "2016897");
    assert_yield(&rows, "unclassified", "37", "794963");
    assert_eq!(rows["All reads"]["pass_reads"], "298");
    assert_eq!(rows["All reads"]["yield_over_time"], "info");
    assert_eq!(
        rows["All reads"]["base_quality"], "n/a",
        "summaries have no bases"
    );
}

#[test]
fn refuses_summaries_mixed_with_reads() {
    let data = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data");
    let output = Command::new(env!("CARGO_BIN_EXE_lrqc"))
        .arg(data.join("sequencing_summary_test.txt"))
        .arg(data.join("ont_dorado.bam"))
        .arg("-o")
        .arg(std::env::temp_dir().join(format!("lrqc-test-mixed-{}", std::process::id())))
        .output()
        .expect("lrqc runs");
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("sequencing summar"));
}

fn data(path: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(path)
}

fn json(out: &Path) -> serde_json::Value {
    serde_json::from_str(&fs::read_to_string(out.join("qc.json")).expect("qc.json written"))
        .expect("valid JSON")
}

#[test]
fn reads_standard_input() {
    let out = std::env::temp_dir().join(format!("lrqc-test-stdin-{}", std::process::id()));
    let _ = fs::remove_dir_all(&out);
    let bam = fs::read(data("ont_dorado.bam")).expect("fixture");
    let mut child = Command::new(env!("CARGO_BIN_EXE_lrqc"))
        .args(["-", "-t", "2", "-o"])
        .arg(&out)
        .stdin(std::process::Stdio::piped())
        .spawn()
        .expect("lrqc runs");
    std::io::Write::write_all(&mut child.stdin.take().expect("stdin"), &bam).expect("pipe");
    assert!(child.wait().expect("lrqc exits").success());
    assert_eq!(json(&out)["groups"][0]["summary"]["bases"], 1_053_839);
}

#[test]
fn screens_for_contaminants() {
    let (rows, out) = run("screen", &["synthetic/ecoli_spike.fastq"], &[]);
    assert_eq!(
        rows["All reads"]["contaminant_pct"], "20.00",
        "10 of 50 reads are E. coli"
    );
    assert_eq!(rows["All reads"]["contaminants"], "fail");
    let hits = &json(&out)["groups"][0]["data"]["contaminants"];
    let ecoli = hits
        .as_array()
        .expect("hits")
        .iter()
        .find(|h| h["reference"]["key"] == "ecoli")
        .expect("E. coli row");
    assert_eq!(ecoli["reads"], 10);
}

#[test]
fn reads_cdna_structure_and_poly_a() {
    let (_, out) = run("cdna", &["synthetic/cdna.fastq"], &[]);
    let data = &json(&out)["groups"][0]["data"];
    // 30 forward and 10 reverse full-length reads of 60; 10 with the SSP only.
    assert_eq!(data["cdna"]["forward_pct"], 50.0);
    assert!((data["cdna"]["reverse_pct"].as_f64().expect("number") - 100.0 / 6.0).abs() < 1e-9);
    assert!((data["cdna"]["ssp_only_pct"].as_f64().expect("number") - 100.0 / 6.0).abs() < 1e-9);
    assert_eq!(
        data["poly_a"]
            .as_array()
            .map(|b| b.iter().filter_map(|c| c.as_u64()).sum::<u64>()),
        Some(60)
    );
}

#[test]
fn counts_kinnex_adapters() {
    let adapters = data("synthetic/kinnex_adapters.fasta");
    let (_, out) = run(
        "kinnex",
        &["synthetic/kinnex.fastq"],
        &["--kinnex-adapters", adapters.to_str().expect("UTF-8 path")],
    );
    let kinnex = json(&out)["groups"][0]["data"]["kinnex"].clone();
    assert_eq!(kinnex, serde_json::json!([0, 0, 0, 4, 0, 6]));
}

#[test]
fn names_samples_and_reports_coverage() {
    let sheet = data("samples.csv");
    let (rows, out) = run(
        "samples",
        &["ont_dorado.bam"],
        &[
            "--sample-sheet",
            sheet.to_str().expect("UTF-8 path"),
            "--genome-size",
            "1m",
        ],
    );
    assert_eq!(rows["barcode01"]["sample"], "patient_A");
    // 508,656 bases over 1 Mb.
    assert_eq!(rows["barcode01"]["coverage"], "0.51");
    let mqc = fs::read_to_string(out.join("lrqc_general_stats_mqc.tsv")).expect("MultiQC file");
    assert!(mqc.contains("plot_type: 'generalstats'") && mqc.contains("patient_A"));
}

#[test]
fn fails_on_request() {
    let thresholds =
        std::env::temp_dir().join(format!("lrqc-test-thresholds-{}.toml", std::process::id()));
    fs::write(
        &thresholds,
        "contaminant_fail = 50.0\ncontaminant_warn = 30.0\n",
    )
    .expect("write");
    let status = |extra: &[&str]| {
        Command::new(env!("CARGO_BIN_EXE_lrqc"))
            .arg(data("synthetic/ecoli_spike.fastq"))
            .args(extra)
            .arg("-o")
            .arg(std::env::temp_dir().join(format!("lrqc-test-failon-{}", std::process::id())))
            .output()
            .expect("lrqc runs")
            .status
            .code()
    };
    assert_eq!(status(&["--fail-on", "fail"]), Some(2));
    // Raising the contaminant levels above 20% leaves nothing failing.
    let raised = [
        "--fail-on",
        "fail",
        "--thresholds",
        thresholds.to_str().expect("UTF-8 path"),
    ];
    assert_eq!(status(&raised), Some(0));
}

#[test]
fn splits_runs_and_compares_reports() {
    let (rows, first) = run("runs", &["ont_dorado.bam", "minknow"], &["--split-runs"]);
    assert_yield(&rows, "PAY87954 · barcode01", "20", "508656");
    assert_yield(&rows, "PAY12345 · barcode01", "20", "422247");
    let (_, second) = run("runs2", &["lima"], &[]);
    let out = std::env::temp_dir().join(format!("lrqc-test-compare-{}", std::process::id()));
    let status = Command::new(env!("CARGO_BIN_EXE_lrqc"))
        .arg("compare")
        .args([&first, &second])
        .arg("-o")
        .arg(&out)
        .status()
        .expect("lrqc runs");
    assert!(status.success());
    let html = fs::read_to_string(out.join("compare.html")).expect("comparison written");
    assert!(html.contains("PAY12345 · barcode01") && html.contains("bc2001--bc2001"));
}

#[test]
fn compares_reports_holding_retired_modules() {
    let (_, first) = run("retired", &["ont_dorado.bam"], &[]);
    // A qc.json from an earlier version, with the pass/fail module since removed.
    let mut old = json(&first);
    for group in old["groups"].as_array_mut().expect("groups") {
        group["modules"]
            .as_array_mut()
            .expect("modules")
            .push(serde_json::json!({"id": "pass_rate", "status": "Pass", "note": "100.0% pass"}));
    }
    let old_path = first.join("old.json");
    fs::write(&old_path, old.to_string()).expect("writes old.json");
    let out = std::env::temp_dir().join(format!("lrqc-test-retired-{}", std::process::id()));
    let status = Command::new(env!("CARGO_BIN_EXE_lrqc"))
        .arg("compare")
        .arg(&old_path)
        .arg(&first)
        .arg("-o")
        .arg(&out)
        .status()
        .expect("lrqc runs");
    assert!(
        status.success(),
        "compare reads a qc.json with a retired module"
    );
}

#[test]
fn writes_the_per_read_table_in_input_order() {
    let (rows, out) = run("perread", &["ont_dorado.bam"], &["--per-read"]);
    let mut text = String::new();
    std::io::Read::read_to_string(
        &mut flate2::read::MultiGzDecoder::new(
            fs::File::open(out.join("reads.tsv.gz")).expect("table written"),
        ),
        &mut text,
    )
    .expect("valid gzip");
    let mut lines = text.lines();
    let header: Vec<&str> = lines.next().expect("header").split('\t').collect();
    let column = |name: &str| header.iter().position(|c| *c == name).expect("column");
    let table: Vec<Vec<&str>> = lines.map(|l| l.split('\t').collect()).collect();
    assert_eq!(table.len(), 60);
    let bases: u64 = table
        .iter()
        .map(|r| r[column("length")].parse::<u64>().expect("length"))
        .sum();
    assert_eq!(bases.to_string(), rows["All reads"]["bases"]);
    let barcode01 = table
        .iter()
        .filter(|r| r[column("group")] == "barcode01")
        .count();
    assert_eq!(barcode01.to_string(), rows["barcode01"]["reads"]);
    assert!(table.iter().all(|r| r.len() == header.len()));
    assert!(
        table[0][column("start_time")].ends_with('Z'),
        "ISO 8601 times"
    );
}

/// A fresh temporary directory for one test.
fn temp_dir(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("lrqc-test-{name}-{}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).expect("temp dir");
    dir
}

/// Runs lrqc on paths outside `tests/data`; returns its output directory.
fn run_paths(name: &str, inputs: &[&Path], extra: &[&str]) -> PathBuf {
    let out = std::env::temp_dir().join(format!("lrqc-test-out-{name}-{}", std::process::id()));
    let _ = fs::remove_dir_all(&out);
    let output = Command::new(env!("CARGO_BIN_EXE_lrqc"))
        .args(inputs)
        .args(extra)
        .arg("-o")
        .arg(&out)
        .args(["-t", "2"])
        .output()
        .expect("lrqc runs");
    assert!(
        output.status.success(),
        "lrqc failed on {inputs:?}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    out
}

fn summary_reads(out: &Path) -> u64 {
    json(out)["groups"][0]["summary"]["reads"]
        .as_u64()
        .expect("reads")
}

#[test]
fn reads_a_whole_minknow_run_folder() {
    // A run folder holds its sequencing summary beside fastq_pass; the summary describes
    // the same reads, so the read files are used and the run is not refused.
    let run = temp_dir("runfolder");
    for (sub, file) in [
        (
            "fastq_pass/barcode01",
            "minknow/fastq_pass/barcode01/PAY12345_pass_barcode01_abc_0.fastq.gz",
        ),
        (
            "fastq_fail/barcode02",
            "minknow/fastq_fail/barcode02/PAY12345_fail_barcode02_abc_0.fastq.gz",
        ),
    ] {
        fs::create_dir_all(run.join(sub)).expect("dirs");
        fs::copy(data(file), run.join(sub).join("reads_0.fastq.gz")).expect("copy");
    }
    fs::copy(
        data("sequencing_summary_test.txt"),
        run.join("sequencing_summary_X.txt"),
    )
    .expect("copy");
    fs::write(run.join("final_summary_X.txt"), "").expect("write");
    // macOS copies leave AppleDouble files beside every read file.
    fs::write(
        run.join("fastq_pass/barcode01/._reads_0.fastq.gz"),
        b"\0\x05\x16\x07",
    )
    .expect("write");
    let out = run_paths("runfolder", &[&run], &[]);
    assert_eq!(summary_reads(&out), 30);
}

#[test]
fn reads_bam_once_when_minknow_wrote_fastq_too() {
    let run = temp_dir("bothformats");
    fs::create_dir_all(run.join("bam_pass")).expect("dirs");
    fs::create_dir_all(run.join("fastq_pass")).expect("dirs");
    fs::copy(data("ont_dorado.bam"), run.join("bam_pass/calls_0.bam")).expect("copy");
    fs::copy(
        data("synthetic/ecoli_spike.fastq"),
        run.join("fastq_pass/calls_0.fastq"),
    )
    .expect("copy");
    let out = run_paths("bothformats", &[&run], &[]);
    assert_eq!(summary_reads(&out), 60, "only the BAM's reads");
}

#[test]
#[cfg(unix)]
fn reads_each_file_once() {
    let dir = temp_dir("dedup");
    fs::create_dir_all(dir.join("a")).expect("dirs");
    fs::copy(
        data("synthetic/ecoli_spike.fastq"),
        dir.join("a/reads.fastq"),
    )
    .expect("copy");
    std::os::unix::fs::symlink("..", dir.join("a/up")).expect("symlink loop");
    let out = run_paths("dedup", &[&dir], &[]);
    assert_eq!(summary_reads(&out), 50);
    let file = dir.join("a/reads.fastq");
    let again = dir.join("a/../a/reads.fastq");
    let out = run_paths("dedup2", &[&file, &again], &[]);
    assert_eq!(summary_reads(&out), 50);
}

#[test]
fn keeps_the_name_of_a_single_barcode() {
    let (rows, _) = run(
        "single",
        &["minknow/fastq_pass"],
        &[
            "--sample-sheet",
            data("samples.csv").to_str().expect("utf-8 path"),
        ],
    );
    assert_yield(&rows, "barcode01", "20", "422247");
    assert_eq!(rows["barcode01"]["sample"], "patient_A");
    assert!(!rows.contains_key("All reads"));
}

#[test]
fn passes_hifi_fastq_without_rq() {
    // SMRT Link's FASTQ export carries no rq; binned base qualities put the estimated read
    // Q below Q20 for some HiFi reads, which must not fail them.
    let dir = temp_dir("hifi");
    let mut fastq = String::new();
    for i in 0..20 {
        // Mostly Q40 with Q3 bases, as Revio bins them: estimated read Q about Q17.
        let quals: String = (0..1000)
            .map(|j| if j % 25 == 0 { '$' } else { 'I' })
            .collect();
        let seq: String = (0..1000)
            .map(|j| ['A', 'C', 'G', 'T'][(j * 7 + i) % 4])
            .collect();
        fastq.push_str(&format!(
            "@m84011_220902_175841_s1/{i}/ccs\n{seq}\n+\n{quals}\n"
        ));
    }
    let path = dir.join("m84011.reads.fastq");
    fs::write(&path, fastq).expect("write");
    let out = run_paths("hifi", &[&path], &["--no-screen"]);
    let json = json(&out);
    let group = &json["groups"][0];
    assert_eq!(group["platform"], "PacBio");
    assert_eq!(group["summary"]["pass_reads"], 20);
    let quality = group["modules"]
        .as_array()
        .expect("modules")
        .iter()
        .find(|m| {
            m["note"]
                .as_str()
                .is_some_and(|n| n.contains("median read Q"))
        })
        .expect("read quality module");
    assert_eq!(quality["status"], "Info", "{quality}");
}

/// A FASTQ of dorado-style reads: `(run, barcode, count, start offset in minutes)`.
fn dorado_fastq(path: &Path, groups: &[(&str, &str, usize, u32)]) {
    let mut fastq = String::new();
    let mut n = 0;
    for &(run, barcode, count, start) in groups {
        for i in 0..count {
            let seq: String = (0..400)
                .map(|j| ['A', 'C', 'G', 'T'][(j * 7 + n) % 4])
                .collect();
            let minute = start + i as u32;
            fastq.push_str(&format!(
                "@r{n}\tqs:f:15\tRG:Z:{run}_dna_r10.4.1_e8.2_400bps_sup@v5.0.0\tBC:Z:{barcode}\tst:Z:2024-03-01T{:02}:{:02}:00.000+00:00\n{seq}\n+\n{}\n",
                minute / 60,
                minute % 60,
                "0".repeat(400)
            ));
            n += 1;
        }
    }
    fs::write(path, fastq).expect("write");
}

#[test]
fn grades_split_runs_within_each_run() {
    let dir = temp_dir("splitruns");
    let path = dir.join("runs.fastq");
    // Run b is a tenth of run a, evenly split; each run is 15% unclassified.
    dorado_fastq(
        &path,
        &[
            ("aaaa1111", "barcode01", 170, 0),
            ("aaaa1111", "barcode02", 170, 0),
            ("aaaa1111", "unclassified", 60, 0),
            ("bbbb2222", "barcode01", 17, 0),
            ("bbbb2222", "barcode02", 17, 0),
            ("bbbb2222", "unclassified", 6, 0),
        ],
    );
    let out = run_paths("splitruns", &[&path], &["--split-runs", "--no-screen"]);
    let json = json(&out);
    for group in json["groups"].as_array().expect("groups") {
        let name = group["name"].as_str().expect("name");
        let balance = group["modules"]
            .as_array()
            .expect("modules")
            .iter()
            .find(|m| m["id"] == "barcode_balance")
            .expect("balance module");
        if name.ends_with("unclassified") {
            assert_eq!(
                balance["status"], "Warn",
                "{name}: 15% of its run: {balance}"
            );
        } else if name.contains("barcode") {
            assert_eq!(
                balance["status"], "Pass",
                "{name} is balanced in its run: {balance}"
            );
        }
    }
}

#[test]
fn starts_every_barcode_at_the_run_start_and_names_multiqc_rows_uniquely() {
    let dir = temp_dir("timeline");
    let path = dir.join("reads.fastq");
    dorado_fastq(
        &path,
        &[
            ("aaaa1111", "barcode01", 120, 0),
            ("aaaa1111", "barcode02", 60, 60),
        ],
    );
    let sheet = dir.join("sheet.csv");
    fs::write(&sheet, "barcode,alias\nbarcode01,same\nbarcode02,same\n").expect("write");
    let out = run_paths(
        "timeline",
        &[&path],
        &[
            "--no-screen",
            "--sample-sheet",
            sheet.to_str().expect("utf-8"),
        ],
    );
    let json = json(&out);
    for group in json["groups"].as_array().expect("groups") {
        let time = group["data"]["time"].as_array().expect("time");
        assert_eq!(time.len(), 120, "one slice per minute, none extra");
        let first = time
            .iter()
            .find(|t| t["reads"].as_u64() > Some(0))
            .expect("reads");
        let expected = if group["name"] == "barcode02" {
            1.0
        } else {
            0.0
        };
        assert_eq!(first["hours"].as_f64(), Some(expected), "{}", group["name"]);
    }
    let stats = fs::read_to_string(out.join("lrqc_general_stats_mqc.tsv")).expect("written");
    let names: Vec<&str> = stats
        .lines()
        .filter(|l| !l.starts_with('#'))
        .skip(1)
        .map(|l| l.split('\t').next().expect("name"))
        .collect();
    let unique: std::collections::HashSet<&&str> = names.iter().collect();
    assert_eq!(
        unique.len(),
        names.len(),
        "MultiQC keeps one row per name: {names:?}"
    );
    for column in ["bases", "pass_pct", "gc_pct"] {
        assert!(
            stats.contains(&format!("#     {column}:")),
            "{column} needs a pconfig entry"
        );
    }
}

/// Copies `from` to `to`, creating its directory, dated a minute ago so `--watch` takes it
/// as settled at once.
fn place(from: &Path, to: &Path) {
    fs::create_dir_all(to.parent().expect("a parent")).expect("creates the directory");
    fs::copy(from, to).expect("copies the file");
    let past = std::time::SystemTime::now() - std::time::Duration::from_secs(60);
    fs::File::options()
        .write(true)
        .open(to)
        .and_then(|f| f.set_modified(past))
        .expect("dates the file");
}

#[test]
fn watch_counts_readable_files_despite_an_unreadable_one() {
    let dir = std::env::temp_dir().join(format!("lrqc-test-watch-bad-{}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    let good = data("minknow/fastq_pass/barcode01/PAY12345_pass_barcode01_abc_0.fastq.gz");
    place(&good, &dir.join("run/fastq_pass/good.fastq.gz"));
    // The first 300 bytes of a gzip file: truncated, as a crashed copy leaves it.
    let bytes = fs::read(&good).expect("reads the fixture");
    fs::write(dir.join("run/fastq_pass/bad.fastq.gz"), &bytes[..300]).expect("writes");
    fs::write(dir.join("run/final_summary_x.txt"), "").expect("writes the final summary");
    let out = dir.join("out");
    let output = Command::new(env!("CARGO_BIN_EXE_lrqc"))
        .arg(dir.join("run"))
        .args(["--watch", "1", "-t", "2", "-o"])
        .arg(&out)
        .output()
        .expect("lrqc runs");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(!output.status.success(), "the unreadable file is an error");
    assert!(stderr.contains("bad.fastq.gz"), "{stderr}");
    assert_eq!(
        json(&out)["groups"][0]["summary"]["reads"],
        20,
        "the good file is counted"
    );
}

#[test]
fn watch_reads_fastq_once_the_bam_of_the_same_run_appears() {
    let dir = std::env::temp_dir().join(format!("lrqc-test-watch-bam-{}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    let run = dir.join("run");
    place(
        &data("minknow/fastq_pass/barcode01/PAY12345_pass_barcode01_abc_0.fastq.gz"),
        &run.join("fastq_pass/barcode01/x.fastq.gz"),
    );
    let out = dir.join("out");
    let mut child = Command::new(env!("CARGO_BIN_EXE_lrqc"))
        .arg(&run)
        .args(["--watch", "1", "-t", "2", "-o"])
        .arg(&out)
        .stderr(std::process::Stdio::null())
        .spawn()
        .expect("lrqc starts");
    // The FASTQ is read on the first pass; the BAM appears after it.
    std::thread::sleep(std::time::Duration::from_millis(2500));
    place(&data("ont_dorado.bam"), &run.join("bam_pass/a.bam"));
    fs::write(run.join("final_summary_x.txt"), "").expect("writes the final summary");
    assert!(child.wait().expect("lrqc ends").success());
    // A one-off run of the finished folder reads only the BAM.
    let once = dir.join("once");
    let status = Command::new(env!("CARGO_BIN_EXE_lrqc"))
        .arg(&run)
        .args(["-t", "2", "-o"])
        .arg(&once)
        .stderr(std::process::Stdio::null())
        .status()
        .expect("lrqc runs");
    assert!(status.success());
    assert_eq!(
        json(&out)["groups"][0]["summary"]["reads"],
        json(&once)["groups"][0]["summary"]["reads"]
    );
}
