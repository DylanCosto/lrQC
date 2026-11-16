# Using lrQC

Everything the [README](../README.md) leaves out: what lrQC reads and how, every check and its
levels, thresholds, watching a run, comparing reports, every option, performance and limits.
[How lrQC works](methods.md) has the methods and their validation.

## Inputs

| Input | Recognized by |
|---|---|
| BAM, aligned or unaligned | `.bam` |
| FASTQ, FASTA | `.fastq`, `.fq`, `.fasta`, `.fa`, `.fna`, each optionally `.gz` (gzip or BGZF); FASTQ may be wrapped over several lines |
| Sequencing summary | `sequencing_summary*.txt` |
| Directory | searched recursively for the above |
| Standard input | `-`; the format is read from the first bytes, gzip or not |

CRAM and POD5 are not read; convert CRAM to BAM first (`samtools view -b`). A summary
describes the same reads as the run's FASTQ or BAM, so lrQC refuses to mix the two when both
are named; in a directory (a MinKNOW run folder), the read files are used and the summary
beside them is skipped. Where MinKNOW wrote both BAM
and FASTQ (`bam_pass/` beside `fastq_pass/`), only the BAM files are read, whether the run
folder or the two output folders are given. Each file is read once, however it is reached
(named twice, through symlinks); hidden files found in directories, such as macOS's `._*`
copies, are skipped. A BGZF file (BAM, bgzipped FASTQ) that lacks its end-of-file block is
truncated, and lrQC stops with an error rather than report part of it.

Tag values no basecaller writes (a `qs` or `rq` that is not a number, negative, infinite or
out of range) are ignored, and the read falls back to what its bases say.

### Barcodes and samples

Each read's barcode comes from the first of these that names one:

1. The record: dorado's `BC:Z` tag (`SQK-NBD114-24_barcode01` → `barcode01`; a `BC` holding
   a barcode's sequence, as the SAM specification has it, is ignored), guppy/MinKNOW's
   `barcode=` in a FASTQ header, a summary's `barcode_arrangement`, or PacBio lima's
   `bc:B:S` index pair. Give lima's barcode FASTA with `--barcode-names` to name pairs
   `bc2001--bc2001`; without it they are `index0--index0`, unless rule 3 names them.
2. The read group: dorado's `…_barcode01` read group IDs (in the BAM header, in a BAM
   record's `RG:Z:` even when the header lacks it, or in the FASTQ header's `RG:Z:` tag),
   SMRT Link's `<hash>/0--0`.
3. The path: `barcode01`, `bc2001--bc2001`, `bc2045` or `unclassified` as a file name
   part or a directory (MinKNOW's `fastq_pass/barcode01/`, lima's
   `*.bc2001--bc2001.bam`). Directories above the file or directory given count only
   when named exactly so (`barcode01/`, not `proj_bc01/`), so the result does not depend
   on whether the path was typed absolute or relative.

`unassigned`, `unbarcoded` and lima's `default--default` become `unclassified`. Reads with no
barcode information form `no_barcode`. A run whose reads all fall in one barcode is reported as
that barcode (with its sample name), without a separate "all reads"; a run without barcodes
is "All reads".

`--sample-sheet` names barcodes by sample: a MinKNOW sample sheet (`barcode`, `alias`
columns) or a SMRT Link biosample CSV (`Barcode`, `Bio Sample`, or the older `Barcode Name`,
`Bio Sample Name`). It is read as CSV, quotes and a spreadsheet's byte order mark included.

`--split-runs` reports each run and barcode apart (`PAY87954 · barcode01`). The run is the
flow cell or movie of the read group (`PU`), the run ID of dorado's read group or guppy's
`runid=`, the movie of a PacBio read name, the flow cell at the front of MinKNOW's file
names, or a summary's `run_id`.

### Pass and fail

A read passes or fails by MinKNOW's `fastq_pass`/`fastq_fail` (or `bam_*`) folders, PacBio's
`hifi_reads`/`fail_reads` files, a summary's `passes_filtering`, or the QC-fail flag. With
none of these, it passes when its read Q reaches `--min-qscore` (default Q10 for ONT, Q20 for
PacBio, the HiFi threshold). Path names count by word: `pass` or `fail` in the name of the
file, or of a directory at or below the path given (`HG002.pass.fastq`); directories above the
path given count only when named exactly as MinKNOW names them (`fastq_pass`, `bam_fail`), so
a `qc_fail_review/` folder higher up fails nothing. PacBio CCS
reads (`…/ccs`) without `rq`, as in SMRT Link's FASTQ export, pass unless `--min-qscore` is
given: CCS output holds HiFi reads, and their read Q is only an estimate (below).

### Read Q

A read's Q is, in order: the basecaller's (`qs:f` tag, summary `mean_qscore_template`),
PacBio's predicted accuracy (`rq`, as −10·log10(1 − rq)), or the mean error probability of
its bases (dorado's method: the bases after the first 60, or every base of a piece dorado split from a longer read; this matched dorado's own `qs` to within 0.01 on 18,235 reads). PacBio's `qs:i` is a query start and is
ignored as a quality. PacBio base qualities are binned (Revio uses seven levels), so a read Q
from them runs several Q off `rq` either way; when most of a group's reads have no `rq`, its
read quality is shown without a verdict.

Reads downloaded from the SRA lose their headers and names, so lrQC cannot tell the platform
from them: give `--platform pacbio` for HiFi reads, or they are graded at ONT's Q10.

## Modules and checks

Every threshold is printed beside its verdict in the report, and every one can be changed
(see [Thresholds](#thresholds)). [How lrQC works](methods.md) explains where the checks
depart from FastQC, how the contaminant screen and barcode search work, and how they were
validated.

| Module | Verdict | Warn | Fail |
|---|---|---|---|
| Basic statistics | pass | – | no reads |
| Yield after filtering | info | reads, bases and coverage at min length × min read Q | |
| Read length distribution | info | zero-length reads | – |
| Read quality distribution | graded | median read Q < Q10 ONT, Q20 PacBio (or `--min-qscore`) | < 3/4 of that |
| Read length vs quality | info | | |
| Per-position base quality | graded | a position bin's median below the read Q warn level, or its lower quartile below the fail level | a median below the fail level |
| Per-position base content | graded | A/T or G/C differ by > 10% beyond the first and last 20 bases | > 20% |
| Per-read GC content | graded | a second GC peak holding > 5% of reads of 500+ bases | > 20% |
| Per-position N content | graded | > 5% N at a position | > 20% |
| Duplication levels | graded | > 20% of sampled reads are copies beyond the first (both 50-base ends match, FastQC's measure) | > 50% |
| Overrepresented end sequences | graded | a 50-base read start or end in > 0.1% of sampled reads (sequences above the warn level or 0.1% are listed) | > 1% |
| Adapter and primer content | graded | > 5% of reads with any adapter or primer within 150 bases of an end, or > 2% with one inside | > 10% / > 5% |
| Controls and contaminants | graded | > 1% of reads from a contaminant | > 5% |
| Barcodes inside reads | graded | > 5% of reads keep their own barcode at an end, or > 1% hold another sample's | > 20% / > 5% |
| Duplicate read IDs | graded | any | > 1% of reads |
| Barcode balance | graded | a barcode below 0.25× the median barcode's reads; > 10% unclassified (with `--split-runs`, within each run) | < 0.1×; > 25% |
| cDNA and poly(A) | info | full-length fraction, strand, poly(A) lengths (`pt`) | |
| Yield over time | info | ONT: yield, read Q, pass rate and read length over the run | |
| Channel activity | info | ONT: `ch` tags or summary | |
| Translocation speed | info | ONT: read length / `du` | |
| Duplex | info | ONT: duplex classes (`dx`, or a duplex summary's `duplex_parent_template`), their yield and Q, the duplex rate | |
| Signal, muxes and split reads | info | ONT: signal shift and scale (`sm`, `sd`), samples per base (`ns`), mux yield (`mx`), reads split by dorado (`pi`) | |
| Read end reasons | info | dorado `er` tags or summaries | |
| PacBio HiFi metrics | info | passes, passes vs read Q, effective coverage, lima score, reads per ZMW, kinetics | |
| Kinnex arrays | info | segmentation adapters per read (`--kinnex-adapters`) | |
| Base modifications | info | MM/ML calls per code, % modified, % modified at CpG | |
| Alignment | graded | < 90% of primary reads mapped | < 70% |

The contaminant screen compares every read with fourteen references built into the binary:
lambda (ONT DCS), yeast ENO2 (ONT RNA calibration strand), ERCC and SIRV spike-ins, phiX,
E. coli K-12, *Mycoplasma hyorhinis*, *orale*, *fermentans* and *arginini*, *Acholeplasma
laidlawii*, the human mitochondrion, the human nuclear genome and cloning vectors (UniVec
Core). Contaminants are graded; controls, the mitochondrion and vectors are reported. When
most reads match a contaminant, the sample is taken to be that organism and the others are
still graded. Human is detected on reads of about 7 kb and up (see
[Human](methods.md#human)).

## Thresholds

`lrqc --print-thresholds` prints every level as TOML. A file given with `--thresholds`
overrides the keys it names; unknown keys, negative levels and a fail level that is not worse
than its warn level are errors.

```toml
contaminant_warn = 0.1
contaminant_fail = 1.0
read_q_warn_ont = 15.0
```

`--fail-on warn` or `--fail-on fail` makes lrQC exit with status 2 when any check of any
group reaches that level, naming the worst on standard error, so pipelines can stop on a bad
run. Status 1 means lrQC itself failed, including a mistyped option.

## Watching a run

`--watch SECONDS` reads files as they appear: each pass reads the files untouched for at
least that long, adds them to the totals and rewrites the outputs. Each output is written to a
temporary file and then renamed, so a browser reload never sees half a page, and while the run
goes on the page refreshes itself at the same interval.

lrQC stops when the run is over and every file is read: when MinKNOW writes
`final_summary_*.txt`, or a PacBio instrument writes `<movie>.transferdone`, in the run folder
(the folder watched, one of its subfolders, or one or two folders above it). It watches read
files, not the sequencing summary, which MinKNOW completes only at the end of the run.

- A file that cannot be read yet (still being copied) is read again on the next pass, without
  holding back the other files of that pass and without counting anything twice. A file still
  unreadable when the run ends is left out, and lrQC exits with an error naming it, after
  writing the report of everything else.
- Where MinKNOW writes both FASTQ and BAM, the BAM is read: if FASTQ was read before the first
  BAM appeared, the run is read again without it, so `--watch` and a one-off run agree.
- PacBio barcodes to search reads for come from the file names as files appear; if one appears
  after other files were read, every file is read again once the run finishes.

## Comparing reports

`lrqc compare` reads the `qc.json` of earlier reports (files or their directories) and writes
`compare.html`. It has one matrix of checks and one summary table over every group of every
report, overlays of each report's first group (all reads, or a one-barcode run's barcode), and
overlays of every barcode across reports. A `qc.json` from an earlier lrQC is compared where
the format allows, with a warning; one from a later lrQC is refused.

## Options

```
lrqc [OPTIONS] <INPUTS>...
lrqc compare [-o DIR] <REPORTS>...

  -o, --output <DIR>             output directory [default: lrqc]
  -t, --threads <N>              worker threads [default: all cores]
      --platform <P>             auto, ont or pacbio [default: auto, from headers and tags]
      --min-qscore <Q>           read Q at which reads pass when the input does not say
      --genome-size <SIZE>       haploid genome size (3.1g, 4.6m, 48502), for coverage
      --sample-sheet <CSV>       MinKNOW or SMRT Link sample sheet
      --barcode-names <FASTA>    lima's barcode FASTA: names bc indices, finds PacBio barcodes
      --adapters <FASTA>         extra adapter or primer sequences to search for
      --kinnex-adapters <FASTA>  skera's Kinnex/MAS adapters, to count array segments
      --thresholds <TOML>        override warn and fail levels
      --print-thresholds         print the default levels and exit
      --fail-on <LEVEL>          exit 2 when a check warns or fails
      --split-runs               report each run (flow cell or movie) and barcode apart
      --no-screen                skip the control and contaminant screen
      --per-read                 also write reads.tsv.gz, one row per read
      --watch <SECONDS>          follow a running sequencing run
      --title <TEXT>             report title [default: the first input's name]
      --max-reads <N>            stop after about N reads, for a quick look
```

## Performance

Measured on a 12th-generation Intel Core i5 (6 performance and 4 efficiency cores, 16
threads) with `--per-read`, reading from NVMe:

| Input | Size | Reads | Bases | Wall | Peak memory | `reads.tsv.gz` | `report.html` |
|---|---|---|---|---|---|---|---|
| HG002 ONT FASTQ, chr20 (BGZF) | 4.4 GB | 188,140 | 3.5 Gb | 17 s | 230 MB | 7 MB | 0.8 MB |
| HG002 RNA004 direct RNA, BAM | 35.5 GB | 16.5 M | 20.1 Gb | 93 s | 1.2 GB | 783 MB | 1.6 MB |
| HG002 HiFi (Revio-era CCS), unaligned BAM | 45.8 GB | 6.5 M | 101 Gb | 3 min 08 s | 1.2 GB | 185 MB | 1.0 MB |
| HG002 ONT chr1–6, aligned BAM | 82.4 GB | 3.2 M | 58.6 Gb | 3 min 24 s | 420 MB | 180 MB | 1.4 MB |

That is 250–400 MB of compressed input, or 200–540 Mb of sequence, per second. The drive was
not the limit (it reads at 1.1 GB/s). Throughput stops improving beyond about 8 threads on this
CPU: wall time is the same at `-t 8`, `-t 12` and `-t 16`. Memory grows with reads (6 bytes per
read for exact N50s, medians and the filter grid, 10–20 per read name for duplicate detection,
about as much per PacBio ZMW) and with groups (about 1 MB each per thread for per-position
counts). The
per-read table takes 28–56 bytes per read; the report stays near 1 MB per group.

Output is identical at any thread count apart from the recorded command line. On HiFi the screen
takes about a quarter of the CPU time; `--no-screen` saves it.

## Limits

- The channel map is in channel-number order, not the flow cell's physical layout.
- End reasons come from dorado's `er` tags (dorado 1.0 and later) or sequencing summaries;
  older BAM and FASTQ files do not carry them.
- MinKNOW's pore occupancy, mux scan, temperature and bias voltage come from files in the run
  folder (`pore_activity_*.csv`, `report_*.json`), which lrQC does not read.
- Kinnex segment counts need skera's adapter FASTA; lrQC has no built-in Kinnex adapters.
- `--max-reads` stops at batch boundaries, so it reads a little more than asked. It counts BAM
  records, secondary and supplementary ones included, and in watch mode applies to each pass.

## Development

```sh
cargo fmt
cargo clippy --all-targets -- -D warnings
cargo test
```

`tests/data` holds small fixtures cut from real HG003 ONT and HG002 PacBio HiFi reads
(dorado `BC` tags, lima split files, a MinKNOW folder, a sequencing summary), whose expected
counts come from samtools and awk. `tests/data/synthetic` holds constructed reads with known
answers: E. coli fragments among random reads, cDNA reads of yeast ENO2 with ONT's primers,
and Kinnex arrays built from skera's documented adapters.
