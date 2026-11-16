# Benchmarks

lrQC 0.9.0 (as 0.1.0, its release candidate) against eight long-read QC tools on five real datasets, for speed and for the
accuracy of the numbers they report. Measured 2026-11-16; scripts in [bench/](bench/).

## Setup

- **Machine**: 2 × Intel Xeon Gold 5418Y (24 cores), 121 GB RAM, Linux. The machine was
  shared with another job (load average 1.5–20 during the runs), so absolute times carry some
  noise; every tool saw the same conditions, and repeated runs agreed within about 2%.
- **Tools** (bioconda): NanoPlot 1.46.2, NanoStat 1.6.0, nanoq 0.10.0, cramino 2.0.0,
  FastQC 0.13.0, seqkit 2.14.0, LongReadSum 1.6.0, fastplong 0.7.1.
- **Threads**: 8 for every tool that takes a thread count (`-t 8`, `-j 8`, `-w 8`).
  nanoq is single-threaded.
- **Runs**: three per tool and dataset; wall time is the median, memory the peak (GNU time).
- **Outputs**: each tool's default report (lrQC: `report.html`, `summary.tsv`, `qc.json`,
  MultiQC files; NanoPlot: its plots and HTML; FastQC: its HTML and zip). `lrqc_noscreen` is
  lrQC with `--no-screen`, for comparison with tools that do no contaminant screening.

| Dataset | Input | Size | Reads | Bases |
|---|---|---|---|---|
| `hg002_ont_fastq` | GIAB HG002 ONT R10.4.1 chr20, FASTQ (plain gzip) | 2.8 GB | 188,140 | 3.52 Gb |
| `hg002_hifi_fastq` | GIAB HG002 Revio HiFi chr20, FASTQ (plain gzip, no `rq`) | 1.5 GB | 195,097 | 3.15 Gb |
| `hg002_hifi_fastq_uncompressed` | the same, uncompressed | 6.3 GB | 195,097 | 3.15 Gb |
| `vibrio_ont_fastq` | *Vibrio* isolate, ONT (SRA), FASTQ (plain gzip) | 3.8 GB | 794,720 | 4.40 Gb |
| `hg001_revio_ubam` | GIAB HG001 Revio HiFi, unaligned BAM (first 150,000 reads) | 1.2 GB | 149,994 | 2.74 Gb |
| `hg001_ont_aligned_bam` | GIAB HG001 ONT dorado sup, aligned BAM with MM/ML (first 250,000 records) | 4.9 GB | 249,801 | 4.38 Gb |

**Truth** comes from [bench/truth.py](bench/truth.py), which reads every record with numpy
(and samtools for BAM, primary records only) and computes each number from its definition.

## Speed

Wall time (median of 3) and peak memory:

| tool | hg002 ONT FASTQ | hg002 HiFi FASTQ | Vibrio ONT FASTQ | HG001 HiFi uBAM | HG001 ONT aligned BAM | HiFi FASTQ, uncompressed |
|---|---|---|---|---|---|---|
| **lrQC** | 25.6 s, 108 MB | 16.0 s, 112 MB | 35.6 s, 115 MB | 9.2 s, 217 MB | 17.9 s, 265 MB | 6.7 s, 175 MB |
| lrQC, no screen | 25.2 s, 65 MB | 15.8 s, 87 MB | 35.4 s, 91 MB | 7.1 s, 192 MB | 14.7 s, 241 MB | 4.7 s, 149 MB |
| cramino | – | – | – | 37.1 s, 12 MB | 8.4 s, 16 MB | – |
| fastplong | 27.2 s, 4,453 MB | 14.7 s, 455 MB | 27.2 s, 1,131 MB | – | – | 9.3 s, 930 MB |
| seqkit | 51.8 s, 30 MB | 31.5 s, 26 MB | 60.6 s, 28 MB | – | – | 11.8 s, 22 MB |
| nanoq | 63.5 s, 9 MB | 47.1 s, 7 MB | 79.9 s, 17 MB | – | – | 28.1 s, 7 MB |
| fastqc | 99.9 s, 8,406 MB | 69.6 s, 5,457 MB | 88.0 s, 5,885 MB | 66.9 s, 7,130 MB | 104.6 s, 4,534 MB | 68.6 s, 6,401 MB |
| nanostat | 143.1 s, 171 MB | 124.0 s, 173 MB | 181.5 s, 305 MB | 76.1 s, 259 MB | 223.3 s, 349 MB | 95.0 s, 265 MB |
| nanoplot | 151.8 s, 256 MB | 132.4 s, 256 MB | 189.9 s, 285 MB | 84.3 s, 275 MB | 251.2 s, 345 MB | 103.5 s, 256 MB |
| longreadsum | 0.7 s, 326 MB † | 0.9 s, 419 MB † | 0.9 s, 423 MB † | 13.2 s, 5,531 MB | 35.6 s, 5,543 MB | 88.1 s, 419 MB |

– : the tool does not read this input (cramino reads only BAM/CRAM; seqkit, nanoq and fastplong
only FASTQ/FASTA). † LongReadSum stops after one read of a gzipped FASTQ (see below), so its
time there is not a full run.

- lrQC ran a full FastQC-style analysis per read, with the contaminant screen, adapter and
  barcode search, in 7–36 s: 5–15× faster than NanoPlot and NanoStat, 2.5–10× faster than
  FastQC, and faster than the summary-only tools seqkit (1.7–2×) and nanoq (2–4×) on every
  input. fastplong (statistics and its own report, trimming off) is within 10% of lrQC on the
  human FASTQs and faster on the Vibrio run of many short reads (27 s against 36 s). cramino,
  which reports summary numbers of BAM records, is faster on the aligned BAM (8 s against
  18 s) and 4× slower on the unaligned one.
- Plain gzip decompresses on one thread, which bounds lrQC on gzipped FASTQ (16 s gzipped,
  6.7 s uncompressed for the same HiFi reads). BGZF (`bgzip`) input decompresses on all cores.
- The screen costs 0.2–3.2 s of these times (`lrqc_noscreen`), most on PacBio and BAM input.
- lrQC used 65–275 MB; FastQC 4.5–8.4 GB, fastplong up to 4.5 GB, LongReadSum 5.5 GB on BAM.

## Accuracy

Every number each tool printed was checked against the truth: reads, bases, mean and
median length, N50, longest read, GC, read Q (median, mean, share at Q20) and bases at Q20.
Read Q is compared under each tool's own definition: the mean of per-read Q (lrQC, nanoq), the
Q of the reads' mean error (NanoPlot, NanoStat), the Q of the mean error over all bases
(seqkit). lrQC uses the basecaller's Q (`qs`, PacBio `rq`) where the input has it, as dorado
and SMRT Link report it; the others compute Q from base qualities.

| tool | values checked | exact | rounded as printed | wrong |
|---|---|---|---|---|
| **lrQC** | 60 | 60 | 0 | 0 |
| lrQC, no screen | 60 | 60 | 0 | 0 |
| nanoplot | 54 | 54 | 0 | 0 |
| nanostat | 54 | 54 | 0 | 0 |
| nanoq | 36 | 33 | 3 | 0 |
| seqkit | 36 | 32 | 4 | 0 |
| cramino | 10 | 9 | 0 | 1 |
| fastplong | 20 | 17 | 3 | 0 |
| fastqc | 36 | 17 | 15 | 4 |
| longreadsum | 33 | 21 | 0 | 12 |

What the wrong values are:

- **LongReadSum** reads one record of a gzipped FASTQ and reports it as the whole file (1
  read, 292 bases for 195,097 reads and 3.15 Gb), with no error; on the ONT FASTQ its JSON
  could not be read at all. On uncompressed FASTQ and on BAM its numbers are right. Its JSON
  also writes `-nan` where a value is undefined, which JSON parsers reject.
- **cramino** reports aligned bases as the yield of an aligned BAM (4.09 Gb of 4.38 Gb of read
  bases), a different definition rather than an error.
- **FastQC** on an aligned BAM (`-f bam_mapped`) counts only the aligned part of each read, so
  its lengths are those of alignments (longest 113,771 bases for a 682,961-base read). It
  prints bases to two digits and GC as a whole percent.
- NanoPlot and NanoStat print Q to one decimal; nanoq and fastplong truncate mean lengths to a
  whole base. These are marked as rounding, not errors.

lrQC's numbers matched the truth on every input. Two lrQC bugs this benchmark found were fixed
before these runs: the share of reads at Q20 had counted reads of Q19.95–19.99 (from read Qs
rounded to tenths), and the medians of read length and read Q took the upper or lower middle
value for an even count instead of the mean of the two.

## Caveats

- Times are from one shared machine; on an idle machine all tools would be somewhat faster.
- The tools do different amounts of work: NanoPlot draws many plots, FastQC runs its modules
  on every base, seqkit and nanoq compute summary statistics only. The comparison is of what
  each tool produces by default.
- Truth uses the same definitions as the tools it checks; where definitions differ (read Q
  means, aligned yield), each tool is checked against its own.

## Per-dataset results
### hg002_ont_fastq

Truth: 188,140 reads, 3,519,431,361 bases.

| metric | truth | lrQC | lrQC, no screen | nanoplot | nanostat | nanoq | seqkit | fastqc | longreadsum | fastplong |
|---|---|---|---|---|---|---|---|---|---|---|
| reads | 188,140 | 188,140 | 188,140 | 188,140 | 188,140 | 188,140 | 188,140 | 188,140 | – | 188,140 |
| bases | 3,519,431,361 | 3,519,431,361 | 3,519,431,361 | 3,519,431,361 | 3,519,431,361 | 3,519,431,361 | 3,519,431,361 | 3,500,000,000 ≈ | – | 3,519,431,361 |
| mean length | 18,706.45 | 18,706.40 | 18,706.40 | 18,706.40 | 18,706.40 | 18,706 | 18,706.40 | 18,706 | – | 18,706 |
| median length | 17,875 | 17,875 | 17,875 | 17,875 | 17,875 | 17,875 | 17,875 | 17,875 | – | – |
| N50 | 29,411 | 29,411 | 29,411 | 29,411 | 29,411 | 29,411 | 29,411 | – | – | – |
| longest | 629,276 | 629,276 | 629,276 | 629,276 | 629,276 | 629,276 | 629,276 | 629,276 | – | – |
| GC % | 43.31 | 43.31 | 43.31 | – | – | – | 43.31 | 43 ≈ | – | 43.31 |
| median read Q | 25.48 (from bases) | 25.50 | 25.50 | 25.50 | 25.50 | 25.48 | – | – | – | – |
| mean read Q | 24.17 (from bases) | 24.17 | 24.17 | 18.50 | 18.50 | 24.17 | 18.57 | – | – | – |
| reads ≥Q20 % | 79.61 (from bases) | 79.61 | 79.61 | 79.60 | 79.60 | 79.61 | – | – | – | – |
| bases ≥Q20 % | 91.32 | – | – | – | – | – | 91 ≈ | – | – | 91.32 |

- longreadsum: output unreadable (IndexError: list index out of range)

| tool | wall (s, median of runs) | runs | peak memory (MB) |
|---|---|---|---|
| lrQC | 25.6 | 3 | 108 |
| lrQC, no screen | 25.2 | 3 | 65 |
| nanoplot | 151.8 | 3 | 256 |
| nanostat | 143.1 | 3 | 171 |
| nanoq | 63.5 | 3 | 9 |
| seqkit | 51.8 | 3 | 30 |
| fastqc | 99.9 | 3 | 8,406 |
| longreadsum | 0.7 | 3 | 326 |
| fastplong | 27.2 | 3 | 4,453 |

### hg002_hifi_fastq

Truth: 195,097 reads, 3,153,107,709 bases.

| metric | truth | lrQC | lrQC, no screen | nanoplot | nanostat | nanoq | seqkit | fastqc | longreadsum | fastplong |
|---|---|---|---|---|---|---|---|---|---|---|
| reads | 195,097 | 195,097 | 195,097 | 195,097 | 195,097 | 195,097 | 195,097 | 195,097 | 1 ✗ (195,097) | 195,097 |
| bases | 3,153,107,709 | 3,153,107,709 | 3,153,107,709 | 3,153,107,709 | 3,153,107,709 | 3,153,107,709 | 3,153,107,709 | 3,100,000,000 ≈ | 292 ✗ (3,153,107,709) | 3,153,107,709 |
| mean length | 16,161.74 | 16,161.70 | 16,161.70 | 16,161.70 | 16,161.70 | 16,161 ≈ | 16,161.70 | 16,161 ≈ | 292 ✗ (16,161.74) | 16,161 ≈ |
| median length | 16,030 | 16,030 | 16,030 | 16,030 | 16,030 | 16,030 | 16,030 | 16,030 | 292 ✗ (16,030) | – |
| N50 | 17,705 | 17,705 | 17,705 | 17,705 | 17,705 | 17,705 | 17,705 | – | 292 ✗ (17,705) | – |
| longest | 42,750 | 42,750 | 42,750 | 42,750 | 42,750 | 42,750 | 42,750 | 42,750 | 292 ✗ (42,750) | – |
| GC % | 43.10 | 43.10 | 43.10 | – | – | – | 43.10 | 43 ≈ | – | 43.10 |
| median read Q | 30.47 (from bases) | 30.50 | 30.50 | 30.50 | 30.50 | 30.47 | – | – | – | – |
| mean read Q | 29.97 (from bases) | 29.97 | 29.97 | 25.80 | 25.80 | 29.97 | 25.24 | – | – | – |
| reads ≥Q20 % | 93.31 (from bases) | 93.31 | 93.31 | 93.30 | 93.30 | 93.31 | – | – | – | – |
| bases ≥Q20 % | 97.57 | – | – | – | – | – | 98 ≈ | – | – | 97.57 |

| tool | wall (s, median of runs) | runs | peak memory (MB) |
|---|---|---|---|
| lrQC | 16.0 | 3 | 112 |
| lrQC, no screen | 15.8 | 3 | 87 |
| nanoplot | 132.4 | 3 | 256 |
| nanostat | 124.0 | 3 | 173 |
| nanoq | 47.1 | 3 | 7 |
| seqkit | 31.5 | 3 | 26 |
| fastqc | 69.6 | 3 | 5,457 |
| longreadsum | 0.9 | 3 | 419 |
| fastplong | 14.7 | 3 | 455 |

### vibrio_ont_fastq

Truth: 794,720 reads, 4,400,886,127 bases.

| metric | truth | lrQC | lrQC, no screen | nanoplot | nanostat | nanoq | seqkit | fastqc | longreadsum | fastplong |
|---|---|---|---|---|---|---|---|---|---|---|
| reads | 794,720 | 794,720 | 794,720 | 794,720 | 794,720 | 794,720 | 794,720 | 794,720 | 1 ✗ (794,720) | 794,720 |
| bases | 4,400,886,127 | 4,400,886,127 | 4,400,886,127 | 4,400,886,127 | 4,400,886,127 | 4,400,886,127 | 4,400,886,127 | 4,400,000,000 ≈ | 87 ✗ (4,400,886,127) | 4,400,886,127 |
| mean length | 5,537.66 | 5,537.70 | 5,537.70 | 5,537.70 | 5,537.70 | 5,537 ≈ | 5,537.70 | 5,537 ≈ | 87 ✗ (5,537.66) | 5,537 ≈ |
| median length | 3,566 | 3,566 | 3,566 | 3,566 | 3,566 | 3,566 | 3,566 | 3,566 | 87 ✗ (3,566) | – |
| N50 | 8,257 | 8,257 | 8,257 | 8,257 | 8,257 | 8,257 | 8,257 | – | 87 ✗ (8,257) | – |
| longest | 162,343 | 162,343 | 162,343 | 162,343 | 162,343 | 162,343 | 162,343 | 162,343 | 87 ✗ (162,343) | – |
| GC % | 46.22 | 46.22 | 46.22 | – | – | – | 46.22 | 46 ≈ | – | 46.22 |
| median read Q | 19.10 (from bases) | 19.10 | 19.10 | 19.10 | 19.10 | 19.10 | – | – | – | – |
| mean read Q | 18.63 (from bases) | 18.63 | 18.63 | 17.40 | 17.40 | 18.63 | 17.09 | – | – | – |
| reads ≥Q20 % | 35.52 (from bases) | 35.52 | 35.52 | 35.50 | 35.50 | 35.52 | – | – | – | – |
| bases ≥Q20 % | 85.22 | – | – | – | – | – | 85 ≈ | – | – | 85.22 |

| tool | wall (s, median of runs) | runs | peak memory (MB) |
|---|---|---|---|
| lrQC | 35.6 | 3 | 115 |
| lrQC, no screen | 35.4 | 3 | 91 |
| nanoplot | 189.9 | 3 | 285 |
| nanostat | 181.5 | 3 | 305 |
| nanoq | 79.9 | 3 | 17 |
| seqkit | 60.6 | 3 | 28 |
| fastqc | 88.0 | 3 | 5,885 |
| longreadsum | 0.9 | 3 | 423 |
| fastplong | 27.2 | 3 | 1,131 |

### hg001_revio_ubam

Truth: 149,994 reads, 2,737,774,882 bases.

| metric | truth | lrQC | lrQC, no screen | nanoplot | nanostat | cramino | fastqc | longreadsum |
|---|---|---|---|---|---|---|---|---|
| reads | 149,994 | 149,994 | 149,994 | 149,994 | 149,994 | 149,994 | 149,994 | 149,994 |
| bases | 2,737,774,882 | 2,737,774,882 | 2,737,774,882 | 2,737,774,882 | 2,737,774,882 | 2,737,774,882 | 2,700,000,000 ≈ | 2,737,774,882 |
| mean length | 18,252.56 | 18,252.60 | 18,252.60 | 18,252.60 | 18,252.60 | 18,252.56 | 18,252 ≈ | 18,252.56 |
| median length | 17,878 | 17,878 | 17,878 | 17,878 | 17,878 | 17,878 | 17,878 | 17,878 |
| N50 | 19,244 | 19,244 | 19,244 | 19,244 | 19,244 | 19,244 | – | 19,244 |
| longest | 43,523 | 43,523 | 43,523 | 43,523 | 43,523 | – | 43,523 | 43,523 |
| GC % | 40.16 | 40.16 | 40.16 | – | – | – | 40 ≈ | 40.16 |
| median read Q | 32.23 (basecaller's) / 27.88 (from bases) | 32.20 | 32.20 | 27.90 | 27.90 | – | – | – |
| mean read Q | 32.52 (basecaller's) / 27.90 (from bases) | 32.52 | 32.52 | 24.40 | 24.40 | – | – | – |
| reads ≥Q20 % | 100 (basecaller's) / 90.28 (from bases) | 100 | 100 | 90.30 | 90.30 | – | – | – |
| bases ≥Q20 % | 96.75 | – | – | – | – | – | – | – |

| tool | wall (s, median of runs) | runs | peak memory (MB) |
|---|---|---|---|
| lrQC | 9.2 | 3 | 217 |
| lrQC, no screen | 7.1 | 3 | 192 |
| nanoplot | 84.3 | 3 | 275 |
| nanostat | 76.1 | 3 | 259 |
| cramino | 37.1 | 3 | 12 |
| fastqc | 66.9 | 3 | 7,130 |
| longreadsum | 13.2 | 3 | 5,531 |

### hg001_ont_aligned_bam

Truth: 249,801 reads, 4,375,698,632 bases.

| metric | truth | lrQC | lrQC, no screen | nanoplot | nanostat | cramino | fastqc | longreadsum |
|---|---|---|---|---|---|---|---|---|
| reads | 249,801 | 249,801 | 249,801 | 249,801 | 249,801 | 249,801 | 249,801 | 249,801 |
| bases | 4,375,698,632 | 4,375,698,632 | 4,375,698,632 | 4,375,698,632 | 4,375,698,632 | 4,091,143,643 ✗ (4,375,698,632) | 4,000,000,000 ✗ (4,375,698,632) | 4,375,698,632 |
| mean length | 17,516.74 | 17,516.70 | 17,516.70 | 17,516.70 | 17,516.70 | 17,516.74 | 16,377 ✗ (17,516.74) | 17,516.74 |
| median length | 14,391 | 14,391 | 14,391 | 14,391 | 14,391 | 14,391 | 12,667 ✗ (14,391) | 14,391 |
| N50 | 28,337 | 28,337 | 28,337 | 28,337 | 28,337 | 28,337 | – | 28,337 |
| longest | 682,961 | 682,961 | 682,961 | 682,961 | 682,961 | – | 113,771 ✗ (682,961) | 682,961 |
| GC % | 42.24 | 42.24 | 42.24 | – | – | – | 42 ≈ | 42.24 |
| median read Q | 25.61 (basecaller's) / 25.21 (from bases) | 25.60 | 25.60 | 25.20 | 25.20 | – | – | – |
| mean read Q | 24.68 (basecaller's) / 24.15 (from bases) | 24.68 | 24.68 | 19.20 | 19.20 | – | – | – |
| reads ≥Q20 % | 82.83 (basecaller's) / 81.55 (from bases) | 82.83 | 82.83 | 81.50 | 81.50 | – | – | – |
| bases ≥Q20 % | 93.15 | – | – | – | – | – | – | – |

| tool | wall (s, median of runs) | runs | peak memory (MB) |
|---|---|---|---|
| lrQC | 17.9 | 3 | 265 |
| lrQC, no screen | 14.7 | 3 | 241 |
| nanoplot | 251.2 | 3 | 345 |
| nanostat | 223.3 | 3 | 349 |
| cramino | 8.4 | 3 | 16 |
| fastqc | 104.6 | 3 | 4,534 |
| longreadsum | 35.6 | 3 | 5,543 |

### hg002_hifi_fastq_uncompressed

Truth: 195,097 reads, 3,153,107,709 bases.

| metric | truth | lrQC | lrQC, no screen | nanoplot | nanostat | nanoq | seqkit | fastqc | longreadsum | fastplong |
|---|---|---|---|---|---|---|---|---|---|---|
| reads | 195,097 | 195,097 | 195,097 | 195,097 | 195,097 | 195,097 | 195,097 | 195,097 | 195,097 | 195,097 |
| bases | 3,153,107,709 | 3,153,107,709 | 3,153,107,709 | 3,153,107,709 | 3,153,107,709 | 3,153,107,709 | 3,153,107,709 | 3,100,000,000 ≈ | 3,153,107,709 | 3,153,107,709 |
| mean length | 16,161.74 | 16,161.70 | 16,161.70 | 16,161.70 | 16,161.70 | 16,161 ≈ | 16,161.70 | 16,161 ≈ | 16,161.74 | 16,161 ≈ |
| median length | 16,030 | 16,030 | 16,030 | 16,030 | 16,030 | 16,030 | 16,030 | 16,030 | 16,030 | – |
| N50 | 17,705 | 17,705 | 17,705 | 17,705 | 17,705 | 17,705 | 17,705 | – | 17,705 | – |
| longest | 42,750 | 42,750 | 42,750 | 42,750 | 42,750 | 42,750 | 42,750 | 42,750 | 42,750 | – |
| GC % | 43.10 | 43.10 | 43.10 | – | – | – | 43.10 | 43 ≈ | 43.10 | 43.10 |
| median read Q | 30.47 (from bases) | 30.50 | 30.50 | 30.50 | 30.50 | 30.47 | – | – | – | – |
| mean read Q | 29.97 (from bases) | 29.97 | 29.97 | 25.80 | 25.80 | 29.97 | 25.24 | – | – | – |
| reads ≥Q20 % | 93.31 (from bases) | 93.31 | 93.31 | 93.30 | 93.30 | 93.31 | – | – | – | – |
| bases ≥Q20 % | 97.57 | – | – | – | – | – | 98 ≈ | – | – | 97.57 |

| tool | wall (s, median of runs) | runs | peak memory (MB) |
|---|---|---|---|
| lrQC | 6.7 | 3 | 175 |
| lrQC, no screen | 4.7 | 3 | 149 |
| nanoplot | 103.5 | 3 | 256 |
| nanostat | 95.0 | 3 | 265 |
| nanoq | 28.1 | 3 | 7 |
| seqkit | 11.8 | 3 | 22 |
| fastqc | 68.6 | 3 | 6,401 |
| longreadsum | 88.1 | 3 | 419 |
| fastplong | 9.3 | 3 | 930 |

