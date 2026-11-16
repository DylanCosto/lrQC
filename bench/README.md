# Benchmarks

Scripts that time lrQC against other long-read QC tools and check every tool's headline
numbers against values computed independently. Results are in [BENCHMARKS.md](../BENCHMARKS.md).

## Setup

The other tools, in one environment:

```sh
micromamba create -n lrqc-bench -c conda-forge -c bioconda \
    nanoplot nanostat nanoq cramino fastqc seqkit longreadsum fastplong multiqc numpy
```

and `samtools` on `PATH` (for BAM truth).

## Running

List the inputs in a JSON file:

```json
[{"name": "hg002_ont_fastq", "path": "/data/HG002.chr20.ont.fastq.gz", "format": "fastq"},
 {"name": "hg001_revio_ubam", "path": "/data/m84036.hifi_reads.bam", "format": "ubam"},
 {"name": "hg001_ont_bam", "path": "/data/HG001.sorted.bam", "format": "bam"}]
```

`format` is `fastq`, `ubam` (unaligned BAM) or `bam` (aligned). Then:

```sh
micromamba activate lrqc-bench
python bench/run_bench.py datasets.json results --threads 8 --reps 3 --lrqc target/release/lrqc
mkdir -p truth
for each dataset: python bench/truth.py <path> > truth/<name>.json
python bench/compare.py datasets.json truth results > results.md
```

`run_bench.py` runs one tool at a time with the same thread count, takes wall time and peak
memory from GNU time, and repeats each run (`--reps`) unless the first takes over 5 minutes.
Only the first run's outputs are kept, for `compare.py`.

`truth.py` reads every record (BAM: primary records, `samtools view -F 0x900`) and computes
read and base counts, mean and median length (the standard median), N50, longest read, GC
percent, and read Q both from the basecaller's tag (`qs`, PacBio `rq`) and from base
qualities, with the means tools use: the mean of per-read Q (lrQC, nanoq), the Q of the mean
per-read error (NanoPlot, NanoStat), and the Q of the mean error over all bases (seqkit's
`AvgQual`). `compare.py` checks each tool against the definition that tool uses, and marks a
value that differs beyond its printed precision with ✗ and the true value.
