# Changelog

## 0.9.0 (2026-11-16)

First public release: per-barcode QC reports for Oxford Nanopore and PacBio long reads.

- Reads BAM (aligned or not), FASTQ and FASTA (gzipped or not), MinKNOW output folders,
  sequencing summaries and standard input; splits by barcode from dorado, MinKNOW, lima and
  SMRT Link, or by path.
- Writes a self-contained `report.html`, `summary.tsv`, `qc.json`, MultiQC custom content and,
  with `--per-read`, one row per read.
- Grades each barcode with pass, warn and fail checks whose levels can be changed with
  `--thresholds`; `--fail-on` turns them into an exit status for pipelines.
- Screens reads for spike-in controls (lambda DCS, yeast ENO2, ERCC, SIRV), common
  contaminants (phiX, E. coli, five Mycoplasma species, human) and vectors, and searches
  them for adapters, primers and every barcode in dorado's kits and PacBio's adapter indexes.
- Follows a run as it goes with `--watch`, and compares runs with `lrqc compare`.
- Tested on edge-case inputs and 27 real datasets (GIAB HG001/HG002/HG003 ONT and Revio HiFi,
  20 bacterial isolates, Candida, Drosophila, public ONT and PacBio 96-barcode runs), and
  benchmarked against eight long-read QC tools (see BENCHMARKS.md).
