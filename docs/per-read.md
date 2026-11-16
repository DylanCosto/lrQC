# The per-read table

`--per-read` writes `reads.tsv.gz`: a header line, then one tab-separated row per read (primary
records only), in input order whatever the thread count. It is BGZF, so `zcat`, `gzip -dc`,
R's `read.delim`, pandas' `read_csv(sep="\t")` and `polars.read_csv` read it directly. Missing
values are `NA`; list columns are comma-separated and empty when there is nothing.

| Column | What |
|---|---|
| `read_id` | the read name as in the input |
| `group` | barcode, or `run · barcode` with `--split-runs` |
| `run` | flow cell, run ID or PacBio movie |
| `platform` | `ont`, `pacbio`, or `NA` when no header or tag says |
| `length` | bases |
| `read_q` | read Q, as in the report |
| `filter` | `pass` or `fail`; `NA` for a read with neither a pass/fail source nor a read Q |
| `gc_pct` | GC percent of the read |
| `channel`, `mux` | ONT `ch`, `mx` |
| `start_time` | ISO 8601 UTC from BAM and FASTQ (`st`); seconds from the run start from summaries |
| `duration` | ONT `du`, seconds |
| `duplex` | `simplex`, `duplex` or `simplex_with_duplex_offspring` (`dx`) |
| `split` | 1 if dorado split the read from a longer one (`pi`) |
| `passes`, `accuracy`, `effective_coverage`, `zmw`, `barcode_score` | PacBio `np`, `rq`, `ec`, `zm`, `bq` |
| `poly_a` | dorado poly(A) estimate (`pt`) |
| `screen` | the screen reference the read matched: `lambda`, `eno2`, `ercc`, `sirv`, `phix`, `ecoli`, `chrM`, `m_hyorhinis`, `m_orale`, `m_fermentans`, `m_arginini`, `a_laidlawii`, `univec` or `human`; empty when it matched none, `NA` when it was not screened |
| `adapters` | adapters and primers found, in read order |
| `barcodes_inside` | barcodes found inside the read |
| `kinnex_adapters` | number of Kinnex adapters found, with `--kinnex-adapters` |
| `mod_calls` | number of ML values (modification calls) |
| `mapq`, `identity` | primary alignment MAPQ (`unmapped` if unmapped) and BLAST-like identity, in percent |

Workers compress the rows of each batch, so the table costs little: on the HG002 ONT FASTQ
below it comes to 7 MB for 188,140 reads. In watch mode each pass appends to the table (gzip
members concatenate), with one header.
