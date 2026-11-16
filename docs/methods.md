# How lrQC works

The details behind the checks in the [README](../README.md): where lrQC departs from FastQC,
how the contaminant screen and the barcode search work, and how they were validated.

## Where lrQC departs from FastQC, and why

- **Positions.** Single bases up to 100, then bins 50, 500, 5,000 and 50,000 bases wide, and
  500,000 wide from 1 Mb to 10 Mb for ultra-long reads, measured from both the start and the
  end of reads. Bins reached by fewer than 0.5% of reads
  are not plotted.
- **Base content** is graded from position 21, and only in bins with 2,000 bases or more.
  Every long-read run has biased first and last bases, and the adapter check grades adapters
  left there. In small groups, single-base bins vary by chance more than the 10% warn level.
- **GC content** is graded on a second peak, not on departure from a normal curve. Reads of
  tens of kilobases average GC over isochores, so a clean human sample has a skewed single
  peak, and FastQC's test fails it. The departure is still shown. Peaks are judged on reads of
  500 bases or more: adapter and barcode fragments (about 45% GC) can be a fifth of a run's
  reads and would pass for a second organism; the adapter and overrepresented-sequence checks
  grade them instead.
- **Duplication and overrepresented sequences** use a sample of up to 10,000 reads per group:
  the reads with the smallest name hashes. The sample depends only on the read names, so it
  is the same at any thread count. Duplicates match on both 50-base ends.
- **Adapters** are found anywhere in a read by shared 12-mers on one alignment diagonal,
  covering at least a quarter of the adapter's informative 12-mers. Low-complexity 12-mers
  are not indexed, so (GA)n and similar repeats do not match. The built-in list is
  dorado's ligation, rapid, RNA004, cDNA, PCS110, MAB114 and DRB004 sequences, the Illumina
  read 1 and read 2 primers (10x, Twist UDI) and the PacBio SMRTbell adapter. Add your own with `--adapters extra.fasta`.
- **CpG methylation** comes from walking each `MM` entry's deltas over the read's C bases and
  checking for a following G. In implicit mode (`C+m.`, or no mode, as jasmine writes), the
  bases an entry skips count as unmodified calls, as the SAM specification and modkit have
  it; in explicit mode (`C+m?`) they were not called. dorado's all-context models (`5mC_5hmC`) call every C, so the
  all-context rate is low (about 8% in human); the CpG rate is the one to compare with the
  usual 70–80%.

## The contaminant screen

Every read is compared with fourteen references built into the binary. References and reads are
sampled the same way, keeping the canonical 19-mers whose hash falls in the lowest 1/24 of the
range (FracMinHash), so the index fits in the CPU cache and a read makes one lookup every 24
bases. The human genome is sampled at 1 in 1,000 instead, a subset of the same hashes (1.8
million k-mers: a 3 MB file, 7 MB in memory). Low-complexity 19-mers (fewer than 8 distinct trinucleotides of 17:
poly-A, (AT)n, (CTT)n) are left out everywhere. A read is assigned to the reference holding the
largest share of its sampled k-mers, at that reference's density, when that is at least 3
k-mers and a share of the read's that depends on its read Q: half the share a read
of that accuracy keeps, `(1 − error)^19`, and never less than 8%. That is 49% at Q30, 41% at
Q20 and 8% at Q10 and below. The index records which references hold each k-mer: E. coli K-12
carries lambda-like prophage remnants, and a read from those stretches goes to the reference it
matches best, or to none when it matches both equally.

| Reference | Accession | Kind |
|---|---|---|
| Lambda phage, ONT's DNA control strand (DCS) | NC_001416.1 | control |
| Yeast ENO2, ONT's RNA calibration strand | NC_001140.6:451327-452640 | control |
| PhiX174 | NC_001422.1 | contaminant |
| E. coli K-12 MG1655, rRNA and tRNA genes masked | NC_000913.3 | contaminant |
| Human mitochondrion | NC_012920.1 | organelle |
| ERCC RNA spike-ins (92 transcripts) | ERCC92 (Thermo Fisher) | control |
| Lexogen SIRV-Set 4: 69 SIRV isoforms on 7 loci and 15 long SIRVs | SIRV-Set 4 (2025-04-02) | control |
| *Mycoplasma hyorhinis* NCTC10130, rRNA and tRNA genes masked | GCF_900476065.1 | contaminant |
| *Mycoplasma orale* NCTC10112, likewise | GCF_900660435.1 | contaminant |
| *Mycoplasma fermentans* M64, likewise | GCF_000186005.1 | contaminant |
| *Mycoplasma arginini* HAZ145_1, likewise | GCF_001547975.1 | contaminant |
| *Acholeplasma laidlawii* PG-8A, likewise | GCF_000018785.1 | contaminant |
| Cloning vectors (UniVec Core), sequence shared with the references above masked | UniVec_Core | vector |
| Human nuclear genome (T2T-CHM13v2.0), 1 in 1,000 19-mers | GCF_009914755.1 | contaminant |

The five Mycoplasma species cause most cell-culture contamination. The spike-ins are controls:
reported, never graded. UniVec shares 29% of its bases with E. coli, lambda, phiX and the
spike-ins (vectors carry E. coli genes), and those are masked so the reads stay with the
reference they come from. UniVec also holds natural plasmid backbones (IncQ RSF1010, IncP RK2):
in ONT's 96-barcode run of bacterial isolates, up to 0.6% of a sample's reads matched them,
natural plasmids, not lab vectors. So vector hits are reported but not graded.

### Human

The human genome is too large to index densely, so it is sampled at 1 in 1,000, and reads are
tested for human only on 19-mers of middling G+C (7 to 14 of 19) with at least 10 distinct
trinucleotides:

- **Why the filter.** The AT-rich k-mers that fill human repeats turn up by chance in AT-rich
  genomes. With no filter, 1% of reads from a 27% GC Borrelia isolate were called human; with
  5 to 14 G+C, 14 reads (0.13% of bases) still were, on scattered single hits, two thirds of
  them on k-mers with 5 or 6 G+C. With 7 to 14, none are, and across two public 96-sample
  bacterial runs human calls fell from 642 reads to 94 (at most 0.09% of a sample), while
  calls on HG002 fell about 3%.
- **Recall depends on length.** A Q20 read keeps about 0.7 sampled human k-mers per kb. On
  HG002 ONT reads (Q20), 98.6% of reads of 8 kb and more were called human, 77% at 4–8 kb,
  39% at 2–4 kb and 4% below 2 kb; 98.4% of HG002 HiFi reads were.
- **What the share is of.** When reads long enough to test (6 or more testable k-mers, about
  7 kb) are at least a fifth of the run and number 100 or more, human's share is of them: on
  HG002 ONT, 98.8% of testable reads, so the sample is taken to be human. Otherwise its share
  is of all reads, a lower bound, since the short reads went untested. That share is graded
  but warns at most: it cannot tell a 1.5 kb amplicon run with 2% long human reads
  (contamination, flagged) from a human RNA run whose short reads are human too (on HG002
  RNA004 only 2% of reads were long enough, and 10% of all reads matched human: a warning
  that says why). The share of testable reads alone would call the amplicon run human, since
  its only long reads are the human ones.
- **Rebuilding.** `lrqc build-human-sketch genome.fna.gz assets/human.sketch.gz` rebuilds the
  sketch; mitochondrial records (`chrM`, `MT`, "mitochondrion") are left out, as the
  mitochondrion is its own reference.

### Related bacteria

Related bacteria share conserved sequence with E. coli, so a fixed bar calls them contaminated.
Tested on two public 96-barcode runs of known bacteria (PacBio's 2021 microbial 96-plex and
ONT's `nomiss_96BC_P2I_SUP_2026`), a fixed 8% bar found 6–14% "E. coli" in Salmonella,
Klebsiella, Citrobacter, Enterobacter, Cronobacter and Vibrio. Two changes remove that:

- **The quality-scaled bar.** A HiFi read of E. coli K-12 shares 99% of its sampled 19-mers
  with the reference (median); Salmonella and Klebsiella HiFi reads share under 30%.
- **Masking E. coli's 22 rRNA and 86 tRNA genes**, with 300 bases either side (70 kb, 1.5% of
  the genome; coordinates from NCBI's NC_000913.3 annotation). These genes are near-identical
  across bacteria, and short ONT reads can fall wholly inside them.

With both, every related species stays under the 1% warn level on both runs (Salmonella
0.6–0.8%, the rest lower), and both runs' E. coli samples (K-12, H10407 and the ONT isolates)
still match at 78–99%.

The bar assumes read Q is roughly right. A basecaller's predicted Q runs a point or two
optimistic (median Q20.2 against 98.6% identity on an HG002 ONT run), and the bar tolerates reads
up to about 6 points worse than their label. On constructed data with honest qualities it found
all 75 spiked lambda reads and 243 of 246 E. coli reads (the rest fall in masked genes), and
every 1 kb window of lambda at 4% errors. Reads shorter than about 300 bases keep too few
sampled k-mers to be assigned, and the 8% floor leaves out human reads holding a short nuclear
copy of mitochondrial DNA (a 1 kb NUMT in a 20 kb read shares about 5%).

On an HG002 RNA004 run the screen found 1.7% SIRV and 0.11% ERCC reads (both, as
SIRV-Set 4 holds ERCCs; synthetic sequence no human read matches by chance), 8.9% human mitochondrial reads, and no Mycoplasma: before the low-complexity
filter, 69 junk reads of poly-A and (AT)n (median Q4) matched AT-rich *M. hyorhinis*.

Only contaminants are graded. When most reads match a contaminant, the sample is taken to be
that organism and the check names it; the other contaminants are still graded, so an E. coli
isolate with phiX, or a human sample with Mycoplasma, still fails. The summary's
`contaminant_pct` is the share of the single worst contaminant other than the sample's own
organism. In a demultiplexed run (two or more barcodes), all reads and unclassified reads pool
every sample, so there the screen and the GC check (whose second peak means a mixed sample,
which a pool is by design) are shown without a verdict.
`--no-screen` turns the screen off.

## Barcodes inside reads

Each read is searched for barcode sequences. ONT barcodes are matched inside their kit's
flanks, as dorado matches them: every kit in dorado's `barcode_kits.cpp` (native, rapid,
PCR, PCR-cDNA, 16S, RLB/RPB, MAB114 amplicon, DRB004 direct RNA and Twist UDI), with flanks cut
to the 12 bases next to the barcode. A bare 24-base barcode is too short: one sequencing error
in its middle breaks every 12-mer but one, and a 15-base chance match among hundreds of
sequences turns up in a few percent of long reads. Twist's 10-base barcodes are matched
exactly. PacBio barcodes come from `--barcode-names`, or else from PacBio's SMRTbell adapter
indexes (16 of bc1001–bc1022, and bc2001–bc2384) for the barcodes the file names or sample sheet name;
they are matched exactly. A read's own barcode within 150 bases of an end means it was not trimmed;
another sample's barcode means cross-talk or a chimera across samples. In the unclassified reads
of a demultiplexed run, barcodes found are reported, not graded. On constructed data with 2% errors it found 97% of
untrimmed barcodes and 98% of foreign ones, with no false hits in 26,000 other reads.

## Direct RNA

Reads basecalled with a direct RNA model (dorado's `rna004_…` models, from read group names in
BAM or FASTQ, or MinKNOW's `basecall_model_version_id=` in FASTQ headers) are treated as RNA. Base content, GC content, duplication and overrepresented
sequences are shown without a verdict: single-stranded transcripts need not balance A with T or
G with C, transcriptome GC has several peaks, and highly expressed transcripts repeat by design.
The CpG methylation figure is left out, and modifications are named as RNA modifications (m6A,
m5C, inosine, pseudouridine, and the 2′-O-methyl Am, Cm, Gm and Um).
