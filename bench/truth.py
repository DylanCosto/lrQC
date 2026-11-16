#!/usr/bin/env python3
"""Reference values for the benchmark, computed independently of every QC tool.

usage: truth.py INPUT [--bam]   (FASTQ, gzipped or not, or BAM via samtools; needs numpy)

Prints JSON: reads, bases, mean and median length (the standard median: the mean of the
two middle lengths for an even count), N50, longest, GC percent of A/C/G/T bases, and read Q
two ways: from the basecaller's tag (`qs`, or PacBio `rq`) where present, and from base
qualities as -10*log10(mean error probability) (dorado's and NanoPlot's definition).
BAM counts primary records only (-F 0x900).
"""
import json
import subprocess
import sys

import numpy as np

ERROR = 10 ** (-np.arange(256, dtype=np.float64) / 10)
GC = np.zeros(256, dtype=np.int64)
ACGT = np.zeros(256, dtype=np.int64)
for b in b"GCgc":
    GC[b] = 1
for b in b"ACGTacgtUu":
    ACGT[b] = 1


def phred(err):
    return -10 * np.log10(np.maximum(err, 1e-6))


def records(path, bam):
    if bam:
        proc = subprocess.Popen(
            ["samtools", "view", "-@4", "-F", "0x900", path], stdout=subprocess.PIPE
        )
        for line in proc.stdout:
            f = line.rstrip(b"\n").split(b"\t")
            seq = b"" if f[9] == b"*" else f[9]
            qual = None if f[10] == b"*" else f[10]
            tag = None
            for t in f[11:]:
                if t.startswith(b"qs:f:"):
                    tag = float(t[5:])
                elif t.startswith(b"rq:f:"):
                    rq = float(t[5:])
                    tag = 60.0 if rq >= 1 else float(phred(1 - rq))
            yield seq, qual, tag
    else:
        proc = subprocess.Popen(["sh", "-c", f'zcat -f "{path}"'], stdout=subprocess.PIPE)
        lines = iter(proc.stdout)
        for header in lines:
            seq = next(lines).rstrip(b"\r\n")
            next(lines)
            qual = next(lines).rstrip(b"\r\n")
            tag = None
            for t in header.split()[1:]:
                if t.startswith(b"qs:f:"):
                    tag = float(t[5:])
                elif t.startswith(b"rq:f:"):
                    rq = float(t[5:])
                    tag = 60.0 if rq >= 1 else float(phred(1 - rq))
            yield seq, qual, tag


def main():
    path = sys.argv[1]
    bam = "--bam" in sys.argv or path.endswith(".bam")
    lengths, q_bases, q_tags = [], [], []
    gc = acgt = bases_q20 = qual_bases = 0
    error_sum = 0.0
    for seq, qual, tag in records(path, bam):
        lengths.append(len(seq))
        if seq:
            s = np.frombuffer(seq, dtype=np.uint8)
            gc += int(GC[s].sum())
            acgt += int(ACGT[s].sum())
        if qual:
            q = np.frombuffer(qual, dtype=np.uint8) - 33
            errors = ERROR[q]
            q_bases.append(float(phred(errors.mean())))
            error_sum += float(errors.sum())
            bases_q20 += int((q >= 20).sum())
            qual_bases += len(q)
        if tag is not None:
            q_tags.append(tag)
    lengths = np.sort(np.array(lengths, dtype=np.int64))
    total = int(lengths.sum())
    desc = lengths[::-1]
    cumulative = np.cumsum(desc)
    n50 = int(desc[np.searchsorted(cumulative * 2, total)]) if total else 0
    out = {
        "reads": int(len(lengths)),
        "bases": total,
        "mean_length": total / len(lengths) if len(lengths) else 0,
        "median_length": float(np.median(lengths)) if len(lengths) else 0,
        "n50": n50,
        "longest": int(lengths[-1]) if len(lengths) else 0,
        "gc_pct": 100 * gc / acgt if acgt else None,
        "bases_q20_pct": 100 * bases_q20 / qual_bases if qual_bases else None,
        # The Q of the mean error over all bases (lrqc's base_q; base-weighted).
        "base_error_q": float(phred(error_sum / qual_bases)) if qual_bases else None,
    }
    for name, values in (("q_from_bases", q_bases), ("q_from_tag", q_tags)):
        if values:
            v = np.array(values)
            out[name] = {
                "reads": int(len(v)),
                "mean": float(v.mean()),
                "median": float(np.median(v)),
                "reads_q10_pct": float(100 * (v >= 10).mean()),
                "reads_q20_pct": float(100 * (v >= 20).mean()),
                # The Q of the reads' mean error probability (NanoPlot's mean read quality).
                "q_of_mean_error": float(phred((10 ** (-v / 10)).mean())),
            }
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
