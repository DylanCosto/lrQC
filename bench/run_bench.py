#!/usr/bin/env python3
"""Times lrqc and other long-read QC tools on the same inputs.

usage: run_bench.py DATASETS.json OUTDIR [--threads 8] [--reps 3] [--tools lrqc,nanoq,...]

DATASETS.json is a list of {"name": ..., "path": ..., "format": "fastq" | "ubam" | "bam"}.
Each tool runs with the same thread count, one run at a time; wall time and peak memory
come from GNU time. A tool whose first run takes over 5 minutes runs once. Results go to
OUTDIR/runs.jsonl, one line per run, and each run's outputs to OUTDIR/<dataset>/<tool>/<rep>.

Tools are looked up on PATH; put the competitors in one environment (see bench/README.md).
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
import time


def commands(tool, fmt, path, out, threads, lrqc):
    """The command for one tool on one input, or None when the tool cannot read it."""
    t = str(threads)
    if tool == "lrqc":
        return [lrqc, path, "-o", out, "-t", t]
    if tool == "lrqc_noscreen":
        return [lrqc, path, "-o", out, "-t", t, "--no-screen"]
    if tool == "nanoplot":
        flag = {"fastq": "--fastq", "ubam": "--ubam", "bam": "--bam"}[fmt]
        return ["NanoPlot", flag, path, "-o", out, "-t", t]
    if tool == "nanostat":
        flag = {"fastq": "--fastq", "ubam": "--ubam", "bam": "--bam"}[fmt]
        return ["NanoStat", flag, path, "-o", out, "-n", "NanoStats.txt", "-t", t]
    if tool == "nanoq" and fmt == "fastq":
        return ["sh", "-c", f'nanoq -i "{path}" -s -j > "{out}/nanoq.json"']
    if tool == "seqkit" and fmt == "fastq":
        return ["sh", "-c", f'seqkit stats -a -T -j {t} "{path}" > "{out}/seqkit.tsv"']
    if tool == "cramino" and fmt in ("ubam", "bam"):
        extra = "--ubam " if fmt == "ubam" else ""
        return ["sh", "-c", f'cramino {extra}-t {t} --format json "{path}" > "{out}/cramino.json"']
    if tool == "fastqc" and fmt in ("fastq", "ubam", "bam"):
        kind = {"fastq": "fastq", "ubam": "bam", "bam": "bam_mapped"}[fmt]
        return ["fastqc", "-t", t, "--memory", "10000", "-f", kind, "--extract", "--outdir", out, path]
    if tool == "longreadsum":
        sub = {"fastq": "fq", "ubam": "bam", "bam": "bam"}[fmt]
        return ["longreadsum", sub, "-i", path, "-o", out, "-t", t, "-g", f"{out}/longreadsum.log"]
    if tool == "fastplong" and fmt == "fastq":
        return [
            "fastplong", "-i", path, "-o", "/dev/null", "-w", t,
            "-A", "-Q", "-L", "--disable_chimera_splitting",
            "-j", f"{out}/fastplong.json", "-h", f"{out}/fastplong.html",
        ]
    return None


TOOLS = ["lrqc", "lrqc_noscreen", "nanoq", "seqkit", "cramino", "fastplong", "longreadsum", "nanostat", "nanoplot", "fastqc"]


def timed(cmd):
    """Runs cmd under GNU time; returns (seconds, peak KB, exit status)."""
    stats = cmd_out = None
    with open(os.devnull, "w") as null:
        proc = subprocess.run(
            ["/usr/bin/time", "-f", "TIME %e %M"] + cmd,
            stdout=null,
            stderr=subprocess.PIPE,
            text=True,
        )
    for line in proc.stderr.splitlines()[::-1]:
        if line.startswith("TIME "):
            _, secs, kb = line.split()
            stats = (float(secs), int(kb))
            break
    cmd_out = proc.stderr[-2000:]
    return stats or (None, None), proc.returncode, cmd_out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("datasets")
    ap.add_argument("outdir")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--reps", type=int, default=3)
    ap.add_argument("--tools", default=",".join(TOOLS))
    ap.add_argument("--lrqc", default="lrqc")
    args = ap.parse_args()
    datasets = json.load(open(args.datasets))
    os.makedirs(args.outdir, exist_ok=True)
    log = open(os.path.join(args.outdir, "runs.jsonl"), "a")
    for ds in datasets:
        for tool in args.tools.split(","):
            for rep in range(args.reps):
                out = os.path.join(args.outdir, ds["name"], tool, str(rep))
                cmd = commands(tool, ds["format"], ds["path"], out, args.threads, args.lrqc)
                if cmd is None:
                    break
                shutil.rmtree(out, ignore_errors=True)
                os.makedirs(out)
                load = os.getloadavg()[0]
                (secs, kb), rc, err = timed(cmd)
                record = {
                    "dataset": ds["name"], "tool": tool, "rep": rep, "wall_s": secs,
                    "peak_mb": kb / 1024 if kb else None, "exit": rc, "load1": load,
                    "threads": args.threads, "cmd": cmd, "time": time.strftime("%H:%M:%S"),
                }
                if rc != 0:
                    record["stderr"] = err
                log.write(json.dumps(record) + "\n")
                log.flush()
                print(f"{ds['name']:22} {tool:12} rep {rep} {secs}s {kb and kb // 1024} MB exit {rc}", flush=True)
                # Keep the first run's outputs for the accuracy comparison; later reps only time.
                if rep > 0:
                    shutil.rmtree(out, ignore_errors=True)
                if rc != 0 or (secs or 0) > 300:
                    break


if __name__ == "__main__":
    sys.exit(main())
