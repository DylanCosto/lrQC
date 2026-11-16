#!/usr/bin/env python3
"""Compares each tool's headline numbers with truth.py's, and summarizes run times.

usage: compare.py DATASETS.json TRUTH_DIR RESULTS_DIR > results.md

Reads the first run's outputs of each tool (RESULTS_DIR/<dataset>/<tool>/0), the truth in
TRUTH_DIR/<dataset>.json, and timings in RESULTS_DIR/runs.jsonl.
"""
import csv
import glob
import json
import math
import os
import re
import statistics
import sys


def loose_json(path):
    text = open(path).read()
    # LongReadSum writes -nan and nan, which JSON does not allow.
    text = re.sub(r"(?<=[:\s\[,])-?nan\b", "null", text)
    return json.loads(text)


def num(text):
    return float(str(text).replace(",", ""))


def si_bases(text):
    """FastQC's '3.5 Gbp'."""
    value, unit = text.split()
    scale = {"bp": 1, "kbp": 1e3, "Mbp": 1e6, "Gbp": 1e9, "Tbp": 1e12}[unit]
    return float(value) * scale


def lrqc(out):
    row = next(csv.DictReader(open(f"{out}/summary.tsv"), delimiter="\t"))
    f = lambda k: None if row[k] == "NA" else float(row[k])
    return {
        "reads": f("reads"), "bases": f("bases"), "mean_length": f("mean_length"),
        "median_length": f("median_length"), "n50": f("n50"), "longest": f("longest"),
        "gc_pct": f("gc_pct"), "median_q": f("median_read_q"), "mean_q": f("mean_read_q"),
        "reads_q20_pct": f("reads_q20_pct"),
    }


def nanostats(path):
    text = open(path).read()
    get = lambda label: num(re.search(rf"^{label}:\s+([\d,.]+)", text, re.M).group(1))
    longest = re.search(r"Top 5 longest reads.*?\n1:\s+(\d+)", text, re.S)
    q20 = re.search(r"^>Q20:\s+\d+ \(([\d.]+)%\)", text, re.M)
    return {
        "reads": get("Number of reads"), "bases": get("Total bases"),
        "mean_length": get("Mean read length"), "median_length": get("Median read length"),
        "n50": get("Read length N50"), "longest": float(longest.group(1)) if longest else None,
        "median_q": get("Median read quality"), "mean_q": get("Mean read quality"),
        "reads_q20_pct": float(q20.group(1)) if q20 else None,
    }


def nanoq(out):
    d = json.load(open(f"{out}/nanoq.json"))
    return {
        "reads": d["reads"], "bases": d["bases"], "mean_length": d["mean_length"],
        "median_length": d["median_length"], "n50": d["n50"], "longest": d["longest"],
        "median_q": d["median_quality"], "mean_q": d["mean_quality"],
        "reads_q20_pct": 100 * d["quality_thresholds"]["20"] / d["reads"] if d["reads"] else None,
    }


def seqkit(out):
    row = next(csv.DictReader(open(f"{out}/seqkit.tsv"), delimiter="\t"))
    return {
        "reads": num(row["num_seqs"]), "bases": num(row["sum_len"]),
        "mean_length": num(row["avg_len"]), "median_length": num(row["Q2"]),
        "n50": num(row["N50"]), "longest": num(row["max_len"]), "gc_pct": num(row["GC(%)"]),
        "mean_q": num(row["AvgQual"]), "bases_q20_pct": num(row["Q20(%)"]),
    }


def cramino(out):
    d = json.load(open(f"{out}/cramino.json"))
    r = d["read_stats"]
    return {
        "reads": d["alignment_stats"]["num_reads"], "bases": round(r["yield_gb"] * 1e9),
        "mean_length": r["mean_length"], "median_length": r["median_length"], "n50": r["n50"],
    }


def fastqc(out):
    path = glob.glob(f"{out}/*_fastqc/fastqc_data.txt")[0]
    rows = dict(
        line.rstrip("\n").split("\t", 1)
        for line in open(path)
        if "\t" in line and not line.startswith(("#", ">>"))
    )
    lengths = rows["Sequence length"].split("-")
    return {
        "reads": num(rows["Total Sequences"]), "bases": si_bases(rows["Total Bases"]),
        "mean_length": num(rows["Mean Length"]) if "Mean Length" in rows else None,
        "median_length": num(rows["Median Length"]) if "Median Length" in rows else None,
        "longest": num(lengths[-1]), "gc_pct": num(rows["%GC"]),
    }


def longreadsum(out):
    path = glob.glob(f"{out}/*_summary.*.json")[0]
    d = loose_json(path)
    return {
        "reads": d["total_num_reads"], "bases": d["total_num_bases"],
        "mean_length": d["mean_read_length"], "median_length": d["median_read_length"],
        "n50": d["n50_read_length"], "longest": d["longest_read_length"],
        "gc_pct": d["gc_percent"],
    }


def fastplong(out):
    d = json.load(open(f"{out}/fastplong.json"))["summary"]["before_filtering"]
    return {
        "reads": d["total_reads"], "bases": d["total_bases"],
        "mean_length": d["read_mean_length"], "gc_pct": 100 * d["gc_content"],
        "bases_q20_pct": 100 * d["q20_rate"],
    }


PARSERS = {
    "lrqc": lrqc, "lrqc_noscreen": lrqc, "nanoplot": lambda o: nanostats(f"{o}/NanoStats.txt"),
    "nanostat": lambda o: nanostats(f"{o}/NanoStats.txt"), "nanoq": nanoq, "seqkit": seqkit,
    "cramino": cramino, "fastqc": fastqc, "longreadsum": longreadsum, "fastplong": fastplong,
}

# Metric, label, tolerance: a value is right when within it of the truth. Lengths and
# counts are exact up to the rounding tools print with; GC and Q to their printed digits.
METRICS = [
    ("reads", "reads", 0), ("bases", "bases", 0), ("mean_length", "mean length", 0.5),
    ("median_length", "median length", 0.5), ("n50", "N50", 0), ("longest", "longest", 0),
    ("gc_pct", "GC %", 0.06), ("median_q", "median read Q", 0.06),
    ("mean_q", "mean read Q", 0.06), ("reads_q20_pct", "reads ≥Q20 %", 0.06),
    ("bases_q20_pct", "bases ≥Q20 %", 0.06),
]

Q_METRICS = ("median_q", "mean_q", "reads_q20_pct")

# Values a tool prints rounded, so a difference within the rounding is not an error: the
# allowed difference, absolute, or relative when it ends in "rel".
PRINTED = {
    ("fastqc", "bases"): (0.05, "rel"),  # "3.5 Gbp"
    ("fastqc", "gc_pct"): (0.5, "abs"),  # whole percent
    ("fastqc", "mean_length"): (1.0, "abs"),
    ("nanoq", "mean_length"): (1.0, "abs"),  # truncated to a whole base
    ("fastplong", "mean_length"): (1.0, "abs"),
    ("seqkit", "bases_q20_pct"): (0.5, "abs"),  # whole percent
}


def within_printing(tool, metric, value, truth):
    allowed = PRINTED.get((tool, metric))
    if allowed is None:
        return False
    limit, kind = allowed
    diff = abs(value - truth)
    return diff <= (limit * abs(truth) if kind == "rel" else limit) + 1e-9

# What each tool's read Q figures are, as truth.py keys: (read Q source, mean definition).
# lrqc takes the basecaller's Q where there is one; NanoPlot's mean read quality is the Q of
# the reads' mean error; seqkit's AvgQual is the Q of the mean error over all bases.
def q_truth(tool, metric, truth):
    tag = truth.get("q_from_tag")
    bases = truth.get("q_from_bases")
    source = tag if (tool.startswith("lrqc") and tag) else bases
    if source is None:
        return None
    if metric == "median_q":
        return source["median"]
    if metric == "reads_q20_pct":
        return source["reads_q20_pct"]
    if tool in ("nanoplot", "nanostat"):
        return source["q_of_mean_error"]
    if tool == "seqkit":
        return truth.get("base_error_q")
    return source["mean"]


def truth_label(metric, truth):
    tag, bases = truth.get("q_from_tag"), truth.get("q_from_bases")
    key = {"median_q": "median", "mean_q": "mean", "reads_q20_pct": "reads_q20_pct"}[metric]
    parts = []
    if tag:
        parts.append(f"{fmt(tag[key])} (basecaller's)")
    if bases:
        parts.append(f"{fmt(bases[key])} (from bases)")
    return " / ".join(parts) or "–"


def fmt(v):
    if v is None:
        return "–"
    if isinstance(v, float) and not v.is_integer():
        return f"{v:,.2f}"
    return f"{int(v):,}"


def main():
    datasets, truth_dir, results = sys.argv[1:4]
    runs = [json.loads(l) for l in open(f"{results}/runs.jsonl")]
    print("# Benchmark results\n")
    for ds in json.load(open(datasets)):
        name = ds["name"]
        truth = json.load(open(f"{truth_dir}/{name}.json"))
        # Tools that use the basecaller's read Q compare against it; the rest against the
        # Q of the bases. lrqc uses the tag when there is one.
        print(f"## {name}\n")
        print(f"Truth: {fmt(truth['reads'])} reads, {fmt(truth['bases'])} bases.\n")
        tools = [
            t for t in PARSERS
            if any(r["dataset"] == name and r["tool"] == t for r in runs)
        ]
        values = {}
        for tool in tools:
            ok = [r for r in runs if r["dataset"] == name and r["tool"] == tool and r["exit"] == 0]
            if not ok:
                values[tool] = None
                continue
            try:
                values[tool] = PARSERS[tool](f"{results}/{name}/{tool}/0")
            except Exception as e:  # noqa: BLE001 - a tool's output we cannot read
                values[tool] = {"error": f"{type(e).__name__}: {e}"}
        print("| metric | truth | " + " | ".join(tools) + " |")
        print("|---|---|" + "---|" * len(tools))
        for metric, label, tol in METRICS:
            if metric in Q_METRICS:
                if not (truth.get("q_from_tag") or truth.get("q_from_bases")):
                    continue
                t_label = truth_label(metric, truth)
            else:
                if truth.get(metric) is None:
                    continue
                t_label = fmt(float(truth[metric]))
            cells = []
            for tool in tools:
                v = (values[tool] or {}).get(metric)
                t = q_truth(tool, metric, truth) if metric in Q_METRICS else truth.get(metric)
                if v is None or t is None:
                    cells.append("–")
                    continue
                right = abs(float(v) - t) <= tol + 1e-9
                if right:
                    cells.append(fmt(float(v)))
                elif within_printing(tool, metric, float(v), t):
                    cells.append(fmt(float(v)) + " ≈")
                else:
                    cells.append(fmt(float(v)) + f" ✗ ({fmt(float(t))})")
            print(f"| {label} | {t_label} | " + " | ".join(cells) + " |")
        errors = {t: v["error"] for t, v in values.items() if v and "error" in v}
        failed = [t for t, v in values.items() if v is None]
        if errors or failed:
            print()
        for tool, err in errors.items():
            print(f"- {tool}: output unreadable ({err})")
        for tool in failed:
            print(f"- {tool}: run failed")
        print("\n| tool | wall (s, median of runs) | runs | peak memory (MB) |")
        print("|---|---|---|---|")
        for tool in tools:
            ok = [r for r in runs if r["dataset"] == name and r["tool"] == tool and r["exit"] == 0]
            if not ok:
                print(f"| {tool} | failed | | |")
                continue
            wall = statistics.median(r["wall_s"] for r in ok)
            peak = max(r["peak_mb"] for r in ok)
            print(f"| {tool} | {wall:,.1f} | {len(ok)} | {peak:,.0f} |")
        print()


if __name__ == "__main__":
    main()
