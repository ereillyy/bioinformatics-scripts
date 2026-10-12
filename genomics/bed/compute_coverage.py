#!/usr/bin/env python
"""Duplex coverage of ONE sample over ONE region (or the whole-genome total).
Usage: compute_coverage.py <sample> <region|total> <regions.json> [--outdir DIR] [--coverage-dir DIR]
regions.json = {label: raw bed path}; the bed is sorted + merged on the fly (process substitution).
Writes <outdir>/<sample>/<region>.json"""
import argparse, json, os, subprocess, time
from pathlib import Path

# footprint = bp covered at all; duplex_bp = depth-weighted (col 4) -- the mutation-burden denominator
SUM_BP = "awk '{fp+=$3-$2; depth+=$4} END{print fp+0, depth+0, NR}'"


def run(cmd):
    r = subprocess.run(f"set -o pipefail; {cmd}", shell=True, executable="/bin/bash",
                       capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"Command failed:\n{cmd}\n\nstderr:\n{r.stderr}")
    return [int(x) for x in r.stdout.split()]


def main(sample, region, regions_json, outdir, coverage_dir):
    t0 = time.time()
    cov = Path(coverage_dir) / f"{sample}_coverage.bed.gz"
    row = {"sample": sample, "region": region}

    if region == "total":
        row["footprint_bp"], row["duplex_bp"], _ = run(f"zcat {cov} | {SUM_BP}")
    else:
        bed = json.load(open(regions_json))[region]
        # same ordering as the coverage bed (lexicographic chrom, then start) so -sorted is valid
        merged = f"LC_ALL=C sort -k1,1 -k2,2n {bed} | bedtools merge -i -"
        row["n_merged_intervals"] = run(f"{merged} | wc -l")[0]
        if row["n_merged_intervals"] == 0:
            raise RuntimeError(f"{region}: merged region bed is empty ({bed})")
        row["footprint_bp"], row["duplex_bp"], n_rows = run(
            f"bedtools intersect -a <(zcat {cov}) -b <({merged}) -sorted | {SUM_BP}")
        if n_rows == 0:  # not expected for any real region; fail loudly rather than record a silent 0
            raise RuntimeError(f"{sample} {region}: intersect returned no rows")

    out = Path(outdir) / sample / f"{region}.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_suffix(".tmp")
    json.dump(row, open(tmp, "w"))
    os.replace(tmp, out)  # atomic: an existing .json always means a finished task
    print(f"{sample} {region}: {row['footprint_bp']:,} / {row['duplex_bp']:,} ({time.time()-t0:.0f}s)", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("sample"); ap.add_argument("region"); ap.add_argument("regions_json")
    ap.add_argument("--outdir", default="notebook_data/deeptools/coverage_per_sample_v2")
    ap.add_argument("--coverage-dir", default="../dupcaller_pipe_outputs/coverage")
    a = ap.parse_args()
    main(a.sample, a.region, a.regions_json, a.outdir, a.coverage_dir)
