#!/usr/bin/env python
"""Trinucleotide context counts of ONE sample over ONE region (or all of it), or of the reference itself.
Usage: get_context.py <sample|hg38> <region|total> <regions.json> --fasta FA --outdir DIR [--coverage-dir DIR] [--chroms chr1 ...]
Counts the trinucleotide centred on each covered position, on the + strand, raw 64 3-mers (NOT pyrimidine-collapsed),
weighted by duplex depth (col 4 of the coverage bed). Sums to compute_coverage's duplex_bp minus positions whose 3-mer
has an N or runs off a chromosome end (reported as n_invalid_duplex_bp).
sample = "hg38": the reference, every position of --chroms with depth 1 (no coverage needed).
region = "total": no region mask. Otherwise regions.json = {label: bed}, merged on the fly (overlaps counted once).
Writes <outdir>/<sample>.json"""
import argparse, io, itertools, json, os, subprocess, time
from pathlib import Path

import numpy as np
import pandas as pd

TRINUCS = ["".join(t) for t in itertools.product("ACGT", repeat=3)]       # index = 16*b0 + 4*b1 + b2
DEFAULT_CHROMS = [f"chr{i}" for i in list(range(1, 23)) + ["X"]]
LOOKUP = np.full(256, 4, dtype=np.uint8)
for i, b in enumerate("ACGT"):
    LOOKUP[ord(b)] = LOOKUP[ord(b.lower())] = i


def read_chrom(fasta, chrom):
    """One chromosome as base codes (A,C,G,T = 0-3, anything else 4), using the .fai index."""
    for line in open(f"{fasta}.fai"):
        name, length, offset, line_bases, line_bytes = line.split()[:5]
        if name == chrom:
            break
    else:
        raise KeyError(f"{chrom} not in {fasta}.fai")
    with open(fasta, "rb") as f:
        f.seek(int(offset))
        n_lines = -(-int(length) // int(line_bases))
        raw = f.read(int(length) + n_lines * (int(line_bytes) - int(line_bases)))
    seq = np.frombuffer(raw.replace(b"\n", b"").replace(b"\r", b""), dtype=np.uint8)[:int(length)]
    return LOOKUP[seq]


def trinuc_codes(bases):
    """Trinucleotide index (0-63) centred on each position; 64 if it has an N or runs off the chromosome end."""
    code = np.full(len(bases), 64, dtype=np.uint8)
    b0, b1, b2 = bases[:-2], bases[1:-1], bases[2:]
    ok = (b0 < 4) & (b1 < 4) & (b2 < 4)
    code[1:-1][ok] = 16 * b0[ok] + 4 * b1[ok] + b2[ok]
    return code


def region_mask(bed, chrom, length):
    """Boolean array: positions of `chrom` inside the bed (intervals merged by construction)."""
    df = pd.read_csv(bed, sep="\t", header=None, usecols=[0, 1, 2], dtype={0: str})
    df = df[df[0] == chrom]
    diff = np.bincount(df[1].values, minlength=length + 1)[:length + 1] - np.bincount(df[2].values, minlength=length + 1)[:length + 1]
    return np.cumsum(diff[:length]) > 0


def main(sample, region, regions_json, fasta, outdir, coverage_dir, chroms):
    t0 = time.time()
    bed = None if region == "total" else json.load(open(regions_json))[region]
    duplex = np.zeros(65)
    for chrom in chroms:
        code = trinuc_codes(read_chrom(fasta, chrom))
        if sample == "hg38":
            pos, depth = np.arange(len(code)), np.ones(len(code))
        else:
            cov = Path(coverage_dir) / f"{sample}_coverage.bed.gz"
            out = subprocess.run(["tabix", str(cov), chrom], capture_output=True, check=True).stdout
            df = pd.read_csv(io.BytesIO(out), sep="\t", header=None, usecols=[1, 2, 3], dtype=np.int64)
            lens = (df[2] - df[1]).values
            starts = np.repeat(df[1].values - (np.cumsum(lens) - lens), lens)
            pos = starts + np.arange(lens.sum())                      # every covered position, runs expanded
            depth = np.repeat(df[3].values, lens).astype(float)
        if bed is not None:
            keep = region_mask(bed, chrom, len(code))[pos]
            pos, depth = pos[keep], depth[keep]
        duplex += np.bincount(code[pos], weights=depth, minlength=65)
        print(f"{sample} {region} {chrom}: {time.time()-t0:.0f}s", flush=True)

    row = {"sample": sample, "region": region, "chroms": chroms,
           "duplex": dict(zip(TRINUCS, duplex[:64].tolist())), "n_invalid_duplex_bp": duplex[64]}
    out = Path(outdir) / f"{sample}.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_suffix(".tmp")
    json.dump(row, open(tmp, "w"))
    os.replace(tmp, out)                                              # atomic: an existing .json is a finished task
    print(f"{sample} {region}: duplex {duplex[:64].sum():,.0f} (+{duplex[64]:,.0f} invalid) ({time.time()-t0:.0f}s)", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("sample"); ap.add_argument("region"); ap.add_argument("regions_json")
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--coverage-dir", default="../dupcaller_pipe_outputs/coverage")
    ap.add_argument("--chroms", nargs="+", default=DEFAULT_CHROMS)
    a = ap.parse_args()
    main(a.sample, a.region, a.regions_json, a.fasta, a.outdir, a.coverage_dir, a.chroms)
