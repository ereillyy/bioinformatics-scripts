#!/usr/bin/env python3

"""
Concatenate files by either:
- exact filename match (from different directories), or
- prefix grouping (by N '.'-separated fields),
keeping R1/R2 separate and submitting each concat as a Slurm job.

- Output filenames are constructed from the group key and detected read, with configurable formatting.
- Supports a --dry-run mode.

Usage:
    python general_concat.py "<input_glob>" <outdir> --group-mode {exact,prefix} [--prefix-fields N] [--read-style {r,R,r_,R_}] [--dry-run]

Arguments:
    <input_glob>        Glob pattern for input files (in quotes)
    <outdir>            Output directory (created if missing)
    --group-mode        'exact' (by filename) or 'prefix' (by N fields, default: exact)
    --prefix-fields     Number of '.'-separated fields to use as prefix for grouping (required if --group-mode prefix)
    --read-style        Output read style: 'r' (default, e.g. r1), 'R' (R1), 'r_' (r_1), 'R_' (R_1)
    --dry-run           Only print planned actions, do not submit jobs
"""

import os
import sys
import glob
import argparse
from collections import defaultdict
import subprocess
import re

def detect_read(filename):
    # Detects R1/R2/r1/r2/R_1/r_1 etc, returns 'r1' or 'r2'
    m = re.search(r'([._-])([rR])[_]?([12])(\D|$)', filename)
    if m:
        return m.group(2).lower() + m.group(3)  # e.g. r1 or r2
    raise ValueError(f"Could not detect read (R1/R2/r1/r2/R_1/r_1) in filename: {filename}")

def format_read(read, style):
    # read: 'r1' or 'r2'
    if style == 'r':
        return read
    elif style == 'R':
        return read.upper()
    elif style == 'r_':
        return read[0] + '_' + read[1]
    elif style == 'R_':
        return read[0].upper() + '_' + read[1]
    else:
        raise ValueError(f"Unknown read style: {style}")

parser = argparse.ArgumentParser(description="Concatenate files by exact name or prefix, submit to Slurm.")
parser.add_argument("input_glob", help="Glob pattern for input files (in quotes)")
parser.add_argument("outdir", help="Output directory")
parser.add_argument("--group-mode", choices=['exact', 'prefix'], default='exact', help="Group by exact filename or prefix fields (default: exact)")
parser.add_argument("--prefix-fields", type=int, help="Number of '.'-separated fields to use as prefix for grouping (required if --group-mode prefix)")
parser.add_argument("--read-style", choices=['r', 'R', 'r_', 'R_'], default='r', help="Output read style (default: r, e.g. r1)")
parser.add_argument("--dry-run", action="store_true", help="Show planned actions, do not submit jobs")
args = parser.parse_args()

input_files = glob.glob(args.input_glob)
if not input_files:
    print(f"No files found for input_glob: {args.input_glob}")
    sys.exit(1)

# Check all files are .fq.gz or .fastq.gz
for f in input_files:
    if not (f.endswith('.fq.gz') or f.endswith('.fastq.gz')):
        print(f"Error: File does not end with .fq.gz or .fastq.gz: {f}")
        sys.exit(1)

groups = defaultdict(lambda: {"files": [], "outfile_name": None})
if args.group_mode == 'exact':
    for f in input_files:
        base = os.path.basename(f)
        try:
            read = detect_read(base)
        except ValueError as e:
            print(e)
            continue
        # Always output as .fq.gz
        if base.endswith('.fastq.gz'):
            outfile_name = base[:-9] + '.fq.gz'
        else:
            outfile_name = base
        key = (base, read)
        groups[key]["files"].append(f)
        groups[key]["outfile_name"] = outfile_name
elif args.group_mode == 'prefix':
    if not args.prefix_fields:
        print("Error: --prefix-fields N must be specified when using --group-mode prefix, corresponding to number of '.'-separated fields to use as prefix")
        sys.exit(1)
    for f in input_files:
        base = os.path.basename(f)
        fields = base.split('.')
        if len(fields) < args.prefix_fields:
            print(f"Skipping file (not enough fields): {base}")
            continue
        prefix = '.'.join(fields[:args.prefix_fields])
        try:
            read = detect_read(base)
        except ValueError as e:
            print(e)
            continue
        # Always output as .fq.gz
        outfile_name = prefix + '_' + read + ".fq.gz"
        key = (prefix, read)
        groups[key]["files"].append(f)
        groups[key]["outfile_name"] = outfile_name

print(f"Found {len(groups)} unique groups.")

# Prepare output directory
if not os.path.exists(args.outdir):
    print(f"Creating output directory: {args.outdir}")
    if not args.dry_run:
        os.makedirs(args.outdir, exist_ok=True)

for (group_key, read), groupinfo in groups.items():
    files = groupinfo["files"]
    orig_outfile_name = groupinfo["outfile_name"]
    # Reformat the read part in outfile_name according to --read-style
    m = re.search(r'(r[_]?1|R[_]?1|r[_]?2|R[_]?2)', orig_outfile_name)
    if m:
        read_num = m.group(0)[-1]
        new_read = format_read('r' + read_num, args.read_style)
        out_base = orig_outfile_name[:m.start()] + new_read + orig_outfile_name[m.end():]
    else:
        # fallback: just append .fq.gz if not found
        out_base = orig_outfile_name
    out_path = os.path.join(args.outdir, out_base)
    files_str = " ".join(sorted(files))
    print(f"\nWill concatenate {len(files)} files into {out_path}:")
    for f in sorted(files):
        print(f"  {f}")
    cmd = f"cat {files_str} > {out_path}"
    slurm_script = f"""#!/bin/bash
#SBATCH --job-name=cat_{str(group_key).replace('.','_')}_{read}
#SBATCH --output={os.path.join(args.outdir, 'slurmlogs', out_base)}.slurm.out
#SBATCH --mem=12G
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00

set -euo pipefail
echo "Concatenating to {out_path}"
echo "Running command: {cmd}"
{cmd}
"""
    if args.dry_run:
        print(f"[DRY RUN] Would submit Slurm job for: {out_base}")
        print(slurm_script)
    else:
        import tempfile
        with tempfile.NamedTemporaryFile("w", delete=False, suffix=".sh") as tf:
            tf.write(slurm_script)
            tf.flush()
            print(f"Submitting Slurm job for {out_base} with script {tf.name}")
            subprocess.run(["sbatch", tf.name])
