"""
Script to concatenate FASTQ files by group, based on filename patterns.

- Scans the specified directory for FASTQ files matching *_r?_*_s?.fq.gz.
- Groups files by sample and read (e.g., sample_r1_flowcell_s1.fq.gz).
- For each group, concatenates all matching files into a single output file.
- Supports a --dry-run mode to preview planned concatenations.

Assumed FASTQ filename structure:
    {sample}_r{read}_{flowcell}_s{lane}.fq.gz
    e.g., SampleA_r1_HJY7KDRXX_s1.fq.gz

    - sample:   Sample name (may include underscores)
    - r1/r2:    Read 1 or Read 2
    - flowcell: Flowcell ID
    - s1:       Lane number

Usage:
    python ci_fastq_concat.py [directory] [--dry-run]

Arguments:
    directory   Target directory to scan (default: current directory)
    --dry-run   Only show planned concatenations, do not perform them

Output:
    Concatenated FASTQ files for each group, named as per the group pattern.
"""

import os
import sys
import glob
from subprocess import call
import argparse

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Concatenate FASTQ files by group.")
parser.add_argument("directory", nargs="?", default=".", help="Target directory (default: current directory)")
parser.add_argument("--dry-run", action="store_true", help="Show planned concatenations without running them")
args = parser.parse_args()

# Change to the target directory
os.chdir(args.directory)

# Find all FASTQ files matching the pattern
fqs = glob.glob(r"*_r?_*_s?.fq.gz")

uniqes = set()
for fq in fqs:
    # Group by all but the last two underscore-separated fields (read, lane)
    name = "_".join(fq.split("_")[:-2])
    print(f"Processing file: {fq}, base name: {name}")
    if fq.endswith(".gz"):
        extens = fq.split(".")[-2].replace("fastq", "fq")+".gz"
        print(f"Detected gzipped file, adding unique {name+'...'+extens}")
    else:
        # error case, abort
        raise ValueError(f"File does not end with .gz: {fq}")
    uniqes.add(name+"..."+extens)

planned = []
for unique in uniqes:
    print(f"Processing group: {unique}")
    group_name = unique.split("...")[0]
    print(f"Group name: {group_name}")
    group_ext = unique.split("...")[-1]
    # Find all files in this group (matching prefix and extension)
    copies = sorted([fq for fq in fqs if fq.startswith(group_name) and fq.replace("fastq", "fq").endswith(group_ext)])
    outname = unique.replace("...", ".")
    arg_text = "cat "+ " ".join(copies) +" > "+ outname
    planned.append((copies, outname, arg_text))

# Print planned actions
print("Planned concatenations:")
for copies, outname, arg_text in planned:
    print(f"{' + '.join(copies)} -> {outname}")

print(f"\nFound {len(fqs)} FASTQ files in directory '{args.directory}':")
print(f"Total groups to concatenate: {len(planned)}")

if args.dry_run:
    print("\nDry run complete. No files were concatenated.")
else:
    for copies, outname, arg_text in planned:
        print(f"Running: {arg_text}")
        call([arg_text], shell=True)
    print("\nAll concatenations complete.")

