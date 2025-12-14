# updated 15 Aug 2025 to include dry run and skip the header of the csv file

"""
Script to rename FASTQ files based on sample information in .contents.csv files.

- Scans the specified directory for *.contents.csv files.
- For each CSV, parses sample and barcode info to construct new filenames.
- Checks for missing files and naming conflicts.
- Supports a --dry-run mode to preview planned renames.

Assumed file name structure:
    CSV filename:      SLX12345.HJY7KDRXX.s_1.contents.csv
        - SLX12345:    Sequencer run ID
        - HJY7KDRXX:   Flowcell ID
        - s_1:         Lane number

    FASTQ filename:    SLX12345.xGenUDI30.HJY7KDRXX.s_1.r_1.fq.gz
        - xGenUDI30:      Barcode sequence or ID
        - r_1/r_2:     Read 1 or Read 2

Assumed CSV structure:
    - Comma-separated, with a header line (skipped by script)
    - At least 4 columns per row:
        [1]: barcode, [3]: sample name
    - Example row: 1,xGenUDI30,ACTGACTG-ACTGACTG,Sample-Name

Usage:
    python genomics_core_rename.py [directory] [--dry-run]

Arguments:
    directory   Target directory to scan (default: current directory)
    --dry-run   Only show planned renames, do not perform them

Typical CSV filename: SLX12345.HJY7KDRXX.s_1.contents.csv
Typical FASTQ filename: SLX12345.xGenUDI30.HJY7KDRXX.s_1.r_1.fq.gz
"""

import os
import sys
import argparse

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Rename FASTQ files based on .contents.csv files.")
parser.add_argument("directory", nargs="?", default=".", help="Target directory (default: current directory)")
parser.add_argument("--dry-run", action="store_true", help="Show planned renames without making changes")
args = parser.parse_args()

dry_run = args.dry_run
target_dir = args.directory

# Change to the target directory
os.chdir(target_dir)

# Find all .contents.csv files in the directory
csvs = [file for file in os.listdir("./") if file.endswith(".contents.csv")]

problems = []
planned = []

for csv in csvs:
    with open(csv) as table:
        next(table)  # skip header
        # Extract SLX, flowcell, and lane from the CSV filename
        slx = csv.split(".")[0]
        flowcell = csv.split(".")[1]
        lane = csv.split("s_")[-1].split(".")[0]

        for line in table:
            cols = line.strip().split(",")
            if len(cols) < 4:
                problems.append(f"Malformed line in {csv}: {line.strip()}")
                continue
            barcode = cols[1]
            sample = cols[3].replace("-", "_")
            # Construct old and new filenames for R1 and R2
            old_name_r1 = f"{slx}.{barcode}.{flowcell}.s_{lane}.r_1.fq.gz"
            old_name_r2 = f"{slx}.{barcode}.{flowcell}.s_{lane}.r_2.fq.gz"
            new_name1 = f"{sample}_r1_{flowcell}_s{lane}.fq.gz"
            new_name2 = f"{sample}_r2_{flowcell}_s{lane}.fq.gz"

            # Check for problems and plan renames
            if not os.path.isfile(old_name_r1):
                problems.append(f"Missing file: {old_name_r1}")
            elif os.path.exists(new_name1):
                problems.append(f"Target file already exists: {new_name1}")
            else:
                planned.append((old_name_r1, new_name1))

            if not os.path.isfile(old_name_r2):
                problems.append(f"Missing file: {old_name_r2}")
            elif os.path.exists(new_name2):
                problems.append(f"Target file already exists: {new_name2}")
            else:
                planned.append((old_name_r2, new_name2))

# Print planned actions and problems
print("Planned renames:")
for old, new in planned:
    print(f"{old} --> {new}")

if problems:
    print("\nProblems detected:")
    for p in problems:
        print(p)
    print("\nNo files will be renamed until problems are resolved.")

if not dry_run and not problems:
    for old, new in planned:
        os.rename(old, new)
        print(f"Renamed {old} --> {new}")
elif dry_run:
    print("\nDry run complete. No files were renamed.")
