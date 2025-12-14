#!/bin/bash

# -----------------------------------------------------------------------------
# Script to demultiplex paired-end FASTQ files using 'idemux' and SLURM.
#
# - Scans a directory for R1 FASTQ files matching *.r_1.lostreads.fq.gz.
# - For each R1, finds the corresponding R2 file (*.r_2.lostreads.fq.gz).
# - Submits a SLURM job to run 'idemux' for each pair, writing output to a subdirectory of the output directory.
# - Supports a --dry-run mode to preview jobs without submitting.
#
# Assumed FASTQ filename structure:
#   {sample}.r_1.lostreads.fq.gz   (R1)
#   {sample}.r_2.lostreads.fq.gz   (R2)
#   - {sample} can be any string (e.g., SLX12345.xGenUDI30.HJY7KDRXX.s_1)
#
# Assumed sample sheet (CSV) structure:
#   - Has a header line (skipped).
#   - At least 4 columns per row:
#       [1]: barcode, [3]: sample name
#   - Example row: 1,xGenUDI30,,Sample-Name
#
# Usage:
#   ./genomics_core_demultiplex.sh <fq_dir> <sample_sheet> <output_dir> [--dry-run]
#
# Example:
#   ./genomics_core_demultiplex.sh data/SLX-26143 data/SLX-26143/sample_barcodes_slx26143.csv data/SLX-26143/idemux --dry-run
#
# Requirements:
#   - 'idemux' must be installed and available in your PATH.
# -----------------------------------------------------------------------------

set -euo pipefail

DRYRUN=0

if [[ $# -lt 3 || $# -gt 4 ]]; then
    echo "Usage: $0 <fq_dir> <sample_sheet> <output_dir> [--dry-run]"
    exit 1
fi

FQDIR="$1"
SAMPLESHEET="$2"
OUTDIR="$3"

if [[ $# -eq 4 ]]; then
    echo
    if [[ "$4" == "--dry-run" ]]; then
        DRYRUN=1
        echo "Dry run mode enabled: jobs will not be submitted."
    else
        echo "Unknown option: $4"
        exit 1
    fi
fi

echo
echo "FASTQ directory: $FQDIR"
echo "Sample sheet: $SAMPLESHEET"
echo "Output directory: $OUTDIR"
echo

mkdir -p "$OUTDIR"

found_pairs=0
submitted=0

for r1 in "$FQDIR"/*.r_1.lostreads.fq.gz; do
    base=$(basename "$r1" .r_1.lostreads.fq.gz)
    r2="${r1/.r_1./.r_2.}"
    outdir_job="$OUTDIR/$base"
    mkdir -p "$outdir_job"
    if [[ -f "$r2" ]]; then
        found_pairs=$((found_pairs+1))
        jobname="demux_${base}"
        slurm_out="$outdir_job/idemux_slurm-%j.out"
        echo "-----------------------------"
        echo "Preparing job for:"
        echo "  R1: $r1"
        echo "  R2: $r2"
        echo "  Output dir: $outdir_job"
        echo "  Slurm output: $slurm_out"
        echo "  Command: sbatch --job-name=\"$jobname\" --cpus-per-task=1 --mem=200G --time=100:00:00 --partition=epyc --output=\"$slurm_out\" --wrap \"idemux --r1 $r1 --r2 $r2 --sample-sheet $SAMPLESHEET --out $outdir_job\"/fastq"
        if [[ $DRYRUN -eq 1 ]]; then
            echo "Dry run: job not submitted."
        else
            sbatch --job-name="$jobname" \
                --cpus-per-task=1 \
                --mem=200G \
                --time=100:00:00 \
                --partition=epyc \
                --output="$slurm_out" \
                --wrap "idemux --r1 $r1 --r2 $r2 --sample-sheet $SAMPLESHEET --out $outdir_job/fastq"
            submitted=$((submitted+1))
        fi
        echo
    else
        echo "WARNING: No matching R2 file for $r1"
    fi
done

echo
echo "Total pairs found: $found_pairs"
if [[ $DRYRUN -eq 1 ]]; then
    echo "Total jobs that would be submitted: $found_pairs"
else
    echo "Total jobs submitted: $submitted"
fi