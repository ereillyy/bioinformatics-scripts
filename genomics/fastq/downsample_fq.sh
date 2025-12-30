#!/bin/bash
# This script down-samples the input fastq files to a specified number of reads.
# It uses seqtk to sample the reads and gzip to compress the output files.
# 
# Usage: ./downsample_fq.sh '<input_pattern>' '<output_pattern>' <num_reads> [--slurm]
#   Use ? as placeholder for 1/2 in the file patterns
#   Output pattern should NOT include .gz extension (will be added automatically)

#  run in analysis_1 conda env (seqtk version 1.4-r122)

# exit on error
set -e

# check for slurm flag
submit_slurm=false
if [[ " $* " =~ " --slurm " ]]; then
    submit_slurm=true
fi

# if --slurm flag is set, submit to SLURM and exit
if [ "$submit_slurm" = true ]; then
    # remove --slurm from arguments
    args=()
    for arg in "$@"; do
        if [ "$arg" != "--slurm" ]; then
            args+=("'$arg'")
        fi
    done
    
    echo "Submitting to SLURM..."
    sbatch --partition=epyc --mem=8G --cpus-per-task=2 --time=2:00:00 --wrap="bash $0 ${args[*]}"
    exit 0
fi

echo "Starting downsample_fq.sh"

# check if correct number of arguments provided
if [ "$#" -lt 3 ]; then
    echo "Error: Incorrect number of arguments"
    echo "Usage: $0 <input_pattern> <output_pattern> <num_reads> [--slurm]"
    echo "  Use ? as placeholder for 1/2 in the file patterns"
    echo "  Output pattern should NOT include .gz extension (will be added automatically)"
    echo "Example: $0 sample_r?.fq.gz 5mreads_sample_r?.fq 5000000"
    echo "Example with SLURM: $0 sample_r?.fq.gz 5mreads_sample_r?.fq 5000000 --slurm"
    exit 1
fi

# check if seqtk is available
if ! command -v seqtk &> /dev/null; then
    echo "Error: seqtk not found. Please ensure it is installed and in your PATH."
    exit 1
fi

# set up input and output files
input_pattern="$1"
output_pattern="$2"
reads="$3"

# remove .gz extension from output pattern if present
output_pattern="${output_pattern%.gz}"

# replace ? with 1 and 2
input_r1="${input_pattern/\?/1}"
input_r2="${input_pattern/\?/2}"
output_r1="${output_pattern/\?/1}"
output_r2="${output_pattern/\?/2}"

echo "Input R1: $input_r1"
echo "Input R2: $input_r2"
echo "Output R1: ${output_r1}.gz"
echo "Output R2: ${output_r2}.gz"
echo "Number of reads: $reads"

# check if input files exist
if [ ! -f "$input_r1" ]; then
    echo "Error: Input file not found: $input_r1"
    exit 1
fi

if [ ! -f "$input_r2" ]; then
    echo "Error: Input file not found: $input_r2"
    exit 1
fi

# check if output files already exist
if [ -f "${output_r1}.gz" ]; then
    echo "Error: Output file already exists: ${output_r1}.gz"
    echo "Please remove or rename existing output files before running."
    exit 1
fi

if [ -f "${output_r2}.gz" ]; then
    echo "Error: Output file already exists: ${output_r2}.gz"
    echo "Please remove or rename existing output files before running."
    exit 1
fi

# downsample the files
# keep same seed=100 to maintain pairing
echo -e "\nDownsampling r1... \n\t running seqtk sample -s100 $input_r1 $reads > $output_r1"
echo -e "Downsampling r2... \n\t running seqtk sample -s100 $input_r2 $reads > $output_r2"

# run seqtk commands in parallel and capture exit codes
seqtk sample -s100 "$input_r1" "$reads" > "$output_r1" &
pid1=$!
seqtk sample -s100 "$input_r2" "$reads" > "$output_r2" &
pid2=$!

# wait for both processes and check their exit codes
wait $pid1
exit1=$?
wait $pid2
exit2=$?

if [ $exit1 -ne 0 ]; then
    echo "Error: seqtk failed for R1 with exit code $exit1"
    exit 1
fi

if [ $exit2 -ne 0 ]; then
    echo "Error: seqtk failed for R2 with exit code $exit2"
    exit 1
fi

echo "Downsampling completed successfully"

# zip the output files
echo -e "\nCompressing output files... \n\t running gzip $output_r1 and $output_r2"
gzip "$output_r1" &
pid1=$!
gzip "$output_r2" &
pid2=$!

# Wait for both gzip processes and check their exit codes
wait $pid1
exit1=$?
wait $pid2
exit2=$?

if [ $exit1 -ne 0 ]; then
    echo "Error: gzip failed for R1 with exit code $exit1"
    exit 1
fi

if [ $exit2 -ne 0 ]; then
    echo "Error: gzip failed for R2 with exit code $exit2"
    exit 1
fi

echo -e "\nDone downsampling to $reads reads."