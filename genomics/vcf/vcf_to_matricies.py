import os
import argparse
import shutil
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
import sigProfilerPlotting as plot

parser = argparse.ArgumentParser(
    description='Generate mutation signature matrices and plots from VCF files',
    epilog='''
Example usage:
  python vcf_to_matricies.py --indir /path/to/vcfs --outdir /path/to/output --genome GRCh38 --percentage

This script processes VCF files to generate mutation signature matrices using SigProfilerMatrixGenerator,
then creates plots for all available mutation types (SBS, DBS, ID).
''',
    formatter_class=argparse.RawDescriptionHelpFormatter
)
parser.add_argument('--indir', '-i', required=True, help='Input directory containing VCF files.')
parser.add_argument('--outdir', '-o', required=True, help='Output directory for results.')
parser.add_argument('--genome', '-g', default='GRCh38', help='Reference genome (e.g., GRCh38, GRCh37, mm10).')
parser.add_argument('--percentage', action='store_true', help='Plot as percentage of total mutation count (default: False).')
args = parser.parse_args()

# Define mutation types for each category
SBS_TYPES = ["6", "24", "96", "384", "1536"]
DBS_TYPES = ["78", "186"]
ID_TYPES = ["83", "415"]

# Configuration for each mutation category
MUTATION_CONFIGS = [
    {"category": "SBS", "types": SBS_TYPES, "plot_func": plot.plotSBS},
    {"category": "DBS", "types": DBS_TYPES, "plot_func": plot.plotDBS},
    {"category": "ID", "types": ID_TYPES, "plot_func": plot.plotID}
]

# Prepare directories
indir = os.path.abspath(args.indir)
outdir = os.path.abspath(args.outdir)
os.makedirs(os.path.join(outdir, "sigproplots"), exist_ok=True)

project_name = 'sigpro'

print(f"[INFO] Starting mutation signature analysis")
print(f"[INFO] Input directory: {indir}")
print(f"[INFO] Output directory: {outdir}")
print(f"[INFO] Reference genome: {args.genome}")
print(f"[INFO] Percentage mode: {args.percentage}")

# Generate matrices using SigProfilerMatrixGenerator
print(f"\n[INFO] Running SigProfilerMatrixGenerator...")
try:
    matGen.SigProfilerMatrixGeneratorFunc(
        project=project_name,
        reference_genome=args.genome,
        path_to_input_files=indir
    )
    print(f"[SUCCESS] Matrix generation complete")
except Exception as e:
    print(f"[ERROR] Matrix generation failed: {e}")
    raise

# Clean up temporary directories and move matrices to output
print(f"\n[INFO] Cleaning up temporary files...")
cleanup_dirs = [
    os.path.join(indir, "input"),
    os.path.join(indir, "logs"),
    os.path.join(indir, "output", "vcf_files")
]

for cleanup_dir in cleanup_dirs:
    if os.path.exists(cleanup_dir):
        print(f"[INFO] Removing {cleanup_dir}")
        shutil.rmtree(cleanup_dir)
    else:
        print(f"[INFO] Directory not found, skipping: {cleanup_dir}")

# Move matrices from indir/output to outdir
matrices_source = os.path.join(indir, "output")
matrices_dest = os.path.join(outdir, "matrices")

if os.path.exists(matrices_source):
    print(f"[INFO] Moving matrices from {matrices_source} to {matrices_dest}")
    shutil.move(matrices_source, matrices_dest)
    print(f"[SUCCESS] Matrices moved to {matrices_dest}")
else:
    print(f"[WARNING] Matrices source directory not found: {matrices_source}")

# Process each mutation category and type
print(f"\n[INFO] Starting plotting for all mutation types...")
total_plots_generated = 0
total_plots_skipped = 0

for config in MUTATION_CONFIGS:
    category = config["category"]
    types = config["types"]
    plot_func = config["plot_func"]
    
    print(f"\n[INFO] Processing {category} mutations...")
    
    for mut_type in types:
        matrix_file = os.path.join(matrices_dest, category, f"{project_name}.{category}{mut_type}.all")
        
        if os.path.exists(matrix_file):
            print(f"[INFO] Found matrix: {category}{mut_type}")
            print(f"[INFO] Generating plot for {category}{mut_type}...")
            
            try:
                plot_func(
                    matrix_path=matrix_file,
                    output_path=os.path.join(outdir, "sigproplots", category),
                    project=project_name,
                    plot_type=mut_type,
                    percentage=args.percentage
                )
                total_plots_generated += 1
                print(f"[SUCCESS] {category}{mut_type} plot generated")
            except Exception as e:
                print(f"[WARNING] Failed to generate {category}{mut_type} plot: {e}")
                total_plots_skipped += 1
        else:
            print(f"[WARNING] Matrix not found: {matrix_file} (skipping)")
            total_plots_skipped += 1

print(f"\n[INFO] Analysis complete!")
print(f"[INFO] Total plots generated: {total_plots_generated}")
print(f"[INFO] Total plots skipped: {total_plots_skipped}")
print(f"[INFO] Results saved to: {outdir}")
