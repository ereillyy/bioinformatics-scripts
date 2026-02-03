#!/usr/bin/env bash
#
# USAGE: sbatch jupyter_slurm.sh
#
#SBATCH --job-name=jupyter
#SBATCH --partition=epyc
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --time=48:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

# Load configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [ -f "$SCRIPT_DIR/config.sh" ]; then
    source "$SCRIPT_DIR/config.sh"
else
    echo "Error: config.sh not found. Please copy config.sh.example to config.sh and configure it."
    exit 1
fi

# Create logs directory in PROJECT_ROOT
mkdir -p "$PROJECT_ROOT/logs"
cd "$PROJECT_ROOT"

# Load Conda
source ~/miniconda3/etc/profile.d/conda.sh
conda activate "$ENV_NAME"

echo "=== SLURM ==="
echo "Job ID:        $SLURM_JOB_ID"
echo "Node:          $(hostname)"
echo "Start time:    $(date)"
echo "Work dir:      $PROJECT_ROOT"
echo "Conda env:     $ENV_NAME"
which python
python --version

PORT=8888
JUPYTER_LOG="$PROJECT_ROOT/logs/jupyter-${SLURM_JOB_ID}.log"

echo "Starting JupyterLab ..."
jupyter lab \
  --no-browser \
  --port=$PORT \
  --ip=0.0.0.0 \
  --port-retries=100 \
  > "$JUPYTER_LOG" 2>&1 &

# Wait for Jupyter to start and write URL to log
echo "Waiting for Jupyter to start..."
sleep 10

# Try to extract URL from log file
URL=""
if [ -f "$JUPYTER_LOG" ]; then
  URL=$(grep -m1 -o "http://.*" "$JUPYTER_LOG" 2>/dev/null || echo "")
fi

echo "=============================================================="
echo "JUPYTER IS RUNNING ON (remote): $(hostname):$PORT"
echo "From your LAPTOP, run this SSH tunnel:"
echo "ssh -J ${REMOTE_USER}@${JUMP_HOST} ${REMOTE_USER}@$(hostname) -L $PORT:localhost:$PORT"
echo ""
if [ -n "$URL" ]; then
  echo "Then open in your browser:"
  echo "   $URL"
else
  echo "URL not found yet. Check the log file for the URL:"
  echo "   $JUPYTER_LOG"
  echo "Or wait a moment and run: grep 'http://' $JUPYTER_LOG"
fi
echo "=============================================================="

wait
