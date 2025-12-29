#!/usr/bin/env bash
#
# USAGE: sbatch rstudio_slurm.sh
#
#SBATCH --job-name=rstudio
#SBATCH --partition=epyc
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --time=48:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

PROJECT_ROOT="/mnt/scratchc/sjlab/ereilly"
CONTAINER="/mnt/scratchc/sjlab/ereilly/containers/bioconductor_3.18.sif"

mkdir -p logs
cd "$PROJECT_ROOT"

echo "=== SLURM ==="
echo "Job ID:        $SLURM_JOB_ID"
echo "Node:          $(hostname)"
echo "Start time:    $(date)"
echo "Work dir:      $PROJECT_ROOT"
echo "Container:     $CONTAINER"

# Check if container exists
if [ ! -f "$CONTAINER" ]; then
    echo "ERROR: Container not found at $CONTAINER"
    exit 1
fi

# Pick a port
PORT=8788

# Create temporary directory for RStudio Server
export RSTUDIO_TMP="${PROJECT_ROOT}/tmp/rstudio-${SLURM_JOB_ID}"
mkdir -p "${RSTUDIO_TMP}/var/lib"
mkdir -p "${RSTUDIO_TMP}/var/run"
mkdir -p "${RSTUDIO_TMP}/tmp"

# Create database config to avoid warnings
cat > "${RSTUDIO_TMP}/database.conf" <<EOF
provider=sqlite
directory=${RSTUDIO_TMP}/var/lib
EOF

echo "Starting RStudio Server..."

# Get the hostname and port
HOSTNAME=$(hostname)

echo "=============================================================="
echo "RSTUDIO SERVER IS RUNNING ON (remote): ${HOSTNAME}:${PORT}"
echo ""
echo "From your LAPTOP, run this SSH tunnel:"
echo "ssh -J reilly01@clust1-sub-1.cri.camres.org reilly01@${HOSTNAME} -L ${PORT}:localhost:${PORT}"
echo ""
echo "Then open in your browser:"
echo "   http://localhost:${PORT}"
echo ""
echo "Note: Login may not be required (auth-none enabled)"
echo "Working directory: ${PROJECT_ROOT}"
echo "=============================================================="
echo ""
echo "To stop the server:"
echo "scancel ${SLURM_JOB_ID}"
echo "=============================================================="

# Start RStudio Server in container (run in foreground to keep job alive)
singularity exec \
    --bind ${PROJECT_ROOT}:${PROJECT_ROOT} \
    --bind ${RSTUDIO_TMP}/var/lib:/var/lib/rstudio-server \
    --bind ${RSTUDIO_TMP}/var/run:/var/run/rstudio-server \
    --bind ${RSTUDIO_TMP}/tmp:/tmp \
    --bind ${RSTUDIO_TMP}/database.conf:/etc/rstudio/database.conf \
    ${CONTAINER} \
    bash -c "rserver \
        --server-daemonize 0 \
        --www-port=${PORT} \
        --www-address=0.0.0.0 \
        --auth-none=1 \
        --server-user=$(whoami) \
        --server-data-dir=${RSTUDIO_TMP}/var/run"