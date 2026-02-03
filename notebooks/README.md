# Notebook Launcher Scripts

Scripts for launching Jupyter and RStudio sessions on SLURM clusters.

## Setup Instructions

### First Time Setup

1. **Copy the configuration template:**
   ```bash
   cd notebooks
   cp config.sh.example config.sh
   ```

2. **Edit config.sh with your settings:**
   ```bash
   nano config.sh
   ```

3. **Configure the following variables:**
   - `PROJECT_ROOT`: Your project directory path
   - `ENV_NAME`: Your Jupyter conda environment name
   - `CONTAINER`: Path to your RStudio Singularity container
   - `REMOTE_USER`: Your cluster username
   - `JUMP_HOST`: Cluster login node hostname

### Usage

**Launch Jupyter Notebook:**
```bash
sbatch jupyter_slurm.sh
```

**Launch RStudio Server:**
```bash
sbatch rstudio_slurm.sh
```
