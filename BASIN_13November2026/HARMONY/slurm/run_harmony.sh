#!/usr/bin/env bash
#SBATCH --job-name=harmony
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# HARMONY SLURM submission script
# --------------------------------
# Submits the HARMONY normalization pipeline (Pipeline 2) to a SLURM cluster.
#
# Usage:
#   From the HARMONY/ directory, run:
#     sbatch slurm/run_harmony.sh templates/config.yaml
#
#   Or with a custom config:
#     sbatch slurm/run_harmony.sh /path/to/your/config.yaml
#
# Resource guidance:
#   Time:   8h is usually sufficient for a fresh run from BAM files.
#           4h is typically enough when reusing baseline bigWigs from a prior run
#           (set reuse_baseline_from_run in the config).
#   CPUs:   4 is a reasonable default. bamCoverage will use SLURM_CPUS_PER_TASK
#           automatically. Increase to 8 for large sample sets.
#   Memory: 32G is sufficient for most runs. Reduce to 16G for reuse runs.

set -euo pipefail
mkdir -p logs
cd "${SLURM_SUBMIT_DIR}"

# -----------------------------------------------------------------------
# Load required modules
# Adjust module names and versions to match your cluster's environment.
# The versions listed here were used for development and validation:
#   python/3.12.4
#   samtools/1.23
#   deeptools/3.5.6
#
# On systems using Lmod, module names are typically identical or similar.
# On systems using Environment Modules, syntax is the same but names may
# differ (e.g., "Python/3.12.4" vs "python/3.12.4").
# -----------------------------------------------------------------------

# Ensure Lmod is available in non-interactive shells
if [ -f /etc/profile.d/lmod.sh ]; then
  source /etc/profile.d/lmod.sh
elif [ -f /usr/share/lmod/lmod/init/bash ]; then
  source /usr/share/lmod/lmod/init/bash
fi

module purge || true
module load python/3.12.4   # adjust version for your cluster
module load samtools/1.23   # adjust version for your cluster
module load deeptools/3.5.6 # adjust version for your cluster

# -----------------------------------------------------------------------
# Set up a per-job temporary directory and clean it up on exit
# -----------------------------------------------------------------------
export TMPDIR=/tmp/$USER/${SLURM_JOB_ID}
mkdir -p "$TMPDIR"
trap 'rm -rf "$TMPDIR"' EXIT

# -----------------------------------------------------------------------
# Run HARMONY
# Defaults to templates/config.yaml if no argument is provided.
# -----------------------------------------------------------------------
CONFIG=${1:-templates/config.yaml}
python3 harmony.py --config "$CONFIG"
