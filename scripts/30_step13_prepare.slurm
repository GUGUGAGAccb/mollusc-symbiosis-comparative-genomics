#!/usr/bin/env bash
#SBATCH --job-name=step13_prepare
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --output=/path/to/mollusc_project/logs/step13_prepare_%j.out
#SBATCH --error=/path/to/mollusc_project/logs/step13_prepare_%j.err

set -euo pipefail
ROOT=/path/to/mollusc_project
python3 "$ROOT/scripts/30_step13_transcriptomics.py" prepare "$SLURM_JOB_ID"
