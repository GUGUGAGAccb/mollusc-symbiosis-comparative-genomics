#!/usr/bin/env bash
#SBATCH --job-name=step13_summarise
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --output=/path/to/mollusc_project/logs/step13_summarise_%j.out
#SBATCH --error=/path/to/mollusc_project/logs/step13_summarise_%j.err

set -euo pipefail
PREP_JOB_ID=${1:?Usage: sbatch script PREP_JOB_ID}
ROOT=/path/to/mollusc_project
python3 "$ROOT/scripts/30_step13_transcriptomics.py" summarise "$PREP_JOB_ID"
