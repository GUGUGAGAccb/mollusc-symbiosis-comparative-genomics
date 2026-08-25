#!/usr/bin/env bash
#SBATCH --job-name=step14_prioritise
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --output=/path/to/mollusc_project/logs/step14_prioritise_%j.out
#SBATCH --error=/path/to/mollusc_project/logs/step14_prioritise_%j.err

set -euo pipefail
ROOT=/path/to/mollusc_project
python3 "$ROOT/scripts/31_step14_integrated_prioritisation.py" "$SLURM_JOB_ID"
