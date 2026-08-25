#!/usr/bin/env bash
#SBATCH --job-name=step17_repro_archive
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=12:00:00
#SBATCH --output=/path/to/mollusc_project/logs/step17_repro_archive_%j.out
#SBATCH --error=/path/to/mollusc_project/logs/step17_repro_archive_%j.err

set -euo pipefail

export PROJECT=/path/to/mollusc_project
export STEP17_OUT=$PROJECT/results_or_reports/reproducibility_result_archiving/final_v1/run_${SLURM_JOB_ID}
mkdir -p "$PROJECT/results_or_reports/reproducibility_result_archiving/final_v1"
python3 "$PROJECT/scripts/34_step17_reproducibility_archive.py"
ln -sfn "run_${SLURM_JOB_ID}" "$PROJECT/results_or_reports/reproducibility_result_archiving/final_v1/latest"

