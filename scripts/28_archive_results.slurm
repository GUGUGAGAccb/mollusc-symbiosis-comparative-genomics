#!/usr/bin/env bash
#SBATCH --job-name=archive_results
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=1-00:00:00
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --output=/path/to/mollusc_project/logs/archive_results_%j.out
#SBATCH --error=/path/to/mollusc_project/logs/archive_results_%j.err

set -euo pipefail

ROOT=/path/to/mollusc_project
bash "$ROOT/scripts/28_sync_results_archive.sh" all
