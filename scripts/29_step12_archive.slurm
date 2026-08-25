#!/usr/bin/env bash
#SBATCH --job-name=step12_archive
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=02:00:00
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --output=/path/to/mollusc_project/logs/step12_archive_%j.out
#SBATCH --error=/path/to/mollusc_project/logs/step12_archive_%j.err

set -euo pipefail

ROOT=/path/to/mollusc_project
BASE="$ROOT/results_or_reports/molecular_evolution/hyphy_busted_absrel_relax_v1"
SOURCE=$(readlink -f "$BASE/latest")
test -s "$SOURCE/COMPLETE"
mkdir -p "$ROOT/result/reproducibility/scripts"
rsync -a "$ROOT/scripts/" "$ROOT/result/reproducibility/scripts/"
"$ROOT/scripts/28_sync_results_archive.sh" one step12_molecular_evolution "$SOURCE" "$(basename "$SOURCE")"
printf 'Step 12 archived from %s\n' "$SOURCE"
