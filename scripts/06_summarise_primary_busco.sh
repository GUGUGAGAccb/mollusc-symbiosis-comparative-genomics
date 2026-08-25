#!/bin/bash
#SBATCH --job-name=primary_busco_summary
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --mem=4G
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --output=/path/to/mollusc_project/logs/primary_busco_summary_%j.out
#SBATCH --error=/path/to/mollusc_project/logs/primary_busco_summary_%j.err

set -euo pipefail

SCRATCH=/path/to/mollusc_project
"${SCRATCH}/conda_envs/orthofinder_env/bin/python3" \
  "${SCRATCH}/scripts/06_summarise_primary_busco.py" \
  --proteome-dir "${SCRATCH}/orthofinder_primary_proteomes" \
  --busco-dir "${SCRATCH}/busco_primary_proteomes" \
  --output "${SCRATCH}/primary_proteome_busco_summary.tsv"

column -t -s $'\t' "${SCRATCH}/primary_proteome_busco_summary.tsv"
