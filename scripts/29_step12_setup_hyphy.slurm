#!/usr/bin/env bash
#SBATCH --job-name=step12_hyphy_setup
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --output=/path/to/mollusc_project/logs/step12_hyphy_setup_%j.out
#SBATCH --error=/path/to/mollusc_project/logs/step12_hyphy_setup_%j.err

set -euo pipefail

ROOT=/path/to/mollusc_project
ENV_PREFIX="$ROOT/conda_envs/hyphy_env"
READY="$ENV_PREFIX/STEP12_READY"
CONDA_EXE=/path/to/miniconda3/bin/conda
TEMP_PKGS="$ROOT/tmp/step12_conda_pkgs_${SLURM_JOB_ID}"
mkdir -p "$ROOT/logs" "$ROOT/tmp"

if [[ -x "$ENV_PREFIX/bin/hyphy" ]]; then
  "$ENV_PREFIX/bin/hyphy" --version
  printf 'HyPhy environment already available\n' > "$READY"
  exit 0
fi

if [[ -d "$ENV_PREFIX" ]]; then
  rm -rf -- "$ENV_PREFIX"
fi
mkdir -p "$TEMP_PKGS"
cleanup() {
  rm -rf -- "$TEMP_PKGS"
}
trap cleanup EXIT

export CONDA_PKGS_DIRS="$TEMP_PKGS"
"$CONDA_EXE" create -y -p "$ENV_PREFIX" -c conda-forge -c bioconda hyphy
"$ENV_PREFIX/bin/hyphy" --version
"$CONDA_EXE" list -p "$ENV_PREFIX" --explicit > "$ENV_PREFIX/step12_conda_explicit.txt"
"$ENV_PREFIX/bin/hyphy" busted --help > "$ENV_PREFIX/busted_help.txt" 2>&1
"$ENV_PREFIX/bin/hyphy" absrel --help > "$ENV_PREFIX/absrel_help.txt" 2>&1
"$ENV_PREFIX/bin/hyphy" relax --help > "$ENV_PREFIX/relax_help.txt" 2>&1
printf 'HyPhy installation complete for Step 12\n' > "$READY"

