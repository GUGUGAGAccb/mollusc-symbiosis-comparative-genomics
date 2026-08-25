#!/bin/bash
#SBATCH --job-name=install_cafe5
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:30:00
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --output=/path/to/mollusc_project/logs/step05_install_cafe5_%j.out
#SBATCH --error=/path/to/mollusc_project/logs/step05_install_cafe5_%j.err

set -euo pipefail

PROJECT=/path/to/mollusc_project
ENV_PREFIX="$PROJECT/conda_envs/cafe5_env"
CONDA=/path/to/miniconda3/bin/conda
export CONDA_PKGS_DIRS="$PROJECT/conda_pkgs"

mkdir -p "$PROJECT/conda_envs" "$CONDA_PKGS_DIRS" "$PROJECT/logs"

if [[ ! -x "$ENV_PREFIX/bin/cafe5" ]]; then
    "$CONDA" create --yes --solver libmamba --prefix "$ENV_PREFIX" \
        --override-channels --channel conda-forge --channel bioconda \
        cafe=5.1.0
fi

"$ENV_PREFIX/bin/cafe5" --help >/dev/null
echo "CAFE5_READY=$ENV_PREFIX"
