#!/usr/bin/env bash
#SBATCH --job-name=step09_pfam_setup
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --output=/path/to/mollusc_project/logs/step09_pfam_setup_%j.out
#SBATCH --error=/path/to/mollusc_project/logs/step09_pfam_setup_%j.err

set -euo pipefail

ROOT=/path/to/mollusc_project
BASE="$ROOT/databases/pfam"
RUN="$BASE/run_${SLURM_JOB_ID}"
URL=https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release
mkdir -p "$RUN" "$ROOT/logs"

module purge
module load apps/hmmer/3.4

for file in Pfam-A.hmm.gz Pfam-A.hmm.dat.gz Pfam-A.clans.tsv.gz Pfam.version.gz md5_checksums; do
  curl --fail --location --retry 4 --continue-at - --output "$RUN/$file" "$URL/$file"
done

cd "$RUN"
md5sum --check --ignore-missing md5_checksums
gzip -dc Pfam-A.hmm.gz > Pfam-A.hmm
gzip -dc Pfam-A.hmm.dat.gz > Pfam-A.hmm.dat
gzip -dc Pfam-A.clans.tsv.gz > Pfam-A.clans.tsv
gzip -dc Pfam.version.gz > Pfam.version
hmmpress -f Pfam-A.hmm
for suffix in h3f h3i h3m h3p; do test -s "Pfam-A.hmm.$suffix"; done
printf '%s\n' "$URL/Pfam-A.hmm.gz" "$URL/Pfam-A.hmm.dat.gz" "$URL/Pfam-A.clans.tsv.gz" "$URL/Pfam.version.gz" > SOURCE_URLS.txt
sha256sum Pfam-A.hmm Pfam-A.hmm.dat Pfam-A.clans.tsv Pfam.version > SHA256SUMS
printf 'complete\n' > READY
ln -sfn "run_${SLURM_JOB_ID}" "$BASE/current"
printf 'Pfam setup complete: %s\n' "$RUN"
