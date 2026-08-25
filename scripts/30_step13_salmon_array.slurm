#!/usr/bin/env bash
#SBATCH --job-name=step13_salmon
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=12G
#SBATCH --time=12:00:00
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --output=/path/to/mollusc_project/logs/step13_salmon_%A_%a.out
#SBATCH --error=/path/to/mollusc_project/logs/step13_salmon_%A_%a.err

set -euo pipefail
PREP_JOB_ID=${1:?Usage: sbatch script PREP_JOB_ID}
ROOT=/path/to/mollusc_project
RUN="$ROOT/results_or_reports/transcriptomic_integration/single_library_salmon_v1/run_${PREP_JOB_ID}"
MANIFEST="$RUN/sample_manifest.tsv"

test -s "$RUN/PREPARED"
test -s "$MANIFEST"
module load salmon/1.10.2-gve3

row=$(sed -n "$((SLURM_ARRAY_TASK_ID + 2))p" "$MANIFEST")
test -n "$row"
IFS=$'\t' read -r array_index group species accession read1 read2 reference candidate reference_n candidate_n read1_bytes read2_bytes output_dir scope <<< "$row"
test "$array_index" -eq "$SLURM_ARRAY_TASK_ID"
for input in "$read1" "$read2" "$reference" "$candidate"; do test -s "$input"; done

mkdir -p "$output_dir"
if [ -s "$output_dir/COMPLETE" ]; then
  printf 'Already complete: %s %s\n' "$species" "$accession"
  exit 0
fi

index="$output_dir/salmon_index"
quant="$output_dir/quant"
if [ ! -s "$index/versionInfo.json" ]; then
  rm -rf "$index"
  salmon index --transcripts "$reference" --index "$index" --kmerLen 31 \
    --keepDuplicates --threads "$SLURM_CPUS_PER_TASK" \
    > "$output_dir/index.stdout.log" 2> "$output_dir/index.stderr.log"
fi

rm -rf "$quant"
salmon quant --index "$index" --libType A --mates1 "$read1" --mates2 "$read2" \
  --validateMappings --seqBias --gcBias --threads "$SLURM_CPUS_PER_TASK" \
  --output "$quant" > "$output_dir/quant.stdout.log" 2> "$output_dir/quant.stderr.log"

test -s "$quant/quant.sf"
test -s "$quant/aux_info/meta_info.json"
printf 'species\t%s\nrun_accession\t%s\ncompleted_job\t%s_%s\n' \
  "$species" "$accession" "$SLURM_ARRAY_JOB_ID" "$SLURM_ARRAY_TASK_ID" > "$output_dir/COMPLETE"
printf 'Completed Salmon quantification: %s %s\n' "$species" "$accession"
