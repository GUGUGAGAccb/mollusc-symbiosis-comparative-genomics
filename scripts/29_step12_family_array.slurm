#!/usr/bin/env bash
#SBATCH --job-name=step12_hyphy
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --output=/path/to/mollusc_project/logs/step12_hyphy_%x_%A_%a.out
#SBATCH --error=/path/to/mollusc_project/logs/step12_hyphy_%x_%A_%a.err

set -euo pipefail

GROUP=${1:?Usage: sbatch --array=0-29 script GROUP PREP_JOB_ID}
REQUESTED_GROUP=$GROUP
PREP_JOB_ID=${2:?Usage: sbatch --array=0-29 script GROUP PREP_JOB_ID}
case "$GROUP" in bivalvia|gastropoda) ;; *) echo "Unsupported group: $GROUP" >&2; exit 2 ;; esac

ROOT=/path/to/mollusc_project
RUN="$ROOT/results_or_reports/molecular_evolution/hyphy_busted_absrel_relax_v1/run_${PREP_JOB_ID}"
MANIFEST="$RUN/manifest_${GROUP}.tsv"
PYTHON="$ROOT/conda_envs/orthofinder_env/bin/python"
HYPHY="$ROOT/conda_envs/hyphy_env/bin/hyphy"

test -s "$RUN/PREPARED"
test -s "$MANIFEST"
test -s "$ROOT/conda_envs/hyphy_env/STEP12_READY"
test -x "$HYPHY"

eval "$("$PYTHON" - "$MANIFEST" "$SLURM_ARRAY_TASK_ID" <<'PY'
import csv, shlex, sys
path, wanted = sys.argv[1], int(sys.argv[2])
with open(path, newline='') as handle:
    rows = list(csv.DictReader(handle, delimiter='\t'))
matches = [row for row in rows if int(row['array_index']) == wanted]
if len(matches) != 1:
    raise SystemExit(f'Expected one manifest row for array index {wanted}, found {len(matches)}')
row = matches[0]
keys = ['group','orthogroup','status','reason','codon_alignment','clean_tree','relax_tree','output_dir','run_relax']
for key in keys:
    print(f"{key.upper()}={shlex.quote(row[key])}")
PY
)"

test "$GROUP" = "$REQUESTED_GROUP"
mkdir -p "$OUTPUT_DIR"
if [[ "$STATUS" != eligible ]]; then
  printf '%s\n' "$REASON" > "$OUTPUT_DIR/SKIPPED"
  printf 'Skipped %s %s: %s\n' "$GROUP" "$ORTHOGROUP" "$REASON"
  exit 0
fi

test -s "$CODON_ALIGNMENT"
test -s "$CLEAN_TREE"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export HYPHY_THREADS="${SLURM_CPUS_PER_TASK}"

successes=0
failures=0
run_method() {
  local method=$1
  shift
  if timeout 21600 "$HYPHY" "$method" "$@" > "$OUTPUT_DIR/${method}.stdout" 2> "$OUTPUT_DIR/${method}.stderr"; then
    successes=$((successes + 1))
  else
    failures=$((failures + 1))
    printf '%s\n' "$method" >> "$OUTPUT_DIR/FAILED_METHODS"
  fi
}

rm -f "$OUTPUT_DIR/COMPLETE" "$OUTPUT_DIR/FAILED" "$OUTPUT_DIR/FAILED_METHODS"
run_method busted \
  --alignment "$CODON_ALIGNMENT" --tree "$CLEAN_TREE" --branches All \
  --code Universal --output "$OUTPUT_DIR/busted.json"
run_method absrel \
  --alignment "$CODON_ALIGNMENT" --tree "$CLEAN_TREE" --branches All \
  --code Universal --output "$OUTPUT_DIR/absrel.json"

if [[ "$RUN_RELAX" == yes ]]; then
  test -s "$RELAX_TREE"
  run_method relax \
    --alignment "$CODON_ALIGNMENT" --tree "$RELAX_TREE" \
    --test test --reference reference --code Universal --output "$OUTPUT_DIR/relax.json"
else
  printf 'RELAX not run because foreground/reference branch counts were insufficient.\n' > "$OUTPUT_DIR/RELAX_SKIPPED"
fi

if (( successes > 0 )); then
  printf 'Step 12 HyPhy family analysis complete; successes=%s failures=%s\n' "$successes" "$failures" > "$OUTPUT_DIR/COMPLETE"
else
  printf 'All requested HyPhy methods failed\n' > "$OUTPUT_DIR/FAILED"
fi
printf 'Finished %s %s: successes=%s failures=%s\n' "$GROUP" "$ORTHOGROUP" "$successes" "$failures"
