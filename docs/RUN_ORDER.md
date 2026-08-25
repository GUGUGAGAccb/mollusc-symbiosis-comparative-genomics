# Final BluePebble Analysis Workflow: Complete Run Order, Purpose, and Execution

This document explains how to reconstruct Steps 01–17 according to the final project logic. The files in `scripts/` retain the analysis logic, parameters, and dependencies from the final Step 17 archive, but all personal project paths, home paths, cluster usernames, and Slurm accounts were replaced with generic placeholders before publication. Read `PORTABILITY_AND_SECURITY.md` and replace every configuration placeholder before the first run.

## 1. Pre-run preparation

### 1.1 Set up the project structure

The sanitised GitHub version uses the following project-root placeholder throughout:

```text
/path/to/mollusc_project
```

Before running the workflow, replace `/path/to/mollusc_project`, `/path/to/miniconda3`, `/path/to/hpc_home`, `/path/to/eggnog.taxa.db`, `HPC_USER`, and `YOUR_SLURM_ACCOUNT` with values authorised on the target cluster. The eggNOG taxonomy database can alternatively be specified through the `EGGNOG_TAXDB` environment variable. Do not recommit a customised version containing personal paths or credentials to a public repository.

Create at least the following directories:

```text
mollusc_project/
├── scripts/
├── logs/
├── orthofinder_proteomes/
├── genome_fastas/
├── busco_downloads/lineages/metazoa_odb10/
├── conda_envs/
├── databases/gene_ontology/
├── results_or_reports/
└── result/
```

Copy the repository scripts into the project and restore executable permissions:

```bash
cd /path/to/mollusc_project
cp /path/to/repository/scripts/* scripts/
chmod u+x scripts/*
mkdir -p logs results_or_reports result
```

### 1.2 Input data

The active `orthofinder_proteomes/` directory must contain exactly 14 non-empty `.faa` files. Their filenames must match the species identifiers embedded in the scripts. The final class-level design is:

- Bivalvia (6): `Acanthocardia_echinata`, `Americardia_media`, `Cerastoderma_edule`, `Fragum_sueziense`, `Tridacna_gigas`, and `Tridacna_maxima`.
- Gastropoda (8): `Aplysia_californica`, `Berghia_stephanieae`, `Bullacta_exarata`, `Elysia_chlorotica`, `Elysia_crispata`, `Elysia_marginata`, `Onchidella_celtica`, and `Plakobranchus_ocellatus`.

Do not retain both `Fragum_fragum.faa` and `Fragum_sueziense.faa` in the active input directory. The former can be moved to a separate historical-input directory, but it must not enter the final 14-species analysis.

### 1.3 Software and databases

The final workflow uses BUSCO, seqkit, OrthoFinder, DIAMOND, FAMSA, FastTree, IQ-TREE, CAFE5, eggNOG-mapper, HMMER, the Gene Ontology OBO file, MAFFT, trimAl, HyPhy, Salmon, Python, R, and ggplot2. Explicitly recorded versions include IQ-TREE 3.1.3, HyPhy 2.5.101, eggNOG-mapper 2.1.12, Pfam 38.2, Salmon 1.10.2, Python 3.13.11, R 4.5.2, and ggplot2 4.0.3. See `metadata/SOFTWARE_AND_DATABASE_VERSIONS.tsv` for the complete record.

### 1.4 Resource policy

Class-level analyses and all jobs from Step 04 onward request no more than 12 CPUs and 16 GB per job. The earliest proteome BUSCO and combined 14-species OrthoFinder scripts predate that limit and retain their executed requests: 8 CPUs/32 GB for BUSCO and 20 CPUs/80 GB for combined OrthoFinder. If every rerun must fit within 12 CPUs/16 GB, benchmark memory on a small target-cluster test before reducing those early jobs; do not assume OrthoFinder can be downscaled without consequences.

## 2. Ordered execution of Steps 01–17

### Step 01 — Proteome Quality Assessment

**Purpose:** Validate FASTA structure, count proteins, identify duplicate identifiers and invalid characters, evaluate BUSCO completeness and single-copy/duplicated/fragmented/missing proportions, select representative transcripts, and define the core species set at BUSCO complete ≥80%.

**Main software:** seqkit, BUSCO with `metazoa_odb10`, and Python.

**1A. Reconstruct the final *Fragum* replacement data, only if *F. sueziense* has not already been prepared:**

```bash
cd /path/to/mollusc_project
bash scripts/12_submit_fragum_replacement.sh
```

This chain runs:

1. `10_download_fragum_sueziense.sh`: downloads and verifies the genome, peptide, and GFF3 files for `GCA_963680895.1` from Ensembl Rapid Release; decompresses them into the appropriate directories; and selects representative proteins.
2. `11_busco_fragum_sueziense.sh`: runs BUSCO on the representative proteome and refreshes the summary.

**1B. Initial proteome quality control:**

```bash
bash scripts/03_submit_busco_qc.sh
```

The wrapper submits the `01_proteome_busco_qc.sh` BUSCO array and, after successful completion, runs `02_summarise_busco.sh`.

**1C. Representative proteins and representative-proteome BUSCO:**

The authoritative wrapper connecting Steps 01 and 02 is:

```bash
bash scripts/09_submit_primary_to_orthofinder.sh
```

Its first three dependent jobs run `04_prepare_primary_proteomes.sh` → the `05_primary_busco_qc.sh` array → `06_summarise_primary_busco.sh`. The final job proceeds into the combined OrthoFinder analysis in Step 02.

**Main outputs:**

- `proteome_busco_summary.tsv`
- `primary_proteome_selection.tsv`
- `primary_proteome_busco_summary.tsv`
- `orthofinder_species_decisions.tsv`
- `orthofinder_primary_proteomes/*.faa`
- `orthofinder_core_complete80/*.faa`

---

### Step 02 — Orthogroup Inference

**Purpose:** Infer orthogroups, orthologous and paralogous relationships, gene-copy-number matrices, single-copy orthologues, and initial OrthoFinder species trees in the combined 14-species baseline and the two class-level datasets.

**Main software:** OrthoFinder, DIAMOND, FAMSA, and FastTree.

**2A. Combined 14-species baseline:**

The previous `09_submit_primary_to_orthofinder.sh` wrapper automatically submits `08_run_orthofinder.sh`. If Step 01 was completed separately, submit it directly:

```bash
sbatch scripts/08_run_orthofinder.sh
```

**2B. Separate Bivalvia and Gastropoda analyses:**

```bash
bash scripts/15_submit_class_analyses.sh
```

For each class, `14_run_class_orthofinder_iqtree.sh` builds a symbolic-link input directory, runs OrthoFinder, and continues to an initial full-ModelFinder IQ-TREE analysis. The two class jobs can run in parallel.

**Main outputs:**

- `Orthogroups/Orthogroups.tsv`
- `Orthogroups/Orthogroups.GeneCount.tsv`
- `Orthogroups/Orthogroups_SingleCopyOrthologues.txt`
- `Species_Tree/SpeciesTree_rooted.txt`
- OrthoFinder gene trees, duplication tables, and hierarchical orthogroups

---

### Step 03 — Species Tree Reconstruction

**Purpose:** Construct strongly supported species trees from concatenated single-copy-orthologue protein alignments, providing the phylogenetic framework for CAFE5 and gene-tree/species-tree reconciliation.

**Main software:** IQ-TREE 3.1.3 with 1,000 ultrafast bootstrap replicates and 1,000 SH-aLRT replicates.

The final choices are:

- Bivalvia: retain the class-level MFP tree completed by `14_run_class_orthofinder_iqtree.sh`.
- Gastropoda: use the restricted-MFP tree.
- Combined 14 species: use the restricted-MFP baseline tree.
- The combined full-MFP entry point in `13_run_iqtree_species_tree.sh` was superseded and is not part of the authoritative final run order.

Submit the two final restricted-MFP analyses:

```bash
combined_iq=$(sbatch --parsable \
  --job-name=iqtree_restricted_combined \
  --export=ALL,ANALYSIS=combined \
  scripts/21_run_restricted_mfp_iqtree.slurm)

gastropoda_iq=$(sbatch --parsable \
  --job-name=iqtree_restricted_gastropoda \
  --export=ALL,ANALYSIS=gastropoda \
  scripts/21_run_restricted_mfp_iqtree.slurm)

printf 'combined=%s\ngastropoda=%s\n' "$combined_iq" "$gastropoda_iq"
```

The restricted ModelFinder candidate set is `Q.INSECT,Q.YEAST,JTT,LG,WAG`, with frequency modes `FU,F`, rate modes `E,I,G,I+G,R`, at most six rate categories, and random seed 12345.

**Main outputs:** `SpeciesTree_IQTREE.treefile`, `.contree`, `.iqtree`, `.log`, and the renamed concatenated alignment.

---

### Step 04 — Orthogroup Distribution Analysis

**Purpose:** Quantify core, soft-core, shell, lineage-specific, and species-specific orthogroups; measure protein assignment and unassigned-protein rates; and describe orthogroup composition within each class.

**Main tools:** Python and the custom distribution-analysis script.

```bash
step04=$(sbatch --parsable scripts/16_step04_orthogroup_distribution.slurm)
echo "$step04"
```

This step depends only on Step 02 and can run in parallel with the later part of Step 03 and preparatory work for Step 05. The authoritative historical result for manuscript tables and figures is `run_18269662`.

**Main output directory:** `results_or_reports/orthogroup_distribution/run_<JOB_ID>`.

---

### Step 05 — Gene-Family Expansion and Contraction Analysis

**Purpose:** Convert the class-level gene-count matrices and species trees into CAFE5 inputs, estimate the family-evolution rate λ, and identify significant branch-level expansions and contractions.

**Main software:** CAFE5. The relative root age is fixed at 100 for exploratory relative-time screening only.

Install the environment before first use:

```bash
cafe_setup=$(sbatch --parsable scripts/18_step05_install_cafe5.slurm)
```

Optionally validate Bivalvia input preparation:

```bash
cafe_check=$(sbatch --parsable --dependency=afterok:"$cafe_setup" \
  scripts/20_step05_validate_prepare.slurm)
```

Submit the two classes in parallel after their respective species trees are available:

```bash
biv_cafe=$(sbatch --parsable \
  --dependency=afterok:"$cafe_setup" \
  --job-name=step05_cafe_bivalvia \
  --export=ALL,GROUP=bivalvia \
  scripts/19_step05_run_cafe_group.slurm)

gas_cafe=$(sbatch --parsable \
  --dependency=afterok:"$cafe_setup":"$gastropoda_iq" \
  --job-name=step05_cafe_gastropoda \
  --export=ALL,GROUP=gastropoda \
  scripts/19_step05_run_cafe_group.slurm)
```

Primary interpretation must use the Base/error-model output. The Gastropoda Gamma model completed, but its diagnostics reported failure rates above 20% for every input family, so its significance results are not suitable as primary evidence.

**Main outputs:** `cafe_input_primary.tsv`, the ultrametric relative-time tree, `cafe_base_error/`, `cafe_gamma_k2/`, `cafe_input_summary.tsv`, and the output manifest.

---

### Step 06 — Candidate Gene-Family Identification

**Purpose:** Extract significant expanded, contracted, and terminal-branch-change families from the CAFE5 Base model, then combine change magnitude and family-level significance to generate candidate sets and protein FASTA files.

**Main tools:** Python.

```bash
step06=$(sbatch --parsable \
  --dependency=afterok:"$biv_cafe":"$gas_cafe" \
  scripts/22_step06_candidate_families.slurm)
```

Key thresholds are family FDR 0.05, high-confidence branch p-value 0.01, candidate branch p-value 0.05, and a high-confidence minimum absolute copy-number change of 2.

**Main output directory:** `results_or_reports/candidate_gene_families/base_model_v1/run_<JOB_ID>`.

---

### Step 07 — Functional Annotation

**Purpose:** Assign eggNOG seed orthologues, functional descriptions, preferred names, GO terms, KEGG terms, EC numbers, COG categories, CAZy terms, and transferred Pfam annotations to candidate proteins, then create orthogroup-level summaries.

**Main software:** eggNOG-mapper 2.1.12 and DIAMOND.

Run both classes in parallel, then combine them:

```bash
biv_ann=$(sbatch --parsable --dependency=afterok:"$step06" \
  --job-name=step07_eggnog_bivalvia \
  scripts/23_step07_functional_annotation.slurm bivalvia)

gas_ann=$(sbatch --parsable --dependency=afterok:"$step06" \
  --job-name=step07_eggnog_gastropoda \
  scripts/23_step07_functional_annotation.slurm gastropoda)

step07_merge=$(sbatch --parsable \
  --dependency=afterok:"$biv_ann":"$gas_ann" \
  scripts/23_step07_combine.slurm)
```

**Main outputs:** `protein_functional_annotations.tsv`, `candidate_family_functional_annotations.tsv`, long-format GO/KEGG/Pfam tables, annotation-coverage summaries, and the combined summary.

---

### Step 08 — Functional Enrichment Analysis

**Purpose:** Perform Fisher exact tests across candidate subsets and apply Benjamini–Hochberg FDR correction for GO BP/MF/CC, KEGG, Pfam, CAZy, and EC terms.

**Main tools:** Python, Fisher exact tests, and Benjamini–Hochberg FDR correction.

```bash
step08=$(sbatch --parsable \
  --dependency=afterok:"$step07_merge" \
  scripts/24_step08_functional_enrichment.slurm)
```

The background comprises CAFE-significant candidate orthogroups with the relevant annotation, not the entire proteome. Results therefore represent relative enrichment within the candidate set and must not be interpreted as whole-proteome enrichment.

**Main output directory:** `results_or_reports/functional_enrichment/orthogroup_level_v1/run_<JOB_ID>`.

---

### Step 09 — Candidate Sequence Validation

**Purpose:** Validate candidate sequences using direct Pfam 38.2 HMM hits, protein length, sequence completeness, exact duplication, taxonomic source, and high-risk functional flags; classify sequences as `pass`, `review`, `review_high`, or `exclude_recommended`.

**Main software:** HMMER 3.4, Pfam 38.2, and Python.

Download, verify, and run `hmmpress` on Pfam first:

```bash
pfam_job=$(sbatch --parsable scripts/25_step09_prepare_pfam.slurm)
```

Validate the two classes in parallel and then combine them:

```bash
biv_val=$(sbatch --parsable --dependency=afterok:"$pfam_job":"$biv_ann" \
  --job-name=step09_validate_bivalvia \
  scripts/25_step09_validate_group.slurm bivalvia)

gas_val=$(sbatch --parsable --dependency=afterok:"$pfam_job":"$gas_ann" \
  --job-name=step09_validate_gastropoda \
  scripts/25_step09_validate_group.slurm gastropoda)

step09_merge=$(sbatch --parsable \
  --dependency=afterok:"$biv_val":"$gas_val" \
  scripts/25_step09_combine.slurm)
```

**Main outputs:** `sequence_validation.tsv`, `family_validation_summary.tsv`, `direct_pfam_domains.tsv`, `exact_duplicate_groups.tsv`, `high_risk_sequence_flags.tsv`, `validated_pass.faa`, `manual_review.faa`, and `exclude_recommended.faa`.

---

### Step 10 — Candidate Gene-Family Phylogenetics

**Purpose:** Select up to 30 validated candidate families per class, align them with MAFFT, trim them with trimAl, and reconstruct IQ-TREE gene trees to examine within-family duplication and lineage-specific diversification.

**Main software:** MAFFT, trimAl, and IQ-TREE.

```bash
step10_prep=$(sbatch --parsable \
  --dependency=afterok:"$step09_merge" \
  scripts/26_step10_prepare.slurm)

biv_trees=$(sbatch --parsable \
  --dependency=afterok:"$step10_prep" \
  --array=0-29%3 \
  --job-name=step10_tree_bivalvia \
  scripts/26_step10_gene_tree_array.slurm bivalvia)

gas_trees=$(sbatch --parsable \
  --dependency=afterok:"$step10_prep" \
  --array=0-29%3 \
  --job-name=step10_tree_gastropoda \
  scripts/26_step10_gene_tree_array.slurm gastropoda)

step10_merge=$(sbatch --parsable \
  --dependency=afterany:"$biv_trees":"$gas_trees" \
  scripts/26_step10_summarise.slurm)
```

Each array task requests 4 CPUs/8 GB, and each array permits at most three concurrent tasks. Supported trees use 1,000 UFBoot and 1,000 SH-aLRT replicates. If support calculation fails, the scripts allow an explicitly labelled ML-only soft fallback; the final summary must report these cases.

**Main outputs:** per-family `trimmed.faa`, `.treefile`, `.iqtree`, and `.contree` files, noting that ML-only fallbacks may lack `.contree`, plus `gene_tree_summary.tsv` and `failed_or_missing_trees.tsv`.

---

### Step 11 — Gene-Tree/Species-Tree Reconciliation

**Purpose:** Reconcile candidate gene trees against the corresponding class-level species trees under a duplication–loss LCA model, localising duplication and loss events without inferring event timing directly from gene-tree topology alone.

**Main tools:** Custom Python DL-LCA reconciliation with duplication cost = 1 and loss cost = 1.

```bash
biv_rec=$(sbatch --parsable --dependency=afterok:"$step10_merge" \
  --job-name=step11_reconcile_bivalvia \
  scripts/27_step11_reconcile_group.slurm bivalvia)

gas_rec=$(sbatch --parsable --dependency=afterok:"$step10_merge" \
  --job-name=step11_reconcile_gastropoda \
  scripts/27_step11_reconcile_group.slurm gastropoda)

step11_merge=$(sbatch --parsable \
  --dependency=afterok:"$biv_rec":"$gas_rec" \
  scripts/27_step11_combine.slurm)
```

**Main outputs:** `family_reconciliation_summary.tsv`, event tables, per-family reconciled trees, and the combined summary.

After Step 11, copy the authoritative Step 01–11 results into the central `result/` directory. This operation copies rather than moves data and does not delete source results:

```bash
archive_01_11=$(sbatch --parsable \
  --dependency=afterok:"$step11_merge" \
  scripts/28_archive_results.slurm)
```

---

### Step 12 — Molecular Evolution Analysis

**Purpose:** Map candidate proteins back to CDS, create codon-aware alignments, and use HyPhy BUSTED, aBSREL, and RELAX to test gene-wide positive selection, branch-specific selection, and changes in selection intensity.

**Main software:** HyPhy 2.5.101, together with the embedded CDS/protein mapping and codon-alignment procedures.

`29_submit_step12.sh` submits only the environment check and input preparation. The recommended explicit dependency chain is:

```bash
hyphy_setup=$(sbatch --parsable scripts/29_step12_setup_hyphy.slurm)
step12_prep=$(sbatch --parsable \
  --dependency=afterok:"$step11_merge" \
  scripts/29_step12_prepare.slurm)

step12_launch=$(sbatch --parsable \
  --dependency=afterok:"$hyphy_setup":"$step12_prep" \
  scripts/29_step12_launch.slurm "$step12_prep")
```

`29_step12_launch.slurm` validates the HyPhy command-line interface and both 30-row class manifests, then automatically submits the two arrays, the summary job, and the archive job.

**Main output directory:** `results_or_reports/molecular_evolution/hyphy_busted_absrel_relax_v1/run_<PREP_JOB_ID>`. The principal outputs are the family-level molecular-evolution table and method-specific summaries.

---

### Step 13 — Optional Transcriptomic Integration

**Purpose:** Use one public paired-end library per represented species to provide candidate-transcript detection and within-library TPM support. This design does not support differential-expression testing or direct comparison of absolute TPM values across species.

**Main software:** Salmon 1.10.2 and Python.

This step depends on candidate CDS generated in Step 12 and requires the six cleaned read pairs listed in `30_step13_transcriptomics.py`:

```bash
bash scripts/30_submit_step13.sh
```

The wrapper automatically submits prepare → Salmon array with indices 0–5 and maximum concurrency 2 → summarise → archive.

The final six run accessions are `ERR15985856`, `ERR15696446`, `SRR34876314`, `ERR12342469`, `ERR10378018`, and `SRR28022337`. The *Fragum fragum* library `ERR12245519` must not be used as an expression proxy for the core species *Fragum sueziense*.

**Main outputs:** `candidate_expression.tsv`, `species_expression_summary.tsv`, `orthogroup_species_expression.tsv`, `orthogroup_transcriptomic_integration.tsv`, and an explicit unavailable-species table.

---

### Step 14 — Integrated Candidate Prioritisation

**Purpose:** Integrate gene-family change, functional annotation, sequence quality, gene-tree support, reconciliation, molecular-evolution evidence, and optional expression evidence to rank candidate families and genes.

**Main tools:** Python and a custom evidence-weighted scoring system.

After confirming that the Step 06–13 `latest/COMPLETE` markers point to the intended results, run:

```bash
bash scripts/31_submit_step14.sh
```

The maximum family-score contributions are family change 25, functional annotation 10, sequence quality 20, gene tree 10, reconciliation 10, molecular evolution 20, and expression support 5. When HyPhy or RNA evidence is unavailable, scores are normalised against the maximum available score; missing evidence must not be interpreted as a zero biological effect. Step 08 enrichment supplies context but is not counted twice in the score.

**Main outputs:** family ranking, gene ranking, evidence matrix, priority tiers, methods and limitations, and a source manifest.

---

### Step 15 — Sensitivity and Bias Analysis

**Purpose:** Test ranking robustness under 13 alternative scenarios involving proteome quality, duplication, annotation source, missing evidence, and weighting choices; quantify candidate-rank stability and potential bias.

**Main tools:** Python.

Do not use historical run `18428504`. The final authoritative logic corresponds to the corrected script and run `18431248`. For a new rerun, pass the new Step 14 directory explicitly:

```bash
STEP14_RUN=/path/to/mollusc_project/results_or_reports/integrated_candidate_prioritisation/evidence_weighted_v1/latest
bash scripts/32_submit_step15.sh "$STEP14_RUN"
```

**Main outputs:** scenario rankings, rank correlations, top-candidate stability, quality-sensitivity and annotation-source-sensitivity tables, and the summary.

---

### Step 16 — Statistical Analysis and Visualisation

**Purpose:** Summarise proteome quality, orthogroup composition, family-size changes, phylogenetic results, functional enrichment, sequence validation, molecular evolution, expression support, integrated rankings, and sensitivity analyses; generate manuscript-ready tables and figures.

**Main software:** R 4.5.2 and ggplot2 4.0.3.

```bash
STEP14_RUN=/path/to/mollusc_project/results_or_reports/integrated_candidate_prioritisation/evidence_weighted_v1/latest
STEP15_RUN=/path/to/mollusc_project/results_or_reports/sensitivity_bias_analysis/ranking_sensitivity_v1/latest
bash scripts/33_submit_step16.sh "$STEP14_RUN" "$STEP15_RUN"
```

The wrapper submits the statistics/visualisation job followed by the central archive job.

**Main outputs:** statistical TSV files, nine PNG figures, a PDF figure collection, figure captions, the integrated report, and a manifest. The final authoritative historical run is `18451973`.

---

### Step 17 — Reproducibility and Result Archiving

**Purpose:** Freeze the final Step 01–16 run selection; record software and database versions, the Slurm ledger, data lineage, analysis decisions, the result manifest, SHA-256 checksums, and sensitive-information scans; validate extractability; and create compact delivery packages.

**Main tools:** Python, tar/gzip, SHA-256, and the central result-synchronisation script.

Confirm that the Step 12–16 archive jobs have completed and that the central `result/` directory contains Steps 01–16, then run:

```bash
bash scripts/34_submit_step17.sh
```

The wrapper submits `34_step17_reproducibility_archive.slurm` and, after successful completion, runs `34_step17_archive.slurm` to copy Step 17 into the central result directory.

**Main outputs:**

- `FINAL_RUN_SELECTION.tsv`
- `DATA_LINEAGE.tsv`
- `ANALYSIS_DECISIONS.tsv`
- `SOFTWARE_AND_DATABASE_VERSIONS.tsv`
- `SLURM_JOB_LEDGER.psv`
- `MASTER_RESULT_MANIFEST.tsv`
- `CENTRAL_RESULT_SHA256SUMS`
- `SECURITY_SCAN.tsv`
- `ARCHIVE_VALIDATION_REPORT.tsv`
- Compact archives containing reproducibility scripts, all TSV files, and visualisations

## 3. Dependencies and safe parallel execution

The strict primary dependency chain is:

```text
Step01 → Step02 → Step03 → Step05 → Step06 → Step07 → Step09 → Step10 → Step11
                                                └────→ Step08
Step10 + Step11 → Step12 → Step13
Step06–Step13 → Step14 → Step15 → Step16 → Step17
```

The following components can safely run in parallel:

- Bivalvia and Gastropoda class jobs in Step 02.
- Combined and Gastropoda restricted-MFP jobs in Step 03.
- Step 04 after Step 02, in parallel with later species-tree analysis.
- The two class jobs in Step 05.
- The two class jobs in Step 07.
- The two class jobs in Step 09.
- The two 30-task arrays in Step 10, each limited to three concurrent tasks.
- The two class jobs in Step 11.
- The two 30-task HyPhy arrays in Step 12, each limited to three concurrent tasks.

Critical ordering constraints are: Step 05 must use the final class-level trees; Step 10 must wait for Step 09; Step 11 must wait for Step 10; Step 13 requires candidate CDS from Step 12; Step 14 must wait for the required evidence from Steps 06–13; and Step 17 must wait for every result archive to complete.

## 4. Job monitoring and completion criteria

```bash
squeue -u HPC_USER
sacct -j JOB_ID --format=JobIDRaw,JobName,State,ExitCode,Elapsed,AllocCPUS,ReqMem,MaxRSS
```

Do not treat Slurm `COMPLETED` alone as evidence that a bioinformatics step succeeded. For every step, also verify that:

1. The corresponding `COMPLETE` or `PREPARED` marker exists.
2. Principal TSV, tree, or figure files are non-empty.
3. The `.err` file contains no traceback, out-of-memory event, quota or no-space error, missing-input error, or tool-level fatal error.
4. The `latest` symbolic link points to the authoritative run.
5. The original result directory still exists after archiving, confirming that the operation copied rather than moved data.

Use `metadata/FINAL_RUN_SELECTION.tsv`, rather than directory dates, to identify authoritative historical result paths.
