# Script status and authoritative entry points

All 74 files in `scripts/` are retained because they form the final Step 17 reproducibility script set. The GitHub copies preserve the analysis logic but replace personal paths/accounts with generic placeholders. The table below prevents helpers and superseded analysis entry points from being mistaken for the final run.

| Pipeline area | Authoritative scripts | Status and use |
|---|---|---|
| Initial proteome QC | `01_proteome_busco_qc.sh`, `02_summarise_busco.sh`, `03_submit_busco_qc.sh` | Final initial BUSCO/redundancy workflow; launch with `03_submit_busco_qc.sh`. |
| Representative proteins and core set | `04_prepare_primary_proteomes.sh`, `04_select_primary_isoforms.py`, `05_primary_busco_qc.sh`, `06_summarise_primary_busco.py`, `06_summarise_primary_busco.sh`, `07_prepare_core_set.py` | Final representative-isoform and BUSCO ≥80% selection logic. |
| Combined OrthoFinder | `08_run_orthofinder.sh`, `09_submit_primary_to_orthofinder.sh` | Final combined 14-species OrthoFinder baseline. The wrapper also runs representative-proteome preparation and BUSCO. |
| *Fragum* replacement | `10_download_fragum_sueziense.sh`, `11_busco_fragum_sueziense.sh`, `12_submit_fragum_replacement.sh` | Final replacement workflow for *F. sueziense* assembly `GCA_963680895.1`. Run before active 14-species analysis when reconstructing the exact dataset. |
| Full-MFP combined tree | `13_run_iqtree_species_tree.sh` | Superseded for the final combined tree by `21_run_restricted_mfp_iqtree.slurm`; retained only for provenance. Do not use for the authoritative combined result. |
| Class-specific OrthoFinder and initial trees | `14_run_class_orthofinder_iqtree.sh`, `15_submit_class_analyses.sh` | Final class split and Bivalvia tree. The Gastropoda full-MFP tree was replaced by the restricted-MFP rerun. |
| Final restricted-MFP species trees | `21_run_restricted_mfp_iqtree.slurm` | Authoritative for `ANALYSIS=combined` and `ANALYSIS=gastropoda`. It does not implement a Bivalvia mode because the completed Bivalvia class tree was retained. |
| Step 04 | `16_step04_orthogroup_distribution.py`, `16_step04_orthogroup_distribution.slurm` | Authoritative orthogroup-distribution analysis. |
| Step 05 | `17_step05_prepare_cafe.py`, `18_step05_install_cafe5.slurm`, `19_step05_run_cafe_group.slurm`, `20_step05_validate_prepare.slurm` | `18` is an environment installer and `20` is a validation helper; `19` is the final Bivalvia/Gastropoda CAFE5 worker. |
| Step 06 | `22_step06_candidate_families.py`, `22_step06_candidate_families.slurm` | Authoritative Base-model candidate extraction. |
| Step 07 | `23_step07_functional_annotation.py`, `23_step07_functional_annotation.slurm`, `23_step07_combine.slurm` | Authoritative eggNOG-mapper 2.1.12 workflow. These are the corrected `-i`-parameter versions. |
| Step 08 | `24_step08_functional_enrichment.py`, `24_step08_functional_enrichment.slurm` | Authoritative orthogroup-level candidate-internal enrichment. |
| Step 09 | `25_step09_candidate_sequence_validation.py`, `25_step09_prepare_pfam.slurm`, `25_step09_validate_group.slurm`, `25_step09_combine.slurm` | Authoritative Pfam 38.2/HMMER validation. |
| Step 10 | `26_step10_candidate_gene_phylogenetics.py`, `26_step10_prepare.slurm`, `26_step10_gene_tree_array.slurm`, `26_step10_summarise.slurm` | Authoritative top-30-per-class gene-tree workflow. |
| Step 11 | `27_step11_reconciliation.py`, `27_step11_reconcile_group.slurm`, `27_step11_combine.slurm` | Authoritative DL-LCA reconciliation. |
| Central result copy | `28_sync_results_archive.sh`, `28_archive_results.slurm` | Operational archive step: copies results without deleting the sources. Not a biological analysis. |
| Step 12 | `29_step12_molecular_evolution.py`, all `29_step12_*.slurm`, `29_submit_step12.sh` | Authoritative HyPhy BUSTED/aBSREL/RELAX workflow. `29_submit_step12.sh` starts setup/preparation only; `29_step12_launch.slurm` validates them and launches arrays. |
| Step 13 | `30_step13_transcriptomics.py`, all `30_step13_*.slurm`, `30_submit_step13.sh` | Authoritative single-library Salmon integration. |
| Step 14 | `31_step14_integrated_prioritisation.py`, `31_step14_prioritise.slurm`, `31_step14_archive.slurm`, `31_submit_step14.sh` | Authoritative evidence-weighted candidate ranking. |
| Step 15 | `32_step15_sensitivity_bias.py`, `32_step15_sensitivity.slurm`, `32_step15_archive.slurm`, `32_submit_step15.sh` | Authoritative corrected sensitivity workflow; final run was `18431248`. |
| Step 16 | `33_step16_statistics_visualisation.R`, `33_step16_statistics_visualisation.slurm`, `33_step16_archive.slurm`, `33_submit_step16.sh` | Authoritative final tables, nine PNG figures and PDF output. |
| Step 17 | `34_step17_reproducibility_archive.py`, `34_step17_reproducibility_archive.slurm`, `34_step17_archive.slurm`, `34_submit_step17.sh` | Authoritative final inventory, checksums, validation and compact delivery archives. |

The complete final-run paths are in `metadata/FINAL_RUN_SELECTION.tsv`. Historical Slurm failures and replacements are preserved in `metadata/SLURM_JOB_LEDGER.psv`; job state alone should not override the explicit final-run selection table.
