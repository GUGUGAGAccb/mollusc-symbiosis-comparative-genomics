# Portability, data and security notes

## Sanitised-snapshot policy

The files in `scripts/` preserve the analysis logic, parameters and workflow relationships of the Step 17 archive. Before GitHub export, user-specific paths, the cluster username and the Slurm account were replaced by explicit placeholders. The unmodified source archive is retained offline and is identified by the SHA-256 value in `metadata/SOURCE_ARCHIVE_SHA256.txt`; it is not included in this repository.

Before rerunning outside the original project, review at least the following values throughout `scripts/`:

| Setting | Published placeholder | Required action |
|---|---|---|
| Project root | `/path/to/mollusc_project` | Replace with the new absolute project root. |
| Conda base | `/path/to/miniconda3` | Replace with the local Conda/Mamba installation or module. |
| eggNOG taxonomy database | `EGGNOG_TAXDB` or `/path/to/eggnog.taxa.db` | Point to the installed `eggnog.taxa.db` before Step 09. |
| Slurm account | `YOUR_SLURM_ACCOUNT` | Replace `#SBATCH --account` with the authorised allocation. |
| Partitions | `compute` and `test` | Map to local partitions/queues. |
| Modules | eggNOG-mapper, HMMER, trimAl, Salmon | Replace with site-specific module or environment names. |
| Species/run paths | several fixed `Results_Aug04` and final-run paths | Update only if reconstructing with new output names. |

The downstream resource policy was at most 12 CPUs and 16 GB per job. The early original BUSCO/combined-OrthoFinder baseline scripts predate that policy and retain their executed requests of 8 CPUs/32 GB or 20 CPUs/80 GB. Class-specific and downstream jobs conform to the later limit. When reproducing on another cluster, benchmark memory before reducing the early baseline jobs.

## Inputs intentionally excluded from Git

The repository does not contain:

- the 14 primary protein FASTA files or their genomes/GFF annotations;
- raw or cleaned RNA-seq FASTQ files and BAM/CRAM files;
- BUSCO lineage downloads;
- Pfam 38.2 HMM files and indexes;
- the Gene Ontology OBO file;
- Conda environments or software binaries;
- OrthoFinder `WorkingDirectory`, `Orthologues`, sequence/alignment collections and other large reconstructible outputs;
- final result tables and figures, which were delivered separately from the script repository.

Their original locations and archive policy are recorded in `metadata/EXTERNAL_AND_LARGE_DATA_INVENTORY.tsv`. Final result provenance and checksums are recorded in `metadata/MASTER_RESULT_MANIFEST.tsv` and `metadata/CENTRAL_RESULT_SHA256SUMS`.

## Required input naming

Active proteomes must be non-empty files named `Species_genus_species.faa` under `orthofinder_proteomes/`. The exact 14 species names are embedded in `14_run_class_orthofinder_iqtree.sh`. Do not retain both *Fragum fragum* and *Fragum sueziense* in the active input directory: archive the superseded *F. fragum* input outside the active 14-proteome set.

Step 12 additionally expects genome/GFF or transcript-derived CDS inputs in the paths defined near the top of `29_step12_molecular_evolution.py`. Step 13 expects six paired cleaned RNA-seq libraries listed in `30_step13_transcriptomics.py`; it is a detection-only integration and not a replicated differential-expression design.

## Credentials and privacy

No password, VPN credential, SSH private key, access token or GitHub token belongs in this repository. Never place credentials in:

- a script or Slurm command line;
- `.env` files committed to Git;
- job logs;
- README examples;
- Git remote URLs.

The published files contain no original BluePebble username, personal project/home path, personal email address or original Slurm account. The repository remains private as an additional precaution. Authentication must be handled by SSH agents, system keychains or GitHub's device/browser login, not by saving passwords in files.

## Recommended pre-push checks

```bash
git status --short
rg -l -i 'password|passwd|credential|secret|token|private[ _-]?key' . \
  --glob '!metadata/SECURITY_SCAN.tsv' \
  --glob '!docs/PORTABILITY_AND_SECURITY.md'
find . -name '*.swp' -o -name '.DS_Store'
```

Any match must be reviewed in context. Documentation that states credentials are excluded is expected; actual secret values are not.
