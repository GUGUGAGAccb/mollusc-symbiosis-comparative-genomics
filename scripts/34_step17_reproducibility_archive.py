#!/usr/bin/env python3
"""Build the final reproducibility record and compact deliverables for Steps 1-16."""

from __future__ import annotations

import csv
import datetime as dt
import hashlib
import os
import re
import shutil
import subprocess
import tarfile
from pathlib import Path


PROJECT = Path(os.environ.get("PROJECT", "/path/to/mollusc_project"))
RESULT = PROJECT / "result"
OUT = Path(os.environ["STEP17_OUT"])
JOB_ID = os.environ.get("SLURM_JOB_ID", "unknown")


def sha256(path: Path, chunk: int = 8 * 1024 * 1024) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            block = handle.read(chunk)
            if not block:
                break
            h.update(block)
    return h.hexdigest()


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8", errors="backslashreplace") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def run_version(label: str, command: list[str], source: str) -> dict:
    try:
        proc = subprocess.run(command, text=True, capture_output=True, timeout=45, check=False)
        text = (proc.stdout + "\n" + proc.stderr).strip().replace("\t", " ")
        version = next((line.strip() for line in text.splitlines() if line.strip()), "no version text")
        status = "captured" if proc.returncode == 0 else f"command_exit_{proc.returncode}"
    except Exception as exc:
        version, status = f"unavailable: {type(exc).__name__}", "not_captured"
    return {"software_or_database": label, "version": version[:500], "source": source, "status": status}


def safe_text(value) -> str:
    """Represent filesystem names containing surrogate-escaped bytes safely in UTF-8."""
    return str(value).encode("utf-8", "backslashreplace").decode("utf-8")


def add_files_to_tar(archive: Path, root: Path, files: list[Path]) -> None:
    with tarfile.open(archive, "w:gz", compresslevel=6) as tf:
        for path in sorted(files):
            tf.add(path, arcname=safe_text(path.relative_to(root)), recursive=False)


def archive_validation(label: str, archive: Path, expected: int) -> dict:
    try:
        with tarfile.open(archive, "r:gz") as tf:
            members = [m for m in tf.getmembers() if m.isfile()]
            # Force decompression and checksum of every archived regular file.
            for member in members:
                stream = tf.extractfile(member)
                if stream is not None:
                    while stream.read(8 * 1024 * 1024):
                        pass
        ok = len(members) == expected
        return {"archive": label, "status": "PASS" if ok else "FAIL",
                "expected_files": expected, "observed_files": len(members),
                "bytes": archive.stat().st_size, "sha256": sha256(archive),
                "note": "all members decompressed successfully" if ok else "member count differs"}
    except Exception as exc:
        return {"archive": label, "status": "FAIL", "expected_files": expected,
                "observed_files": 0, "bytes": archive.stat().st_size if archive.exists() else 0,
                "sha256": sha256(archive) if archive.exists() else "", "note": repr(exc)}


def main() -> None:
    if OUT.exists():
        raise SystemExit(f"Output already exists: {OUT}")
    OUT.mkdir(parents=True)
    created = dt.datetime.now(dt.timezone.utc).isoformat()

    # Final-run selection. Historical runs remain available but are not treated as final.
    final_runs = [
        ("step01", "proteome_quality", "step01_proteome_quality/current", "final BUSCO/core-species decision set"),
        ("step02", "orthogroup_inference", "step02_orthogroup_inference", "quota-safe OrthoFinder final tables and compressed large tables"),
        ("step03", "species_tree", "step03_species_tree", "class trees and combined restricted-MFP baseline"),
        ("step04", "orthogroup_distribution", "step04_orthogroup_distribution/run_18269662", "final distribution run"),
        ("step05", "gene_family_evolution", "step05_gene_family_evolution", "Bivalvia and Gastropoda CAFE Base-model results"),
        ("step06", "candidate_gene_families", "step06_candidate_gene_families/run_18271359", "final Base-model candidate extraction"),
        ("step07", "functional_annotation", "step07_functional_annotation", "final eggNOG-mapper results"),
        ("step08", "functional_enrichment", "step08_functional_enrichment/run_18272100", "candidate-internal orthogroup enrichment"),
        ("step09", "candidate_sequence_validation", "step09_candidate_sequence_validation", "Pfam 38.2 direct validation"),
        ("step10", "candidate_family_phylogenetics", "step10_candidate_family_phylogenetics/run_18274798", "60 completed supported trees"),
        ("step11", "gene_tree_species_tree_reconciliation", "step11_gene_tree_species_tree_reconciliation", "final DL-LCA reconciliation"),
        ("step12", "molecular_evolution", "step12_molecular_evolution/run_18288232", "final HyPhy results"),
        ("step13", "transcriptomic_integration", "step13_transcriptomic_integration/run_18420304", "six single-library Salmon datasets; no DE"),
        ("step14", "integrated_candidate_prioritisation", "step14_integrated_candidate_prioritisation/run_18427340", "final evidence-weighted ranking"),
        ("step15", "sensitivity_bias_analysis", "step15_sensitivity_bias_analysis/run_18431248", "corrected baseline-consistent run; run_18428504 is historical"),
        ("step16", "statistical_analysis_visualisation", "step16_statistical_analysis_visualisation/run_18451973", "final statistics and nine-figure set"),
    ]
    final_rows = []
    for step, analysis, rel, note in final_runs:
        path = RESULT / rel
        final_rows.append({"step": step, "analysis": analysis, "selected_path": str(path),
                           "exists": "yes" if path.exists() else "no", "selection_note": note})
    write_tsv(OUT / "FINAL_RUN_SELECTION.tsv", final_rows,
              ["step", "analysis", "selected_path", "exists", "selection_note"])
    if any(r["exists"] == "no" for r in final_rows):
        raise SystemExit("One or more selected final result paths are missing")

    # Inventory all central-result files before Step 17 is copied into it.
    manifest = []
    checksum_lines = []
    for path in sorted(RESULT.rglob("*")):
        rel = path.relative_to(RESULT)
        if rel.parts and rel.parts[0].startswith("step17_"):
            continue
        stat = path.lstat()
        if path.is_symlink():
            manifest.append({"relative_path": safe_text(rel), "type": "symlink", "bytes": stat.st_size,
                             "mtime_utc": dt.datetime.fromtimestamp(stat.st_mtime, dt.timezone.utc).isoformat(),
                             "sha256": "", "link_target": safe_text(os.readlink(path))})
        elif path.is_file():
            digest = sha256(path)
            manifest.append({"relative_path": safe_text(rel), "type": "file", "bytes": stat.st_size,
                             "mtime_utc": dt.datetime.fromtimestamp(stat.st_mtime, dt.timezone.utc).isoformat(),
                             "sha256": digest, "link_target": ""})
            checksum_lines.append(f"{digest}  {safe_text(rel)}")
    write_tsv(OUT / "MASTER_RESULT_MANIFEST.tsv", manifest,
              ["relative_path", "type", "bytes", "mtime_utc", "sha256", "link_target"])
    (OUT / "CENTRAL_RESULT_SHA256SUMS").write_text("\n".join(checksum_lines) + "\n", encoding="utf-8", errors="backslashreplace")

    # Software/database record.
    env = PROJECT / "conda_envs"
    versions = [
        run_version("BUSCO", [str(env / "busco_env/bin/busco"), "--version"], str(env / "busco_env")),
        run_version("HMMER hmmscan", [str(env / "busco_env/bin/hmmscan"), "-h"], str(env / "busco_env")),
        run_version("OrthoFinder", [str(env / "orthofinder_env/bin/orthofinder"), "--version"], str(env / "orthofinder_env")),
        run_version("IQ-TREE", [str(env / "orthofinder_env/bin/iqtree"), "--version"], str(env / "orthofinder_env")),
        run_version("CAFE5", [str(env / "cafe5_env/bin/cafe5"), "-h"], str(env / "cafe5_env")),
        run_version("HyPhy", [str(env / "hyphy_env/bin/hyphy"), "--version"], str(env / "hyphy_env")),
        run_version("Python", ["python3", "--version"], "login/batch base environment"),
        run_version("R", [str(env / "busco_env/bin/Rscript"), "--version"], str(env / "busco_env")),
        run_version("ggplot2", [str(env / "busco_env/bin/Rscript"), "-e", "cat(as.character(packageVersion('ggplot2')))"], str(env / "busco_env")),
    ]
    pfam = PROJECT / "databases/pfam/current/Pfam.version"
    versions.append({"software_or_database": "Pfam", "version": pfam.read_text(encoding="utf-8").strip().replace("\n", "; ") if pfam.exists() else "missing",
                     "source": str(pfam), "status": "captured" if pfam.exists() else "missing"})
    versions.extend([
        {"software_or_database": "eggNOG-mapper", "version": "2.1.12", "source": "Step 7 result path and run metadata", "status": "recorded_from_final_results"},
        {"software_or_database": "Salmon", "version": "see Step 13 job logs/result metadata", "source": "Step 13", "status": "recorded_by_reference"},
        {"software_or_database": "MAFFT/trimAl", "version": "see Step 10 logs", "source": "Step 10", "status": "recorded_by_reference"},
    ])
    write_tsv(OUT / "SOFTWARE_AND_DATABASE_VERSIONS.tsv", versions,
              ["software_or_database", "version", "source", "status"])

    # Slurm ledger: root jobs only, including historical failures and replacements.
    ledger_cmd = ["sacct", "-X", "-u", "HPC_USER", "-S", "2026-08-03", "-P", "-n",
                  "--format=JobIDRaw,JobName%60,Account,Partition,State,ExitCode,Elapsed,AllocCPUS,ReqMem,MaxRSS,Start,End"]
    ledger = subprocess.run(ledger_cmd, text=True, capture_output=True, check=True, timeout=120).stdout
    header = "JobIDRaw|JobName|Account|Partition|State|ExitCode|Elapsed|AllocCPUS|ReqMem|MaxRSS|Start|End\n"
    (OUT / "SLURM_JOB_LEDGER.psv").write_text(header + ledger, encoding="utf-8")

    decisions = [
        ("species_scope", "Analyse Bivalvia (6) and Gastropoda (8) separately; retain 14-species combined tree as baseline."),
        ("proteome_threshold", "Core species require BUSCO complete >=80%; Fragum_fragum was replaced by Fragum_sueziense for paired genome/proteome availability and quality."),
        ("resources", "Downstream jobs capped at 12 CPUs and 16 GB memory."),
        ("iqtree", "Restricted ModelFinder set used to accelerate protein-tree inference while retaining UFBoot1000 and SH-aLRT1000."),
        ("cafe", "Relative root age 100 is exploratory; no absolute divergence-time claim. Gastropoda Gamma significance is unreliable; Base model is primary."),
        ("annotation", "eggNOG transfer and direct Pfam evidence kept distinct; annotation-source sensitivity quantified."),
        ("sequence_validation", "Exclude/review flags are retained for audit; high-risk families require manual review."),
        ("transcriptomics", "One selected library per represented species; detection only, no differential expression or cross-species TPM comparison."),
        ("prioritisation", "Missing RNA/HyPhy is not a biological zero; scores are normalised over available evidence."),
        ("sensitivity", "Corrected Step 15 run 18431248 is authoritative; historical run 18428504 is not final."),
        ("archive", "Large reconstructible OrthoFinder directories are not duplicated; final large tables use lossless tar.gz where appropriate."),
    ]
    write_tsv(OUT / "ANALYSIS_DECISIONS.tsv",
              [{"decision": k, "record": v} for k, v in decisions], ["decision", "record"])

    lineage = [
        (1, "Downloaded primary proteomes", "BUSCO + redundancy metrics", "core/supplementary species decisions", "Steps 2-3"),
        (2, "Core primary proteomes", "OrthoFinder", "orthogroups, copy-number matrices, gene trees", "Steps 3-6"),
        (3, "Single-copy orthologue alignments", "IQ-TREE restricted MFP", "supported species trees", "Steps 5,11"),
        (4, "Orthogroup tables", "distribution classification", "core/shell/species-specific summaries", "Step 16"),
        (5, "Gene-count matrices + species trees", "CAFE5 Base model", "significant family changes", "Step 6"),
        (6, "CAFE changes + sequences", "candidate-family selection", "candidate families and proteins", "Steps 7-14"),
        (7, "Candidate proteins", "eggNOG-mapper", "functional annotations", "Steps 8,14"),
        (8, "Candidate family annotations", "Fisher/BH enrichment", "candidate-internal enrichment", "Steps 14,16"),
        (9, "Candidate proteins + Pfam 38.2", "HMMER/rule validation", "pass/review/risk calls", "Steps 10,14"),
        (10, "Validated family sequences", "MAFFT + trimAl + IQ-TREE", "60 candidate gene trees", "Steps 11-12"),
        (11, "Gene trees + species trees", "DL-LCA reconciliation", "duplication/loss events", "Step 14"),
        (12, "Candidate CDS + branches", "BUSTED/aBSREL/RELAX", "selection/relaxation evidence", "Step 14"),
        (13, "Selected RNA-seq libraries", "Salmon", "candidate detection/TPM", "Step 14"),
        (14, "Steps 6-13 evidence", "evidence-weighted ranking", "ranked families and genes", "Steps 15-16"),
        (15, "Step 14 ranking", "13 sensitivity scenarios", "stability and bias reports", "Step 16"),
        (16, "Steps 1-15 summaries", "selected statistics + ggplot2", "tables, nine PNGs, PDF", "Step 17"),
        (17, "Steps 1-16 central results", "inventory/checksum/package validation", "final reproducibility archive", "delivery"),
    ]
    write_tsv(OUT / "DATA_LINEAGE.tsv", [{"step": a, "input": b, "process": c, "output": d, "downstream": e} for a,b,c,d,e in lineage],
              ["step", "input", "process", "output", "downstream"])

    large = [
        ("OrthoFinder WorkingDirectory and reconstructible sequence/alignment trees", "original orthofinder_results directories", "not duplicated to avoid inode/quota pressure"),
        ("Raw/processed RNA-seq FASTQ and BAM", "project rnaseq/source paths", "not included in compact deliverables"),
        ("Primary genomes and proteomes", "project genome/proteome directories", "referenced by provenance; not duplicated"),
        ("Pfam 38.2 database", str(PROJECT / "databases/pfam/current"), "retained on BluePebble; not included in compact package"),
        ("Conda environments", str(PROJECT / "conda_envs"), "retained on BluePebble; versions recorded instead of copying"),
    ]
    write_tsv(OUT / "EXTERNAL_AND_LARGE_DATA_INVENTORY.tsv",
              [{"data": a, "location": b, "archive_policy": c} for a,b,c in large],
              ["data", "location", "archive_policy"])

    cleanup = [{"date": "2026-08-16", "target": "/path/to/hpc_home/project/sra_cache/SRR5760179/SRR5760179.sra",
                "bytes_removed": 7883386705, "reason": "recoverable SRA download cache; paired FASTQ outputs retained",
                "recoverability": "redownload from SRA accession SRR5760179", "other_data_removed": "empty cache/tmp directories only"}]
    write_tsv(OUT / "HOME_CLEANUP_RECORD.tsv", cleanup,
              ["date", "target", "bytes_removed", "reason", "recoverability", "other_data_removed"])

    # Secret-risk scan of scripts selected for reproducibility. Never encode or search for any actual credential.
    scripts = sorted((PROJECT / "scripts").glob("*"))
    text_scripts = [p for p in scripts if p.is_file() and p.stat().st_size <= 5_000_000]
    patterns = {
        "private_key_header": re.compile(r"BEGIN (?:RSA |OPENSSH |EC )?PRIVATE KEY"),
        "password_assignment": re.compile(r"(?i)\b(?:password|passwd|pwd)\s*[:=]\s*['\"]?\S+"),
        "sshpass_usage": re.compile(r"(?i)\bsshpass\b"),
        "token_assignment": re.compile(r"(?i)\b(?:api[_-]?key|access[_-]?token|secret[_-]?key)\s*[:=]\s*['\"]?\S+"),
    }
    findings = []
    for path in text_scripts:
        try:
            text = path.read_text(encoding="utf-8", errors="replace")
        except Exception:
            continue
        for name, regex in patterns.items():
            if regex.search(text):
                findings.append({"file": str(path), "pattern": name, "status": "BLOCK_ARCHIVE"})
    if findings:
        write_tsv(OUT / "SECURITY_SCAN.tsv", findings, ["file", "pattern", "status"])
        raise SystemExit("Potential credential material detected in scripts; archive blocked")
    write_tsv(OUT / "SECURITY_SCAN.tsv",
              [{"file": "all selected project scripts", "pattern": "credential-risk patterns", "status": "PASS_no_matches"}],
              ["file", "pattern", "status"])

    # Compact deliverables requested for local export.
    tsv_files = [p for p in RESULT.rglob("*.tsv") if "step17_" not in p.relative_to(RESULT).parts]
    visual_root = RESULT / "step16_statistical_analysis_visualisation/run_18451973"
    visual_files = sorted(list(visual_root.glob("*.png")) + list(visual_root.glob("*.pdf")) +
                          [visual_root / "figure_captions.tsv", visual_root / "STEP16_REPORT.md"])
    visual_files = [p for p in visual_files if p.exists()]
    visual_archive = OUT / "mollusc_project_visualisations_step16.tar.gz"
    tsv_archive = OUT / "mollusc_project_all_tsv_tables_step01_step16.tar.gz"
    add_files_to_tar(visual_archive, RESULT, visual_files)
    add_files_to_tar(tsv_archive, RESULT, tsv_files)

    # Reproducibility source bundle: scripts plus project-level records created so far.
    scripts_archive = OUT / "reproducibility_scripts.tar.gz"
    add_files_to_tar(scripts_archive, PROJECT, text_scripts)

    validations = [
        archive_validation("visualisations", visual_archive, len(visual_files)),
        archive_validation("all_tsv_tables", tsv_archive, len(tsv_files)),
        archive_validation("reproducibility_scripts", scripts_archive, len(text_scripts)),
    ]
    write_tsv(OUT / "ARCHIVE_VALIDATION_REPORT.tsv", validations,
              ["archive", "status", "expected_files", "observed_files", "bytes", "sha256", "note"])
    if any(r["status"] != "PASS" for r in validations):
        raise SystemExit("Archive validation failed")

    readme = f"""# Mollusc comparative-genomics project: final reproducibility archive

Created: {created}
Step 17 Slurm job: {JOB_ID}

This directory closes Steps 1-16 by recording final-run selection, checksums, software/database
versions, Slurm job history, analysis decisions, data lineage, large-data policy, security scan and
validated compact deliverables. It does not alter upstream biological results.

Start here:
1. FINAL_RUN_SELECTION.tsv — authoritative final outputs for every step.
2. DATA_LINEAGE.tsv — input/process/output dependencies.
3. ANALYSIS_DECISIONS.tsv — decisions and limitations that materially affect interpretation.
4. MASTER_RESULT_MANIFEST.tsv and CENTRAL_RESULT_SHA256SUMS — integrity and provenance.
5. ARCHIVE_VALIDATION_REPORT.tsv — proof that compact archives can be fully decompressed.
6. mollusc_project_visualisations_step16.tar.gz — nine PNG figures, PDF booklet, captions/report.
7. mollusc_project_all_tsv_tables_step01_step16.tar.gz — all {len(tsv_files)} archived TSV tables.

Security: passwords, SSH private keys, access tokens and VPN credentials are intentionally excluded.
Raw sequencing data, databases, Conda environments and reconstructible OrthoFinder working files
remain in their documented BluePebble locations rather than being duplicated.
"""
    (OUT / "FINAL_PROJECT_README.md").write_text(readme, encoding="utf-8")

    methods = """# Final methods and limitations index

The authoritative per-step methods remain with each archived result. This index records the
cross-cutting limitations required for interpretation:

- Core proteomes were selected using BUSCO complete >=80%; quality-threshold sensitivity is in Step 15.
- Bivalvia and Gastropoda were analysed separately for within-class family evolution; the combined
  14-species tree is a phylogenetic baseline, not a combined CAFE test.
- CAFE uses relative root age 100 and supports exploratory family-change screening, not absolute dating.
- Gastropoda Gamma-model family significance is unreliable; Base-model results are primary.
- Functional annotations are predictions. eggNOG transfer and direct Pfam support are retained separately.
- Candidate trees support downstream reconciliation, but duplication/loss timing remains topology-dependent.
- HyPhy eligibility requires adequate codon data; missing tests are not evidence of no selection.
- RNA-seq has one selected library per represented species; no differential expression was performed.
- Step 14 scores prioritise follow-up and are not probabilities of biological truth.
- Step 15 stable-high status is intentionally stringent; sensitive families require evidence-specific review.
- Protein/family count tests in Step 16 use nested observations and are descriptive, not fully independent replication.
"""
    (OUT / "FINAL_METHODS_AND_LIMITATIONS.md").write_text(methods, encoding="utf-8")

    # Package Step 17 documents and script archive, excluding the two larger user deliverables.
    bundle_candidates = [p for p in OUT.iterdir() if p.is_file() and p.name not in {
        visual_archive.name, tsv_archive.name, "PROJECT_ARCHIVE_COMPLETE", "COMPLETE",
        "STEP17_SHA256SUMS", "result_manifest.tsv", "mollusc_project_reproducibility_step17.tar.gz"}]
    repro_bundle = OUT / "mollusc_project_reproducibility_step17.tar.gz"
    add_files_to_tar(repro_bundle, OUT, bundle_candidates)
    validations.append(archive_validation("step17_reproducibility_bundle", repro_bundle, len(bundle_candidates)))
    write_tsv(OUT / "ARCHIVE_VALIDATION_REPORT.tsv", validations,
              ["archive", "status", "expected_files", "observed_files", "bytes", "sha256", "note"])
    if validations[-1]["status"] != "PASS":
        raise SystemExit("Reproducibility bundle validation failed")

    # Final local-file integrity record.
    artifact_files = sorted(p for p in OUT.iterdir() if p.is_file() and p.name not in {"STEP17_SHA256SUMS", "result_manifest.tsv", "PROJECT_ARCHIVE_COMPLETE", "COMPLETE"})
    (OUT / "STEP17_SHA256SUMS").write_text("\n".join(f"{sha256(p)}  {p.name}" for p in artifact_files) + "\n", encoding="utf-8")
    (OUT / "PROJECT_ARCHIVE_COMPLETE").write_text(f"Step 17 archive complete: job {JOB_ID}\n", encoding="utf-8")
    (OUT / "COMPLETE").write_text("Step 17 reproducibility and result archiving complete\n", encoding="utf-8")

    all_files = sorted(p for p in OUT.iterdir() if p.is_file())
    result_rows = [{"file": p.name, "bytes": p.stat().st_size,
                    "sha256": sha256(p) if p.name != "result_manifest.tsv" else "self_excluded"} for p in all_files]
    write_tsv(OUT / "result_manifest.tsv", result_rows, ["file", "bytes", "sha256"])
    print(f"Step 17 completed: {OUT}")


if __name__ == "__main__":
    main()
