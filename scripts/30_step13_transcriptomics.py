#!/usr/bin/env python3
"""Prepare and summarise Step 13 single-library transcriptomic integration."""

from __future__ import annotations

import csv
import json
import os
import sys
from collections import defaultdict
from pathlib import Path
from statistics import mean


ROOT = Path("/path/to/mollusc_project")
BASE = ROOT / "results_or_reports/transcriptomic_integration/single_library_salmon_v1"
STEP12 = ROOT / "results_or_reports/molecular_evolution/hyphy_busted_absrel_relax_v1/run_18288232"

SAMPLES = [
    ("bivalvia", "Acanthocardia_echinata", "ERR15985856"),
    ("bivalvia", "Americardia_media", "ERR15696446"),
    ("gastropoda", "Bullacta_exarata", "SRR34876314"),
    ("gastropoda", "Elysia_crispata", "ERR12342469"),
    ("bivalvia", "Tridacna_gigas", "ERR10378018"),
    ("bivalvia", "Tridacna_maxima", "SRR28022337"),
]


def fasta_ids(path: Path) -> list[str]:
    identifiers: list[str] = []
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                identifiers.append(line[1:].split(None, 1)[0].strip())
    return identifiers


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def replace_symlink(link: Path, target: Path) -> None:
    if link.is_symlink():
        link.unlink()
    elif link.exists():
        raise RuntimeError(f"Refusing to replace non-symlink path: {link}")
    link.symlink_to(target.name)


def prepare(run_id: str) -> None:
    run = BASE / f"run_{run_id}"
    run.mkdir(parents=True, exist_ok=True)
    (run / "salmon_results").mkdir(exist_ok=True)

    manifest_rows: list[dict] = []
    available = {species: (group, accession) for group, species, accession in SAMPLES}
    core_rows = read_tsv(STEP12 / "cds_source_summary.tsv")
    audit_rows: list[dict] = []

    for row in core_rows:
        species = row["species"]
        if species in available:
            group, accession = available[species]
            status = "available_single_library"
            reason = "one paired-end RNA-seq library; expression support only"
        else:
            group = "bivalvia" if species in {
                "Acanthocardia_echinata", "Americardia_media", "Cerastoderma_edule",
                "Fragum_sueziense", "Tridacna_gigas", "Tridacna_maxima"
            } else "gastropoda"
            accession = ""
            status = "unavailable"
            reason = "no matching cleaned RNA-seq library in project"
            if species == "Fragum_sueziense":
                reason = "Fragum_fragum library exists but is not the core species Fragum_sueziense"
        audit_rows.append({
            "group": group,
            "species": species,
            "availability_status": status,
            "run_accession": accession,
            "reason": reason,
        })

    # Record the non-core congener explicitly without treating it as expression evidence.
    audit_rows.append({
        "group": "supplementary_not_integrated",
        "species": "Fragum_fragum",
        "availability_status": "available_taxon_mismatch",
        "run_accession": "ERR12245519",
        "reason": "not substituted for Fragum_sueziense",
    })

    for index, (group, species, accession) in enumerate(SAMPLES):
        read_dir = ROOT / "rnaseq_to_protein" / species / "clean_reads"
        r1 = read_dir / f"{accession}_1.clean.fastq.gz"
        r2 = read_dir / f"{accession}_2.clean.fastq.gz"
        reference = ROOT / "predicted_proteomes" / f"{species}.transcriptome_cds.fna"
        candidate = STEP12 / "cds_by_species" / f"{species}.candidate_cds.fna"
        required = [r1, r2, reference, candidate]
        missing = [str(path) for path in required if not path.is_file() or path.stat().st_size == 0]
        if missing:
            raise FileNotFoundError("Missing required Step 13 inputs: " + ", ".join(missing))

        reference_ids = fasta_ids(reference)
        candidate_ids = fasta_ids(candidate)
        duplicate_reference_ids = len(reference_ids) - len(set(reference_ids))
        absent = sorted(set(candidate_ids) - set(reference_ids))
        if duplicate_reference_ids:
            raise ValueError(f"{species}: {duplicate_reference_ids} duplicated reference identifiers")
        if absent:
            raise ValueError(f"{species}: {len(absent)} candidate IDs absent from complete CDS reference")

        output_dir = run / "salmon_results" / species
        manifest_rows.append({
            "array_index": index,
            "group": group,
            "species": species,
            "run_accession": accession,
            "read1": r1,
            "read2": r2,
            "reference_cds": reference,
            "candidate_cds": candidate,
            "reference_transcripts": len(reference_ids),
            "candidate_transcripts": len(candidate_ids),
            "read1_bytes": r1.stat().st_size,
            "read2_bytes": r2.stat().st_size,
            "output_dir": output_dir,
            "analysis_scope": "single_library_expression_support",
        })

    write_tsv(
        run / "sample_manifest.tsv",
        [
            "array_index", "group", "species", "run_accession", "read1", "read2",
            "reference_cds", "candidate_cds", "reference_transcripts", "candidate_transcripts",
            "read1_bytes", "read2_bytes", "output_dir", "analysis_scope",
        ],
        manifest_rows,
    )
    write_tsv(
        run / "data_availability.tsv",
        ["group", "species", "availability_status", "run_accession", "reason"],
        audit_rows,
    )

    methods = """Step 13: optional transcriptomic integration

Scope
-----
Six core species have one matching paired-end public RNA-seq library each. Salmon 1.10.2
is used to quantify cleaned reads against each species' complete transcript-derived CDS
reference. Candidate transcript estimates are extracted by exact sequence identifier and
summarised by species and orthogroup.

Interpretation limits
---------------------
1. Each species has only one selected library, with no biological replication or balanced
   condition design. Differential expression testing is therefore not performed.
2. TPM is interpreted within a library as transcript-detection/expression-support evidence;
   absolute TPM is not compared directly among species.
3. For these six species, proteomes were predicted from transcriptome assemblies that include
   the selected RNA-seq data. Expression support is consequently partly circular and must not
   be treated as independent validation of gene existence.
4. Fragum_fragum ERR12245519 is not used as a proxy for the core species Fragum_sueziense.
5. Missing species are reported explicitly and are not assigned zero expression.

Expression-support bins
-----------------------
undetected: TPM = 0; trace: 0 < TPM < 0.1; low: 0.1 <= TPM < 1;
moderate: 1 <= TPM < 10; high: TPM >= 10.
"""
    (run / "METHODS_AND_LIMITATIONS.txt").write_text(methods)
    (run / "PREPARED").write_text(f"prepared_run\t{run_id}\nsamples\t{len(manifest_rows)}\n")
    replace_symlink(BASE / "latest_prepared", run)
    print(f"Prepared Step 13 run: {run}")


def as_float(value: str | None) -> float | None:
    if value is None or value == "":
        return None
    try:
        return float(value)
    except ValueError:
        return None


def support_bin(tpm: float) -> str:
    if tpm == 0:
        return "undetected"
    if tpm < 0.1:
        return "trace"
    if tpm < 1:
        return "low"
    if tpm < 10:
        return "moderate"
    return "high"


def summarise(run_id: str) -> None:
    run = BASE / f"run_{run_id}"
    if not (run / "PREPARED").is_file():
        raise FileNotFoundError(f"Prepared run not found: {run}")
    manifest = read_tsv(run / "sample_manifest.tsv")

    validation: dict[tuple[str, str], dict[str, str]] = {}
    for group in ("bivalvia", "gastropoda"):
        path = ROOT / f"results_or_reports/candidate_sequence_validation/{group}/pfam_hmmer_v1/latest/sequence_validation.tsv"
        for row in read_tsv(path):
            validation[(row["species"], row["gene_id"])] = row

    family_evolution: dict[tuple[str, str], dict[str, str]] = {}
    for row in read_tsv(STEP12 / "family_molecular_evolution_results.tsv"):
        family_evolution[(row["group"], row["orthogroup"])] = row

    expression_rows: list[dict] = []
    species_rows: list[dict] = []

    for sample in manifest:
        quant = Path(sample["output_dir"]) / "quant/quant.sf"
        complete = Path(sample["output_dir"]) / "COMPLETE"
        if not quant.is_file() or quant.stat().st_size == 0 or not complete.is_file():
            raise FileNotFoundError(f"Incomplete Salmon result for {sample['species']}: {quant}")
        quant_rows = read_tsv(quant)
        estimates = {row["Name"]: row for row in quant_rows}
        candidate_ids = fasta_ids(Path(sample["candidate_cds"]))

        detected_01 = 0
        detected_1 = 0
        detected_10 = 0
        matched = 0
        for gene_id in candidate_ids:
            estimate = estimates.get(gene_id)
            if estimate:
                matched += 1
                length = float(estimate["Length"])
                effective_length = float(estimate["EffectiveLength"])
                tpm = float(estimate["TPM"])
                reads = float(estimate["NumReads"])
                quant_status = "quantified"
            else:
                length = effective_length = tpm = reads = 0.0
                quant_status = "missing_from_quantification"
            detected_01 += tpm >= 0.1
            detected_1 += tpm >= 1
            detected_10 += tpm >= 10
            val = validation.get((sample["species"], gene_id), {})
            orthogroup = val.get("orthogroup", "")
            evo = family_evolution.get((sample["group"], orthogroup), {})
            expression_rows.append({
                "group": sample["group"],
                "species": sample["species"],
                "run_accession": sample["run_accession"],
                "gene_id": gene_id,
                "orthogroup": orthogroup,
                "validation_status": val.get("validation_status", "not_in_step09_table"),
                "validation_flags": val.get("validation_flags", ""),
                "length": f"{length:.0f}",
                "effective_length": f"{effective_length:.3f}",
                "tpm": f"{tpm:.6f}",
                "num_reads": f"{reads:.6f}",
                "expression_support": support_bin(tpm),
                "quantification_status": quant_status,
                "best_event_direction": evo.get("best_event_direction", ""),
                "best_event_change": evo.get("best_event_change", ""),
                "selection_score": evo.get("selection_score", ""),
                "busted_q_bh": evo.get("busted_q_bh", ""),
                "relax_q_bh": evo.get("relax_q_bh", ""),
                "relax_interpretation": evo.get("relax_interpretation", ""),
            })

        all_detected = sum(float(row["TPM"]) >= 0.1 for row in quant_rows)
        species_rows.append({
            "group": sample["group"],
            "species": sample["species"],
            "run_accession": sample["run_accession"],
            "reference_transcripts": len(quant_rows),
            "reference_detected_tpm_ge_0.1": all_detected,
            "candidate_transcripts": len(candidate_ids),
            "candidate_matched_in_quant": matched,
            "candidate_detected_tpm_ge_0.1": detected_01,
            "candidate_detected_tpm_ge_1": detected_1,
            "candidate_detected_tpm_ge_10": detected_10,
            "candidate_detection_fraction_tpm_ge_0.1": f"{detected_01 / len(candidate_ids):.6f}" if candidate_ids else "",
        })

    expression_fields = [
        "group", "species", "run_accession", "gene_id", "orthogroup", "validation_status",
        "validation_flags", "length", "effective_length", "tpm", "num_reads",
        "expression_support", "quantification_status", "best_event_direction",
        "best_event_change", "selection_score", "busted_q_bh", "relax_q_bh",
        "relax_interpretation",
    ]
    write_tsv(run / "candidate_expression.tsv", expression_fields, expression_rows)
    species_fields = list(species_rows[0].keys())
    write_tsv(run / "species_expression_summary.tsv", species_fields, species_rows)

    grouped: dict[tuple[str, str, str], list[dict]] = defaultdict(list)
    for row in expression_rows:
        if row["orthogroup"]:
            grouped[(row["group"], row["orthogroup"], row["species"])].append(row)
    og_species_rows: list[dict] = []
    for (group, orthogroup, species), rows in sorted(grouped.items()):
        tpms = [float(row["tpm"]) for row in rows]
        og_species_rows.append({
            "group": group,
            "orthogroup": orthogroup,
            "species": species,
            "candidate_copy_count": len(rows),
            "detected_copy_count_tpm_ge_0.1": sum(tpm >= 0.1 for tpm in tpms),
            "sum_tpm": f"{sum(tpms):.6f}",
            "max_tpm": f"{max(tpms):.6f}",
            "mean_tpm": f"{mean(tpms):.6f}",
        })
    write_tsv(
        run / "orthogroup_species_expression.tsv",
        ["group", "orthogroup", "species", "candidate_copy_count", "detected_copy_count_tpm_ge_0.1", "sum_tpm", "max_tpm", "mean_tpm"],
        og_species_rows,
    )

    family_grouped: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for row in og_species_rows:
        family_grouped[(row["group"], row["orthogroup"])].append(row)
    integrated_rows: list[dict] = []
    for key, rows in sorted(family_grouped.items()):
        group, orthogroup = key
        evo = family_evolution.get(key, {})
        integrated_rows.append({
            "group": group,
            "orthogroup": orthogroup,
            "species_with_expression_data": len(rows),
            "species_with_detected_candidate": sum(int(row["detected_copy_count_tpm_ge_0.1"]) > 0 for row in rows),
            "candidate_copy_count": sum(int(row["candidate_copy_count"]) for row in rows),
            "detected_copy_count_tpm_ge_0.1": sum(int(row["detected_copy_count_tpm_ge_0.1"]) for row in rows),
            "max_within_library_tpm": f"{max(float(row['max_tpm']) for row in rows):.6f}",
            "best_event_direction": evo.get("best_event_direction", ""),
            "best_event_change": evo.get("best_event_change", ""),
            "selection_score": evo.get("selection_score", ""),
            "busted_q_bh": evo.get("busted_q_bh", ""),
            "absrel_significant_branches_raw": evo.get("absrel_significant_branches_raw", ""),
            "relax_q_bh": evo.get("relax_q_bh", ""),
            "relax_interpretation": evo.get("relax_interpretation", ""),
        })
    write_tsv(
        run / "orthogroup_transcriptomic_integration.tsv",
        [
            "group", "orthogroup", "species_with_expression_data", "species_with_detected_candidate",
            "candidate_copy_count", "detected_copy_count_tpm_ge_0.1", "max_within_library_tpm",
            "best_event_direction", "best_event_change", "selection_score", "busted_q_bh",
            "absrel_significant_branches_raw", "relax_q_bh", "relax_interpretation",
        ],
        integrated_rows,
    )

    de_rows = [{
        "group": row["group"],
        "species": row["species"],
        "run_accession": row["run_accession"],
        "status": "not_performed",
        "reason": "single library; no biological replicates or condition contrast",
    } for row in manifest]
    write_tsv(
        run / "differential_expression_status.tsv",
        ["group", "species", "run_accession", "status", "reason"],
        de_rows,
    )

    step13_rows: list[dict] = []
    for group in ("bivalvia", "gastropoda"):
        group_species = [row for row in species_rows if row["group"] == group]
        step13_rows.append({
            "group": group,
            "libraries": len(group_species),
            "core_species_with_data": len(group_species),
            "candidate_transcripts": sum(int(row["candidate_transcripts"]) for row in group_species),
            "candidate_detected_tpm_ge_0.1": sum(int(row["candidate_detected_tpm_ge_0.1"]) for row in group_species),
            "candidate_detected_tpm_ge_1": sum(int(row["candidate_detected_tpm_ge_1"]) for row in group_species),
            "candidate_detected_tpm_ge_10": sum(int(row["candidate_detected_tpm_ge_10"]) for row in group_species),
            "differential_expression": "not_performed_single_library",
        })
    write_tsv(
        run / "step13_summary.tsv",
        ["group", "libraries", "core_species_with_data", "candidate_transcripts", "candidate_detected_tpm_ge_0.1", "candidate_detected_tpm_ge_1", "candidate_detected_tpm_ge_10", "differential_expression"],
        step13_rows,
    )

    manifest_rows = []
    for path in sorted(run.iterdir()):
        if path.is_file() and path.name not in {"COMPLETE", "result_manifest.tsv"}:
            manifest_rows.append({"file": path.name, "bytes": path.stat().st_size})
    write_tsv(run / "result_manifest.tsv", ["file", "bytes"], manifest_rows)
    (run / "COMPLETE").write_text(f"completed_run\t{run_id}\nlibraries\t{len(manifest)}\n")
    replace_symlink(BASE / "latest", run)
    print(f"Completed Step 13 summary: {run}")


def main() -> None:
    if len(sys.argv) != 3 or sys.argv[1] not in {"prepare", "summarise"}:
        raise SystemExit("Usage: 30_step13_transcriptomics.py {prepare|summarise} RUN_ID")
    if sys.argv[1] == "prepare":
        prepare(sys.argv[2])
    else:
        summarise(sys.argv[2])


if __name__ == "__main__":
    main()
