#!/usr/bin/env python3
"""Step 12: codon-aware molecular evolution analysis preparation and summary.

The preparation stage maps the exact protein identifiers used by Step 10 to
their source CDS records, validates translation, constructs protein-guided
codon alignments, prunes the corresponding gene trees, and creates HyPhy
manifests.  The summary stage parses BUSTED, aBSREL and RELAX JSON outputs and
applies Benjamini-Hochberg correction across families/branches.
"""

from __future__ import annotations

import argparse
import copy
import csv
import gzip
import io
import json
import math
import re
import shutil
import zipfile
from collections import Counter, defaultdict
from pathlib import Path

from Bio import Phylo
from Bio.Seq import Seq


SPECIES_SOURCES = {
    "Acanthocardia_echinata": ("transdecoder", "predicted_proteomes/Acanthocardia_echinata.transcriptome_cds.fna"),
    "Americardia_media": ("transdecoder", "predicted_proteomes/Americardia_media.transcriptome_cds.fna"),
    "Bullacta_exarata": ("transdecoder", "predicted_proteomes/Bullacta_exarata.transcriptome_cds.fna"),
    "Elysia_crispata": ("transdecoder", "predicted_proteomes/Elysia_crispata.transcriptome_cds.fna"),
    "Tridacna_gigas": ("transdecoder", "predicted_proteomes/Tridacna_gigas.transcriptome_cds.fna"),
    "Tridacna_maxima": ("transdecoder", "predicted_proteomes/Tridacna_maxima.transcriptome_cds.fna"),
    "Aplysia_californica": ("ncbi_zip", "genomes/Gastropoda/Aplysia_californica/GCF_000002075.1.zip"),
    "Berghia_stephanieae": ("ncbi_zip", "genomes/Gastropoda/Berghia_stephanieae/GCF_034508935.2.zip"),
    "Elysia_chlorotica": ("ncbi_zip", "genomes/Gastropoda/Elysia_chlorotica/GCA_003991915.1.zip"),
    "Elysia_marginata": ("ncbi_zip", "genomes/Gastropoda/Elysia_marginata/GCA_019649035.1.zip"),
    "Plakobranchus_ocellatus": ("ncbi_zip", "genomes/Gastropoda/Plakobranchus_ocellatus/GCA_019648995.1.zip"),
    "Cerastoderma_edule": (
        "ensembl_gff",
        "genomes/Bivalvia/Cerastoderma_edule/Cerastoderma_edule_GCA_963989375.1_Ensembl.gff3.gz",
        "genome_fastas/Cerastoderma_edule.fna",
    ),
    "Fragum_sueziense": (
        "ensembl_gff",
        "source_downloads/Fragum_sueziense/GCA_963680895.1/Fragum_sueziense-GCA_963680895.1-2024_01-genes.gff3.gz",
        "genome_fastas/Fragum_sueziense.fna",
    ),
    "Onchidella_celtica": (
        "ensembl_gff",
        "genomes/Gastropoda/Onchidella_celtica/Onchidella_celtica_GCA_963931925.1_Ensembl.gff3.gz",
        "genome_fastas/Onchidella_celtica.fna",
    ),
}

MANIFEST_FIELDS = [
    "array_index", "group", "rank", "orthogroup", "status", "reason",
    "aligned_protein", "source_tree", "codon_alignment", "clean_tree", "relax_tree",
    "output_dir", "input_sequences", "mapped_sequences", "retained_sequences",
    "retained_species", "codon_sites", "cds_mapping_fraction", "foreground_node",
    "foreground_species", "test_branches", "reference_branches", "run_relax",
    "best_event_direction", "best_event_change", "selection_score",
]


def read_tsv(path: Path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8", errors="backslashreplace") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def iter_fasta(handle):
    name = None
    description = ""
    chunks = []
    for raw in handle:
        if isinstance(raw, bytes):
            raw = raw.decode("utf-8", errors="replace")
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                yield name, description, "".join(chunks).upper()
            description = line[1:]
            name = description.split()[0]
            chunks = []
        else:
            if name is None:
                raise ValueError("Sequence appeared before a FASTA header")
            chunks.append(line)
    if name is not None:
        yield name, description, "".join(chunks).upper()


def read_fasta(path: Path):
    with path.open() as handle:
        return {name: seq for name, _, seq in iter_fasta(handle)}


def write_fasta(path: Path, records):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for name, seq in records:
            handle.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                handle.write(seq[i:i + 80] + "\n")


def parse_attrs(text):
    attrs = {}
    for item in text.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            attrs[key] = value
    return attrs


def target_base_map(targets):
    mapping = defaultdict(list)
    for target in targets:
        mapping[target.rsplit(".", 1)[0]].append(target)
    return mapping


def extract_transdecoder(path: Path, targets):
    found = {}
    with path.open() as handle:
        for name, _, seq in iter_fasta(handle):
            if name in targets:
                found[name] = seq
    return found


def extract_ncbi_zip(path: Path, targets):
    found = {}
    protein_re = re.compile(r"\[protein_id=([^\]]+)\]")
    with zipfile.ZipFile(path) as archive:
        members = [name for name in archive.namelist() if name.endswith("/cds_from_genomic.fna")]
        if len(members) != 1:
            raise ValueError(f"Expected one cds_from_genomic.fna in {path}, found {members}")
        with archive.open(members[0]) as raw:
            stream = io.TextIOWrapper(raw, encoding="utf-8", errors="replace")
            for _, description, seq in iter_fasta(stream):
                match = protein_re.search(description)
                if match and match.group(1) in targets:
                    found[match.group(1)] = seq
    return found


def reverse_complement(seq):
    return str(Seq(seq).reverse_complement())


def extract_ensembl_gff(gff_path: Path, genome_path: Path, targets):
    base_targets = target_base_map(targets)
    fragments = defaultdict(list)
    needed_seqids = set()
    opener = gzip.open if gff_path.suffix == ".gz" else open
    with opener(gff_path, "rt") as handle:
        for raw in handle:
            if not raw or raw.startswith("#"):
                continue
            cols = raw.rstrip("\n").split("\t")
            if len(cols) != 9 or cols[2] != "CDS":
                continue
            attrs = parse_attrs(cols[8])
            base = attrs.get("protein_id", "")
            if base not in base_targets:
                continue
            version = attrs.get("version", "")
            full = f"{base}.{version}" if version else base
            candidates = base_targets[base]
            target = full if full in candidates else (base if base in candidates else candidates[0])
            seqid, start, end, strand = cols[0], int(cols[3]) - 1, int(cols[4]), cols[6]
            fragments[target].append([seqid, start, end, strand, None])
            needed_seqids.add(seqid)

    by_seqid = defaultdict(list)
    for protein, rows in fragments.items():
        for row in rows:
            by_seqid[row[0]].append((protein, row))

    seen_seqids = set()
    with genome_path.open() as handle:
        for record_id, description, sequence in iter_fasta(handle):
            aliases = {record_id}
            match = re.search(r"chromosome:\s*([^,\s]+)", description)
            if match:
                aliases.add(match.group(1))
            matched = aliases & needed_seqids
            for seqid in matched:
                seen_seqids.add(seqid)
                for _, row in by_seqid[seqid]:
                    row[4] = sequence[row[1]:row[2]]

    missing_seqids = needed_seqids - seen_seqids
    if missing_seqids:
        raise ValueError(f"Genome/GFF sequence ID mismatch for {gff_path}: {sorted(missing_seqids)[:10]}")

    found = {}
    for protein, rows in fragments.items():
        strands = {row[3] for row in rows}
        seqids = {row[0] for row in rows}
        if len(strands) != 1 or len(seqids) != 1 or any(row[4] is None for row in rows):
            continue
        strand = rows[0][3]
        rows.sort(key=lambda row: row[1], reverse=(strand == "-"))
        pieces = [reverse_complement(row[4]) if strand == "-" else row[4] for row in rows]
        found[protein] = "".join(pieces)
    return found


def species_gene(label):
    if "__" not in label:
        return "", label
    return label.split("__", 1)


def clean_cds(seq):
    return re.sub(r"[^ACGTN]", "N", seq.upper().replace("U", "T"))


def translated(seq):
    seq = clean_cds(seq)
    if len(seq) % 3:
        return None
    return str(Seq(seq).translate(table=1)).rstrip("*")


def compatible_translation(protein, cds):
    protein = protein.replace("-", "").replace("*", "").upper()
    aa = translated(cds)
    if aa is None or "*" in aa or len(aa) != len(protein):
        return False, aa or ""
    for index, (observed, expected) in enumerate(zip(aa, protein)):
        if observed == expected or observed == "X" or expected == "X":
            continue
        if index == 0 and expected == "M":
            continue
        return False, aa
    return True, aa


def backtranslate(aligned_protein, cds):
    cds = clean_cds(cds)
    if cds[-3:] in {"TAA", "TAG", "TGA"}:
        cds = cds[:-3]
    codons = [cds[i:i + 3] for i in range(0, len(cds), 3)]
    result = []
    index = 0
    for aa in aligned_protein:
        if aa == "-":
            result.append("---")
        else:
            if index >= len(codons):
                raise ValueError("Protein alignment consumes more codons than available")
            result.append(codons[index])
            index += 1
    if index != len(codons):
        raise ValueError("Protein alignment did not consume the full CDS")
    return result


def trim_codon_alignment(codon_lists, min_valid_fraction=0.50):
    if not codon_lists:
        return {}, 0
    names = list(codon_lists)
    length = len(codon_lists[names[0]])
    stop_codons = {"TAA", "TAG", "TGA"}
    keep = []
    required = max(3, math.ceil(len(names) * min_valid_fraction))
    for index in range(length):
        valid = 0
        for name in names:
            codon = codon_lists[name][index]
            if codon != "---" and "N" not in codon and codon not in stop_codons:
                valid += 1
        if valid >= required:
            keep.append(index)
    trimmed = {
        name: "".join(codons[index] for index in keep)
        for name, codons in codon_lists.items()
    }
    return trimmed, len(keep)


def prune_tree(source: Path, retained, output: Path):
    tree = Phylo.read(str(source), "newick")
    retained = set(retained)
    observed = {tip.name for tip in tree.get_terminals()}
    for name in sorted(observed - retained):
        tree.prune(name)
    remaining = {tip.name for tip in tree.get_terminals()}
    if remaining != retained:
        raise ValueError(f"Pruned tree mismatch; missing={sorted(retained - remaining)[:5]}")
    for clade in tree.get_nonterminals():
        clade.name = None
        clade.confidence = None
    tree.rooted = False
    output.parent.mkdir(parents=True, exist_ok=True)
    Phylo.write(tree, str(output), "newick", plain=False)
    return tree


def write_relax_tree(tree, test_species, output: Path):
    relax = copy.deepcopy(tree)
    test_count = 0
    reference_count = 0
    for tip in relax.get_terminals():
        species, _ = species_gene(tip.name)
        if species in test_species:
            tip.name = f"{tip.name}{{test}}"
            test_count += 1
        else:
            tip.name = f"{tip.name}{{reference}}"
            reference_count += 1
    Phylo.write(relax, str(output), "newick", plain=False)
    return test_count, reference_count


def event_species(group, node, branch_lookup, observed_species):
    base = re.sub(r"<[^>]*>$", "", node or "")
    if (group, base) in branch_lookup:
        return set(filter(None, branch_lookup[(group, base)].split(";")))
    if base in observed_species:
        return {base}
    return set()


def make_cds_database(root: Path, targets_by_species, output_dir: Path):
    rows = []
    all_cds = {}
    cds_dir = output_dir / "cds_by_species"
    cds_dir.mkdir(parents=True, exist_ok=True)
    for species in sorted(targets_by_species):
        targets = set(targets_by_species[species])
        source = SPECIES_SOURCES.get(species)
        if source is None:
            found = {}
            source_type = "missing_configuration"
            source_path = ""
        else:
            source_type = source[0]
            source_path = root / source[1]
            if source_type == "transdecoder":
                found = extract_transdecoder(source_path, targets)
            elif source_type == "ncbi_zip":
                found = extract_ncbi_zip(source_path, targets)
            elif source_type == "ensembl_gff":
                found = extract_ensembl_gff(source_path, root / source[2], targets)
            else:
                raise ValueError(f"Unsupported source type {source_type}")
        found = {name: clean_cds(seq) for name, seq in found.items()}
        all_cds[species] = found
        write_fasta(cds_dir / f"{species}.candidate_cds.fna", sorted(found.items()))
        rows.append({
            "species": species,
            "source_type": source_type,
            "source_path": str(source_path),
            "requested_candidate_ids": len(targets),
            "extracted_candidate_cds": len(found),
            "missing_candidate_cds": len(targets - set(found)),
            "coverage_fraction": f"{len(found) / len(targets):.6f}" if targets else "",
        })
    write_tsv(output_dir / "cds_source_summary.tsv", [
        "species", "source_type", "source_path", "requested_candidate_ids",
        "extracted_candidate_cds", "missing_candidate_cds", "coverage_fraction",
    ], rows)
    return all_cds


def prepare(args):
    root = Path(args.root)
    step10 = Path(args.step10).resolve()
    step11 = Path(args.step11).resolve()
    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)

    gene_rows = [row for row in read_tsv(step10 / "gene_tree_summary.tsv") if row.get("status") == "complete"]
    reconciliation = {
        (row["group"], row["orthogroup"]): row
        for row in read_tsv(step11 / "family_reconciliation_summary.tsv")
    }
    branch_lookup = {
        (row["group"], row["species_node"]): row["descendant_species"]
        for row in read_tsv(step11 / "species_tree_branches.tsv")
    }

    family_inputs = []
    targets_by_species = defaultdict(set)
    for row in gene_rows:
        treefile = Path(row["treefile"])
        aligned = treefile.parent / "aligned.faa"
        if not aligned.is_file():
            continue
        proteins = read_fasta(aligned)
        for label in proteins:
            species, gene = species_gene(label)
            targets_by_species[species].add(gene)
        family_inputs.append((row, aligned, treefile, proteins))

    all_cds = make_cds_database(root, targets_by_species, output)
    input_root = output / "inputs"
    result_root = output / "hyphy_results"
    preparation_rows = []
    manifests = {"bivalvia": [], "gastropoda": []}

    family_inputs.sort(key=lambda item: (item[0]["group"], int(float(item[0]["rank"]))))
    group_indices = Counter()
    for row, aligned_path, treefile, proteins in family_inputs:
        group = row["group"]
        orthogroup = row["orthogroup"]
        rank = int(float(row["rank"]))
        rec = reconciliation.get((group, orthogroup), {})
        family_dir = input_root / group / orthogroup
        family_dir.mkdir(parents=True, exist_ok=True)
        codon_lists = {}
        validation_rows = []
        mapped = 0
        for label, aligned_protein in proteins.items():
            species, gene = species_gene(label)
            cds = all_cds.get(species, {}).get(gene)
            status = "missing_cds"
            detail = ""
            if cds is not None:
                mapped += 1
                ok, observed = compatible_translation(aligned_protein, cds)
                if ok:
                    try:
                        codon_lists[label] = backtranslate(aligned_protein, cds)
                        status = "retained"
                    except ValueError as exc:
                        status = "backtranslation_failed"
                        detail = str(exc)
                else:
                    status = "translation_mismatch"
                    detail = f"protein_aa={len(aligned_protein.replace('-', ''))};translated_aa={len(observed)}"
            validation_rows.append({
                "tree_label": label, "species": species, "gene_id": gene,
                "status": status, "detail": detail,
            })

        trimmed, codon_sites = trim_codon_alignment(codon_lists)
        retained = [name for name, seq in trimmed.items() if seq and seq.count("---") / max(1, codon_sites) <= 0.50]
        trimmed = {name: trimmed[name] for name in retained}
        observed_species = {species_gene(name)[0] for name in retained}
        mapping_fraction = mapped / len(proteins) if proteins else 0.0
        reasons = []
        if len(retained) < args.min_sequences:
            reasons.append(f"retained_sequences<{args.min_sequences}")
        if len(observed_species) < args.min_species:
            reasons.append(f"retained_species<{args.min_species}")
        if codon_sites < args.min_codon_sites:
            reasons.append(f"codon_sites<{args.min_codon_sites}")
        if mapping_fraction < args.min_mapping_fraction:
            reasons.append(f"cds_mapping_fraction<{args.min_mapping_fraction}")
        status = "eligible" if not reasons else "skipped"

        codon_path = family_dir / f"{orthogroup}.codon.fasta"
        tree_path = family_dir / f"{orthogroup}.clean.nwk"
        relax_path = family_dir / f"{orthogroup}.relax.nwk"
        write_tsv(family_dir / "sequence_cds_validation.tsv",
                  ["tree_label", "species", "gene_id", "status", "detail"], validation_rows)
        run_relax = "no"
        test_count = 0
        ref_count = 0
        foreground = set()
        if status == "eligible":
            write_fasta(codon_path, sorted(trimmed.items()))
            tree = prune_tree(treefile, retained, tree_path)
            foreground = event_species(group, rec.get("cafe_best_event_node", ""), branch_lookup, observed_species)
            test_count, ref_count = write_relax_tree(tree, foreground, relax_path)
            run_relax = "yes" if test_count >= 2 and ref_count >= 2 else "no"
        else:
            codon_path = Path("")
            tree_path = Path("")
            relax_path = Path("")

        array_index = group_indices[group]
        group_indices[group] += 1
        manifest_row = {
            "array_index": array_index,
            "group": group,
            "rank": rank,
            "orthogroup": orthogroup,
            "status": status,
            "reason": ";".join(reasons),
            "aligned_protein": str(aligned_path),
            "source_tree": str(treefile),
            "codon_alignment": str(codon_path) if status == "eligible" else "",
            "clean_tree": str(tree_path) if status == "eligible" else "",
            "relax_tree": str(relax_path) if status == "eligible" else "",
            "output_dir": str(result_root / group / orthogroup),
            "input_sequences": len(proteins),
            "mapped_sequences": mapped,
            "retained_sequences": len(retained),
            "retained_species": len(observed_species),
            "codon_sites": codon_sites,
            "cds_mapping_fraction": f"{mapping_fraction:.6f}",
            "foreground_node": rec.get("cafe_best_event_node", ""),
            "foreground_species": ";".join(sorted(foreground)),
            "test_branches": test_count,
            "reference_branches": ref_count,
            "run_relax": run_relax,
            "best_event_direction": rec.get("cafe_best_event_direction", row.get("best_event_direction", "")),
            "best_event_change": rec.get("cafe_best_event_change", row.get("best_event_change", "")),
            "selection_score": rec.get("selection_score", row.get("selection_score", "")),
        }
        manifests[group].append(manifest_row)
        preparation_rows.append(manifest_row)

    for group, rows in manifests.items():
        rows.sort(key=lambda row: int(row["array_index"]))
        write_tsv(output / f"manifest_{group}.tsv", MANIFEST_FIELDS, rows)
    write_tsv(output / "family_preparation_summary.tsv", MANIFEST_FIELDS, preparation_rows)

    methods = f"""Step 12 molecular-evolution input preparation

Protein identifiers are mapped to exact source CDS records from six TransDecoder outputs, five NCBI
cds_from_genomic archives, and three Ensembl GFF3/genome pairs. Each CDS is translated with the universal
nuclear genetic code and must match its Step 10 protein (X is treated as an ambiguity; an alternative first
codon is allowed when the annotated protein starts with M). Protein MAFFT alignments are back-translated to
codons, and columns retained by at least 50% of sequences are kept.

Family eligibility thresholds: at least {args.min_sequences} retained sequences, {args.min_species} species,
{args.min_codon_sites} codons, and CDS mapping fraction >= {args.min_mapping_fraction:.2f}. BUSTED tests
gene-wide episodic positive selection; aBSREL tests individual branches. RELAX compares terminal gene branches
from the CAFE best-event descendant species against remaining terminal branches, only when both sets contain
at least two branches. RELAX K>1 indicates intensified selection and K<1 indicates relaxed selection; RELAX is
not itself a positive-selection test.
"""
    (output / "METHODS_AND_LIMITATIONS.txt").write_text(methods)
    (output / "PREPARED").write_text("Step 12 inputs prepared successfully\n")


def find_key(obj, patterns):
    patterns = [pattern.lower() for pattern in patterns]
    matches = []

    def walk(value, path=()):
        if isinstance(value, dict):
            for key, child in value.items():
                new_path = path + (str(key),)
                key_lower = str(key).lower()
                if any(pattern == key_lower or pattern in key_lower for pattern in patterns):
                    matches.append((new_path, child))
                walk(child, new_path)
        elif isinstance(value, list):
            for index, child in enumerate(value):
                walk(child, path + (str(index),))

    walk(obj)
    return matches


def first_number(obj, patterns):
    for _, value in find_key(obj, patterns):
        if isinstance(value, (int, float)) and math.isfinite(float(value)):
            return float(value)
    return None


def load_json(path: Path):
    if not path.is_file():
        return None
    try:
        with path.open() as handle:
            return json.load(handle)
    except (json.JSONDecodeError, OSError):
        return None


def bh_adjust(values):
    indexed = [(index, value) for index, value in enumerate(values) if value is not None and math.isfinite(value)]
    indexed.sort(key=lambda item: item[1])
    adjusted = [None] * len(values)
    running = 1.0
    total = len(indexed)
    for reverse_rank, (index, value) in enumerate(reversed(indexed), 1):
        rank = total - reverse_rank + 1
        running = min(running, value * total / rank)
        adjusted[index] = min(1.0, running)
    return adjusted


def parse_absrel_branches(data, group, orthogroup):
    rows = []
    if not isinstance(data, dict):
        return rows
    branch_attributes = data.get("branch attributes", {})
    if isinstance(branch_attributes, dict) and branch_attributes:
        partition = next(iter(branch_attributes.values()))
        if isinstance(partition, dict):
            for branch, attrs in partition.items():
                if not isinstance(attrs, dict):
                    continue
                pvalue = None
                for key, value in attrs.items():
                    if "corrected p-value" in key.lower() and isinstance(value, (int, float)):
                        pvalue = float(value)
                        break
                omega = None
                for key, value in attrs.items():
                    if "omega" in key.lower() and isinstance(value, (int, float)):
                        omega = float(value)
                        break
                rows.append({
                    "group": group, "orthogroup": orthogroup, "branch": branch,
                    "absrel_corrected_p": pvalue, "omega_summary": omega,
                })
    return rows


def summarise(args):
    run = Path(args.run_dir).resolve()
    preparation = read_tsv(run / "family_preparation_summary.tsv")
    family_rows = []
    branch_rows = []
    for prep in preparation:
        group = prep["group"]
        orthogroup = prep["orthogroup"]
        out = Path(prep["output_dir"])
        busted = load_json(out / "busted.json")
        absrel = load_json(out / "absrel.json")
        relax = load_json(out / "relax.json")
        busted_test = busted.get("test results", {}) if isinstance(busted, dict) else {}
        relax_test = relax.get("test results", {}) if isinstance(relax, dict) else {}
        busted_p = first_number(busted_test, ["p-value"]) if busted_test else None
        busted_lrt = first_number(busted_test, ["lrt"]) if busted_test else None
        relax_p = first_number(relax_test, ["p-value"]) if relax_test else None
        relax_k = first_number(relax_test, ["relaxation or intensification parameter", "selection intensity parameter", "k (selection intensity)"]) if relax_test else None
        branches = parse_absrel_branches(absrel, group, orthogroup)
        branch_rows.extend(branches)
        family_rows.append({
            **prep,
            "execution_status": "complete" if (out / "COMPLETE").is_file() else ("failed" if (out / "FAILED").is_file() else prep["status"]),
            "busted_p": busted_p,
            "busted_lrt": busted_lrt,
            "absrel_tested_branches": len(branches),
            "absrel_significant_branches_raw": sum(1 for row in branches if row["absrel_corrected_p"] is not None and row["absrel_corrected_p"] <= 0.05),
            "relax_k": relax_k,
            "relax_p": relax_p,
            "busted_json": str(out / "busted.json") if busted else "",
            "absrel_json": str(out / "absrel.json") if absrel else "",
            "relax_json": str(out / "relax.json") if relax else "",
        })

    busted_q = bh_adjust([row["busted_p"] for row in family_rows])
    relax_q = bh_adjust([row["relax_p"] for row in family_rows])
    for row, bq, rq in zip(family_rows, busted_q, relax_q):
        row["busted_q_bh"] = bq
        row["relax_q_bh"] = rq
        if row["relax_k"] is None:
            row["relax_interpretation"] = ""
        elif row["relax_k"] > 1:
            row["relax_interpretation"] = "intensified"
        elif row["relax_k"] < 1:
            row["relax_interpretation"] = "relaxed"
        else:
            row["relax_interpretation"] = "unchanged"

    branch_q = bh_adjust([row["absrel_corrected_p"] for row in branch_rows])
    for row, qvalue in zip(branch_rows, branch_q):
        row["global_q_bh"] = qvalue

    family_fields = MANIFEST_FIELDS + [
        "execution_status", "busted_p", "busted_lrt", "busted_q_bh",
        "absrel_tested_branches", "absrel_significant_branches_raw",
        "relax_k", "relax_p", "relax_q_bh", "relax_interpretation",
        "busted_json", "absrel_json", "relax_json",
    ]
    write_tsv(run / "family_molecular_evolution_results.tsv", family_fields, family_rows)
    write_tsv(run / "absrel_branch_results.tsv",
              ["group", "orthogroup", "branch", "absrel_corrected_p", "global_q_bh", "omega_summary"], branch_rows)
    significant_families = [
        row for row in family_rows
        if (row["busted_q_bh"] is not None and row["busted_q_bh"] <= 0.05)
        or (row["relax_q_bh"] is not None and row["relax_q_bh"] <= 0.05)
    ]
    significant_branches = [row for row in branch_rows if row["global_q_bh"] is not None and row["global_q_bh"] <= 0.05]
    write_tsv(run / "significant_family_results.tsv", family_fields, significant_families)
    write_tsv(run / "significant_absrel_branches.tsv",
              ["group", "orthogroup", "branch", "absrel_corrected_p", "global_q_bh", "omega_summary"], significant_branches)

    summary_rows = []
    for group in ("bivalvia", "gastropoda"):
        subset = [row for row in family_rows if row["group"] == group]
        b_sig = sum(row["busted_q_bh"] is not None and row["busted_q_bh"] <= 0.05 for row in subset)
        r_sig = [row for row in subset if row["relax_q_bh"] is not None and row["relax_q_bh"] <= 0.05]
        summary_rows.append({
            "group": group,
            "families_prepared": len(subset),
            "families_eligible": sum(row["status"] == "eligible" for row in subset),
            "families_completed": sum(row["execution_status"] == "complete" for row in subset),
            "busted_fdr_significant": b_sig,
            "absrel_global_fdr_significant_branches": sum(row["group"] == group for row in significant_branches),
            "relax_fdr_significant": len(r_sig),
            "relax_significant_intensified": sum(row["relax_interpretation"] == "intensified" for row in r_sig),
            "relax_significant_relaxed": sum(row["relax_interpretation"] == "relaxed" for row in r_sig),
        })
    summary_fields = [
        "group", "families_prepared", "families_eligible", "families_completed",
        "busted_fdr_significant", "absrel_global_fdr_significant_branches",
        "relax_fdr_significant", "relax_significant_intensified", "relax_significant_relaxed",
    ]
    write_tsv(run / "step12_summary.tsv", summary_fields, summary_rows)

    methods = (run / "METHODS_AND_LIMITATIONS.txt").read_text()
    methods += """
HyPhy result interpretation

BUSTED family P-values and RELAX family P-values were independently adjusted across all tested families using
Benjamini-Hochberg FDR. aBSREL corrected branch P-values were additionally adjusted across all tested branches
from all families. A significant BUSTED result indicates evidence of episodic diversifying selection somewhere
in the gene tree, not necessarily on the CAFE focal lineage. aBSREL localizes candidate branches but is sensitive
to alignment and gene-tree error. RELAX K>1 means intensified selection and K<1 means relaxed selection; neither
direction alone demonstrates adaptive positive selection. Families with multiple optimal roots in Step 11 and
families dominated by partial transcript-derived sequences require additional manual scrutiny.
"""
    (run / "METHODS_AND_LIMITATIONS.txt").write_text(methods)

    manifest_rows = []
    for path in sorted(run.rglob("*")):
        if path.is_file() and path.name != "result_manifest.tsv":
            manifest_rows.append({"relative_path": str(path.relative_to(run)), "bytes": path.stat().st_size})
    write_tsv(run / "result_manifest.tsv", ["relative_path", "bytes"], manifest_rows)
    (run / "COMPLETE").write_text("Step 12 molecular evolution analysis complete\n")


def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)

    prep = subparsers.add_parser("prepare")
    prep.add_argument("--root", required=True)
    prep.add_argument("--step10", required=True)
    prep.add_argument("--step11", required=True)
    prep.add_argument("--output-dir", required=True)
    prep.add_argument("--min-sequences", type=int, default=5)
    prep.add_argument("--min-species", type=int, default=3)
    prep.add_argument("--min-codon-sites", type=int, default=100)
    prep.add_argument("--min-mapping-fraction", type=float, default=0.70)
    prep.set_defaults(func=prepare)

    summary = subparsers.add_parser("summarise")
    summary.add_argument("--run-dir", required=True)
    summary.set_defaults(func=summarise)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
