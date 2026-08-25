#!/usr/bin/env python3
"""Select one longest protein per explicitly identifiable gene.

The script is deliberately conservative: sequences are collapsed only when
their FASTA headers expose a reliable gene relationship. Otherwise each
sequence ID is retained as its own group.
"""

from __future__ import annotations

import argparse
import collections
import pathlib
import re
from dataclasses import dataclass


@dataclass
class Record:
    seq_id: str
    header: str
    sequence: str
    gene_key: str
    rule: str


def read_fasta(path: pathlib.Path):
    header = None
    chunks: list[str] = []
    with path.open() as handle:
        for raw in handle:
            line = raw.rstrip("\r\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:]
                chunks = []
            elif line.strip():
                if header is None:
                    raise ValueError(f"Sequence before first header in {path}")
                chunks.append(line.strip())
    if header is not None:
        yield header, "".join(chunks)


def infer_gene_key(header: str, seq_id: str) -> tuple[str, str]:
    match = re.search(r"(?:^|\s)gene:([^\s]+)", header)
    if match:
        return f"ensembl_gene:{match.group(1)}", "ensembl_gene"

    match = re.search(r"_g(\d+)_i\d+(?:\.p\d+)?(?:\b|$)", seq_id)
    if match:
        return f"transdecoder_gene:g{match.group(1)}", "transdecoder_gene"

    match = re.search(r"\b(LOC\d+)\b", header)
    if match:
        return f"refseq_locus:{match.group(1)}", "refseq_locus"

    for pattern, name in (
        (r"\b(EGW\d+_\d+)\b", "egw_locus"),
        (r"\b(ElyMa_\d+)\b", "elyma_locus"),
        (r"\b(PoB_\d+)\b", "pob_locus"),
    ):
        match = re.search(pattern, header)
        if match:
            return f"{name}:{match.group(1)}", name

    return f"unique_id:{seq_id}", "unique_id"


def wrap(sequence: str, width: int = 60):
    for start in range(0, len(sequence), width):
        yield sequence[start : start + width]


def process_file(input_path: pathlib.Path, output_path: pathlib.Path, map_path: pathlib.Path):
    best: dict[str, Record] = {}
    removed: list[tuple[Record, str, str]] = []
    rule_sequences = collections.Counter()
    original = 0

    for header, sequence in read_fasta(input_path):
        original += 1
        seq_id = header.split()[0]
        gene_key, rule = infer_gene_key(header, seq_id)
        rule_sequences[rule] += 1
        current = Record(seq_id, header, sequence, gene_key, rule)
        previous = best.get(gene_key)
        if previous is None:
            best[gene_key] = current
            continue

        current_rank = (len(current.sequence), current.seq_id)
        previous_rank = (len(previous.sequence), previous.seq_id)
        if current_rank > previous_rank:
            removed.append((previous, current.seq_id, "shorter_or_tied_isoform"))
            best[gene_key] = current
        else:
            removed.append((current, previous.seq_id, "shorter_or_tied_isoform"))

    output_path.parent.mkdir(parents=True, exist_ok=True)
    map_path.parent.mkdir(parents=True, exist_ok=True)

    representatives = sorted(best.values(), key=lambda record: record.seq_id)
    with output_path.open("w") as out:
        for record in representatives:
            out.write(f">{record.header}\n")
            for line in wrap(record.sequence):
                out.write(f"{line}\n")

    with map_path.open("w") as out:
        out.write("removed_id\trepresentative_id\tgene_key\trule\treason\tremoved_length\n")
        for record, representative_id, reason in sorted(removed, key=lambda item: item[0].seq_id):
            out.write(
                f"{record.seq_id}\t{representative_id}\t{record.gene_key}\t"
                f"{record.rule}\t{reason}\t{len(record.sequence)}\n"
            )

    return {
        "species": input_path.stem,
        "original_proteins": original,
        "representative_proteins": len(representatives),
        "collapsed_isoforms": original - len(representatives),
        "explicitly_grouped_sequences": original - rule_sequences["unique_id"],
        "unique_id_sequences": rule_sequences["unique_id"],
        "rule_counts": ";".join(f"{key}:{value}" for key, value in sorted(rule_sequences.items())),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True, type=pathlib.Path)
    parser.add_argument("--output-dir", required=True, type=pathlib.Path)
    parser.add_argument("--mapping-dir", required=True, type=pathlib.Path)
    parser.add_argument("--summary", required=True, type=pathlib.Path)
    args = parser.parse_args()

    inputs = sorted(args.input_dir.glob("*.faa"))
    if not inputs:
        raise SystemExit(f"No .faa files found in {args.input_dir}")

    rows = []
    for input_path in inputs:
        rows.append(
            process_file(
                input_path,
                args.output_dir / input_path.name,
                args.mapping_dir / f"{input_path.stem}.removed_isoforms.tsv",
            )
        )

    args.summary.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "species",
        "original_proteins",
        "representative_proteins",
        "collapsed_isoforms",
        "explicitly_grouped_sequences",
        "unique_id_sequences",
        "rule_counts",
    ]
    with args.summary.open("w") as out:
        out.write("\t".join(columns) + "\n")
        for row in rows:
            out.write("\t".join(str(row[column]) for column in columns) + "\n")


if __name__ == "__main__":
    main()
