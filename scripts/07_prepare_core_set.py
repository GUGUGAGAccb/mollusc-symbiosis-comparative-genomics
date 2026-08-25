#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import pathlib


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary", type=pathlib.Path, required=True)
    parser.add_argument("--proteome-dir", type=pathlib.Path, required=True)
    parser.add_argument("--core-dir", type=pathlib.Path, required=True)
    parser.add_argument("--decision-table", type=pathlib.Path, required=True)
    parser.add_argument("--min-complete", type=float, default=80.0)
    args = parser.parse_args()

    args.core_dir.mkdir(parents=True, exist_ok=True)
    for existing in args.core_dir.glob("*.faa"):
        existing.unlink()

    decisions = []
    with args.summary.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            complete = float(row["complete_pct"])
            include = complete >= args.min_complete
            reason = f"BUSCO complete {complete:.1f}% {'>=' if include else '<'} {args.min_complete:.1f}%"
            decisions.append((row["species"], complete, "core" if include else "supplementary", reason))
            if include:
                source = args.proteome_dir / f"{row['species']}.faa"
                if not source.is_file():
                    raise FileNotFoundError(source)
                (args.core_dir / source.name).symlink_to(source)

    if sum(1 for row in decisions if row[2] == "core") < 4:
        raise SystemExit("Fewer than four species passed the BUSCO threshold")

    with args.decision_table.open("w") as out:
        out.write("species\tcomplete_pct\tanalysis_role\treason\n")
        for row in decisions:
            out.write("\t".join(map(str, row)) + "\n")


if __name__ == "__main__":
    main()
