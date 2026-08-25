#!/usr/bin/env python3

from __future__ import annotations

import argparse
import pathlib
import re


PATTERN = re.compile(
    r"C:(?P<C>[0-9.]+)%\[S:(?P<S>[0-9.]+)%,D:(?P<D>[0-9.]+)%\],"
    r"F:(?P<F>[0-9.]+)%,M:(?P<M>[0-9.]+)%,n:(?P<n>\d+)"
)


def protein_count(path: pathlib.Path) -> int:
    with path.open() as handle:
        return sum(1 for line in handle if line.startswith(">"))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--proteome-dir", type=pathlib.Path, required=True)
    parser.add_argument("--busco-dir", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    args = parser.parse_args()

    columns = ["species", "proteins", "complete_pct", "single_copy_pct", "duplicated_pct", "fragmented_pct", "missing_pct", "busco_n"]
    rows = []
    for fasta in sorted(args.proteome_dir.glob("*.faa")):
        species = fasta.stem
        summaries = sorted((args.busco_dir / species).glob("**/short_summary*.txt"))
        values = {key: "NA" for key in ("C", "S", "D", "F", "M", "n")}
        if summaries:
            text = summaries[0].read_text()
            match = PATTERN.search(text.replace(" ", ""))
            if match:
                values.update(match.groupdict())
        rows.append([species, str(protein_count(fasta)), values["C"], values["S"], values["D"], values["F"], values["M"], values["n"]])

    with args.output.open("w") as out:
        out.write("\t".join(columns) + "\n")
        for row in rows:
            out.write("\t".join(row) + "\n")


if __name__ == "__main__":
    main()
