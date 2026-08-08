#!/usr/bin/env python3
"""
maniFasta: partition a combined proteome/species fetch list by fetch backend.

Reads the staged combined mod_I fetch list (as written by runBuild's
source-plan staging) and splits each row to the right fetcher's input based
on its `fetch_source` value: per-row column wins, else --default-source.
NCBI rows go to --ncbi-out, UniProt rows go to --uniprot-out.

Mirrors mod_II_partition_protein_sources.py, with one difference: Module I
(proteome) sources only ever have two valid backends. UniParc and PDB are
accession-level concepts with no proteome equivalent (see source_planner.py's
validate_rows, which already rejects fetch_source=uniparc/pdb on mod_I rows
at the registry-options level) so an unrecognized per-row value here is a
hard configuration error, not a "not implemented yet" bucket.

A file with no fetch_source column (or all-NCBI) is a safe pass-through:
everything lands in --ncbi-out, exactly as before.
"""
import argparse
import csv
import sys
from pathlib import Path

VALID_FETCH_SOURCES = {"ncbi", "uniprot"}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-i", "--input", required=True, type=Path)
    ap.add_argument("--default-source", default="ncbi")
    ap.add_argument("--ncbi-out", required=True, type=Path)
    ap.add_argument("--uniprot-out", required=True, type=Path)
    args = ap.parse_args()

    default = args.default_source.strip().lower()
    if default not in VALID_FETCH_SOURCES:
        sys.exit(f"[ERROR] --default-source must be one of {sorted(VALID_FETCH_SOURCES)}")

    with args.input.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE)
        fields = reader.fieldnames or []
        cols = {c.lower(): c for c in fields}
        if not ({"species", "accession", "taxonid"} & set(cols)):
            sys.exit(f"[ERROR] {args.input}: needs a 'species', 'accession', or "
                     f"'taxonID' column. Found: {fields}")
        fs_col = cols.get("fetch_source")

        out = {"ncbi": args.ncbi_out, "uniprot": args.uniprot_out}
        handles, writers, counts = {}, {}, {}
        # Open every destination path (so each gets a header even if it
        # receives no rows -- downstream "-s file && grep -c . file > 1"
        # gating in runBuild_db.sh relies on the header always being there).
        for path in out.values():
            handles[path] = path.open("w", newline="")
            writers[path] = csv.DictWriter(handles[path], fieldnames=fields,
                                           delimiter="\t", lineterminator="\n")
            writers[path].writeheader()

        for r in reader:
            src = ((r.get(fs_col) or "").strip().lower() if fs_col else "") or default
            if src not in VALID_FETCH_SOURCES:
                disp = (r.get(cols.get("species", "")) or r.get(cols.get("accession", ""))
                        or r.get(cols.get("taxonid", "")) or "<row>")
                sys.exit(f"[ERROR] {args.input}: row {disp!r} has fetch_source={src!r}; "
                         f"Module I (proteome) sources support only "
                         f"{sorted(VALID_FETCH_SOURCES)} (UniParc/PDB are Module II "
                         f"accession-level backends only)")
            writers[out[src]].writerow(r)
            counts[src] = counts.get(src, 0) + 1

        for h in handles.values():
            h.close()

    for src in sorted(counts, key=lambda s: -counts[s]):
        print(f"[INFO] partition: {src:8s} {counts[src]:6d}", file=sys.stderr)
    if not counts:
        print("[INFO] partition: input had no data rows", file=sys.stderr)


if __name__ == "__main__":
    main()
