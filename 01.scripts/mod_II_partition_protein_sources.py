#!/usr/bin/env python3
"""
maniFasta: partition a combined protein accession list by fetch backend.

Reads the staged combined protein list (as written by runBuild's source-plan
staging) and splits each row to the right fetcher's input based on its
`fetch_source` value: per-row column wins, else --default-source. NCBI rows go
to --ncbi-out, UniProt rows to --uniprot-out, UniParc rows to --uniparc-out,
and PDB rows to --pdb-out. Anything still unimplemented is written to
--pending-out and reported rather than silently dropped. Output rows preserve
the input columns verbatim.

A file with no fetch_source column (or all-NCBI) is a safe pass-through:
everything lands in --ncbi-out, exactly as before.
"""
import argparse
import csv
import sys
from pathlib import Path

VALID_FETCH_SOURCES = {"ncbi", "uniprot", "uniparc", "pdb"}
IMPLEMENTED = {"ncbi", "uniprot", "uniparc", "pdb"}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-i", "--input", required=True, type=Path)
    ap.add_argument("--default-source", default="ncbi")
    ap.add_argument("--ncbi-out", required=True, type=Path)
    ap.add_argument("--uniprot-out", required=True, type=Path)
    ap.add_argument("--uniparc-out", required=True, type=Path)
    ap.add_argument("--pdb-out", required=True, type=Path)
    ap.add_argument("--pending-out", required=True, type=Path)
    args = ap.parse_args()

    default = args.default_source.strip().lower()
    if default not in VALID_FETCH_SOURCES:
        sys.exit(f"[ERROR] --default-source must be one of {sorted(VALID_FETCH_SOURCES)}")

    with args.input.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE)
        fields = reader.fieldnames or []
        cols = {c.lower(): c for c in fields}
        if "accession" not in cols and "fetch_id" not in cols:
            sys.exit(f"[ERROR] {args.input}: needs an 'accession' or 'fetch_id' column")
        fs_col = cols.get("fetch_source")
        out = {"ncbi": args.ncbi_out, "uniprot": args.uniprot_out,
               "uniparc": args.uniparc_out, "pdb": args.pdb_out}
        handles, writers, counts = {}, {}, {}
        # Open every distinct destination path (including pending) so each gets
        # a header even when it receives no rows.
        for path in list(out.values()) + [args.pending_out]:
            if path not in handles:
                handles[path] = path.open("w", newline="")
                writers[path] = csv.DictWriter(handles[path], fieldnames=fields,
                                               delimiter="\t", lineterminator="\n")
                writers[path].writeheader()
        for r in reader:
            src = ((r.get(fs_col) or "").strip().lower() if fs_col else "") or default
            if src not in VALID_FETCH_SOURCES:
                acc = r.get(cols.get("accession", "")) or r.get(cols.get("fetch_id", ""))
                sys.exit(f"[ERROR] {args.input}: accession {acc!r} has fetch_source={src!r}; "
                         f"choose {sorted(VALID_FETCH_SOURCES)}")
            dest = out.get(src, args.pending_out)
            writers[dest].writerow(r)
            counts[src] = counts.get(src, 0) + 1
        for h in handles.values():
            h.close()

    for src in sorted(counts, key=lambda s: -counts[s]):
        tag = "" if src in IMPLEMENTED else "  [DECLARED ONLY - no fetcher yet, not fetched]"
        print(f"[INFO] partition: {src:8s} {counts[src]:6d}{tag}", file=sys.stderr)
    pend = sum(c for s, c in counts.items() if s not in IMPLEMENTED)
    if pend:
        print(f"[WARN] {pend} accession(s) routed to {args.pending_out.name} await a "
              f"fetcher; they are NOT in this build.", file=sys.stderr)


if __name__ == "__main__":
    main()
