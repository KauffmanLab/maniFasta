"""
mod_IV_homd_filter.py
──────────────
Filter a downloaded HOMD taxon TSV by body site(s).
No third-party dependencies — uses only the Python standard library.

Usage
-----
  # See every available body site in the file:
  python mod_IV_homd_filter.py --file taxa.tsv --list-sites

  # Filter to one or more sites (partial, case-insensitive):
  python mod_IV_homd_filter.py --file taxa.tsv --sites Oral
  python mod_IV_homd_filter.py --file taxa.tsv --sites Oral Nasal Skin

  # Save the filtered result:
  python mod_IV_homd_filter.py --file taxa.tsv --sites Oral --output oral_taxa.tsv
"""

import re
import csv
import argparse
from pathlib import Path

_ABUNDANCE_RE = re.compile(r"\s*\(Abundance:[^)]*\)", re.IGNORECASE)


def load_table(filepath: Path) -> tuple[str, list[str], list[dict]]:
    """
    Read the HOMD taxon TSV, preserving the metadata line on row 0.
    Returns (metadata_line, headers, rows) where each row is a dict keyed by column name.
    """
    with filepath.open(encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        metadata = "\t".join(next(reader))  # "HOMD.org Taxon Data::..."
        headers = [h.strip() for h in next(reader)]
        rows = [dict(zip(headers, row)) for row in reader]
    return metadata, headers, rows


def _find_body_site_col(headers: list[str]) -> str:
    for col in headers:
        if "body site" in col.lower():
            return col
    raise RuntimeError(
        f"Could not find a 'Body Site(s)' column.\nAvailable columns: {headers}"
    )


def _parse_sites(raw: str) -> list[str]:
    """
    Split one Body Site(s) cell into a clean list of site names.

    "Oral (Abundance: Low) | Nasal (Abundance: Medium)"
        → ["Oral", "Nasal"]
    """
    parts = [p.strip() for p in raw.split("|") if p.strip()]
    return [_ABUNDANCE_RE.sub("", p).strip() for p in parts if p.strip()]


def get_all_body_sites(headers: list[str], rows: list[dict]) -> list[str]:
    """
    Return every unique body site present in the file, sorted alphabetically,
    with abundance annotations removed.

    Reads directly from the data so new HOMD site categories appear
    automatically without any code changes.
    """
    col = _find_body_site_col(headers)
    sites: set[str] = set()
    for row in rows:
        sites.update(_parse_sites(row.get(col, "")))
    return sorted(s for s in sites if s)


def filter_by_sites(headers: list[str], rows: list[dict], requested: list[str]) -> list[dict]:
    """
    Return only the rows where the taxon is present at ANY of the requested
    body sites. Matching is case-insensitive and substring-based, so "Oral"
    matches "Oral Cavity" if that ever appears in future HOMD updates.
    """
    col = _find_body_site_col(headers)
    requested_lower = [r.lower() for r in requested]

    def row_matches(row: dict) -> bool:
        for site in _parse_sites(row.get(col, "")):
            site_lower = site.lower()
            if any(req in site_lower or site_lower in req for req in requested_lower):
                return True
        return False

    return [row for row in rows if row_matches(row)]


def write_tsv(filepath: Path, metadata: str, headers: list[str], rows: list[dict]):
    with filepath.open("w", encoding="utf-8", newline="") as f:
        f.write(metadata + "\n")
        writer = csv.DictWriter(f, fieldnames=headers, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


# ── CLI ────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Filter a HOMD taxon TSV by body site.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python mod_IV_homd_filter.py --file taxa.tsv --list-sites
  python mod_IV_homd_filter.py --file taxa.tsv --sites Oral
  python mod_IV_homd_filter.py --file taxa.tsv --sites Oral Nasal Skin
  python mod_IV_homd_filter.py --file taxa.tsv --sites Oral --output oral_taxa.tsv
        """,
    )

    parser.add_argument(
        "--file", required=True,
        help="Path to the downloaded HOMD taxon TSV file."
    )
    parser.add_argument(
        "--list-sites", action="store_true",
        help="Print every available body site in the file and exit."
    )
    parser.add_argument(
        "--sites", nargs="*", metavar="SITE",
        help="One or more body sites to filter by (space-separated). "
             "Use 'all' or omit a value to include all taxa. "
             "A taxon is included if it matches ANY of the listed sites."
    )
    parser.add_argument(
        "--output", metavar="FILE",
        help="Optional: save the filtered table to this TSV file."
    )

    args = parser.parse_args()

    path = Path(args.file)
    if not path.exists():
        parser.error(f"File not found: {path}")

    metadata, headers, rows = load_table(path)

    # --list-sites
    if args.list_sites:
        sites = get_all_body_sites(headers, rows)
        print(f"\nAvailable body sites ({len(sites)} total):\n")
        for site in sites:
            print(f"  {site}")
        print()
        return

    # No sites specified
    if not args.sites:
        print("No body site specified — skipping. "
              "Use --list-sites to see options, or --sites SITE [SITE ...] to filter.")
        return

    # --sites all or --sites with no arguments → return everything
    use_all = args.sites == [] or any(s.lower() == "all" for s in args.sites)
    filtered = rows if use_all else filter_by_sites(headers, rows, args.sites)
    site_label = "all" if use_all else str(args.sites)

    print(f"\nMatched {len(filtered)} / {len(rows)} taxa for sites: {site_label}\n")

    col_id      = next((h for h in headers if "HMT" in h), headers[0])
    col_genus   = next((h for h in headers if h.lower() == "genus"), None)
    col_species = next((h for h in headers if h.lower() == "species"), None)
    col_sites   = _find_body_site_col(headers)

    display_cols = [c for c in [col_id, col_genus, col_species, col_sites] if c]
    header_line  = "\t".join(display_cols)
    print(header_line)
    print("-" * len(header_line))
    for row in filtered:
        print("\t".join(row.get(c, "") for c in display_cols))

    if args.output:
        out = Path(args.output)
        write_tsv(out, metadata, headers, filtered)
        print(f"\nSaved -> {out}")


if __name__ == "__main__":
    main()
