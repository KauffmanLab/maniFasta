#!/usr/bin/env python3
"""
mod_IV_select_representatives.py
─────────────────────────
Dereplicate HOMD genome assemblies to ONE representative per taxonomic group.

Given:
  * a (body-site-filtered) HOMD taxonomy TSV   -> rank value per HMT ID
  * GCA_ID_info.txt                            -> GCA  -> HMT mapping
  * an assembly-stats TSV  (gca, contigs, size)-> quality metrics per GCA

it groups every candidate assembly by the requested taxonomic rank
(species | genus | family) and picks the best assembly in each group:

    1. fewest contigs            (most contiguous assembly)
    2. tie -> largest size       (most bases)
    3. tie -> smallest GCA       (deterministic, reproducible)

Assemblies with no stats are ranked last within their group and only kept
if they are the sole member. stdlib only (matches the rest of maniFasta).
"""

import sys
import csv
import re
import argparse
from pathlib import Path


def info(m): print(f"[INFO] {m}", file=sys.stderr)
def warn(m): print(f"[WARN] {m}", file=sys.stderr)
def die(m):
    print(f"ERROR: {m}", file=sys.stderr)
    raise SystemExit(1)


# eHOMD changed HMT-ID zero-padding between releases (HMT-047 -> HMT-0047) while
# PROKKA/V11.03/GCA_ID_info.txt kept the 3-digit form. Canonicalise both sides of
# the taxonomy <-> GCA_ID_info join so either padding works.
_HMT_RE = re.compile(r"^HMT-?0*([0-9]+)$", re.IGNORECASE)


def norm_hmt(hmt: str) -> str:
    hmt = (hmt or "").strip()
    m = _HMT_RE.match(hmt)
    return f"HMT-{int(m.group(1)):04d}" if m else hmt


# ── taxonomy parsing ──────────────────────────────────────────────────────────

def _read_homd_tax(tax_file: Path):
    """Yield (header, rows) tolerating the leading 'HOMD.org Taxon Data::' banner."""
    header, rows = None, []
    with tax_file.open(encoding="utf-8") as fh:
        for cols in csv.reader(fh, delimiter="\t"):
            if not cols or cols[0].startswith("HOMD.org Taxon Data::"):
                continue
            if header is None:
                header = [c.strip() for c in cols]
                continue
            rows.append(cols)
    if header is None:
        die(f"no header found in taxonomy file: {tax_file}")
    return header, rows


def _col_index(header, *candidates):
    low = [h.lower() for h in header]
    for cand in candidates:
        if cand in low:
            return low.index(cand)
    return None


def build_hmt_rank_map(tax_file: Path, rank: str):
    """
    Return hmt_id -> grouping_key for the requested rank.

    For 'species' the key is composited with genus ("Genus species") so two
    different genera that happen to share an epithet are never merged.
    Blank rank values are keyed per-HMT so unclassified taxa stay singletons.
    """
    header, rows = _read_homd_tax(tax_file)

    idx_hmt = _col_index(header, "hmt-id", "hmt_id", "hmt id")
    if idx_hmt is None:
        idx_hmt = next((i for i, h in enumerate(header) if "hmt" in h.lower()), None)
    if idx_hmt is None:
        die(f"no HMT column in taxonomy header: {header}")

    idx_rank = _col_index(header, rank.lower())
    if idx_rank is None:
        die(f"taxonomy file has no '{rank}' column. Columns present: {header}\n"
            f"       (for 'family' you may need to augment the table with lineage first)")

    idx_genus = _col_index(header, "genus") if rank == "species" else None

    hmt_rank = {}
    for cols in rows:
        if len(cols) <= max(i for i in (idx_hmt, idx_rank, idx_genus) if i is not None):
            continue
        hmt = norm_hmt(cols[idx_hmt].strip())
        if not hmt:
            continue
        val = cols[idx_rank].strip()
        if not val:
            hmt_rank[hmt] = f"__unclassified__{hmt}"
            continue
        if idx_genus is not None:
            genus = cols[idx_genus].strip()
            val = f"{genus} {val}".strip()
        hmt_rank[hmt] = val
    if not hmt_rank:
        die(f"no HMT->{rank} pairs parsed from {tax_file}")
    return hmt_rank


# ── GCA -> HMT (mirrors the parser already in mod_IV_runFetch_HOMD_Proteomes.sh) ──────

def build_gca_hmt(gca_info: Path):
    pairs, seen = [], set()
    with gca_info.open(encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("GCA-ID"):
                continue
            cols = line.split()
            if len(cols) < 2:
                continue
            gca, hmt = cols[0].strip(), norm_hmt(cols[1].strip())
            if gca and gca not in seen:
                pairs.append((gca, hmt))
                seen.add(gca)
    if not pairs:
        die(f"no GCA->HMT pairs parsed from {gca_info}")
    return pairs


# ── stats parsing ─────────────────────────────────────────────────────────────

def build_stats(stats_file: Path):
    """Read a TSV with header columns gca, contigs, size -> gca -> (contigs, size)."""
    stats = {}
    with stats_file.open(encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fields = {f.lower(): f for f in (reader.fieldnames or [])}
        for need in ("gca", "contigs", "size"):
            if need not in fields:
                die(f"stats file missing '{need}' column; has {reader.fieldnames}")
        for row in reader:
            gca = row[fields["gca"]].strip()
            try:
                stats[gca] = (int(row[fields["contigs"]]), int(row[fields["size"]]))
            except (ValueError, TypeError):
                continue
    return stats


# ── selection ─────────────────────────────────────────────────────────────────

def select(hmt_rank, gca_pairs, stats):
    BIG = float("inf")
    groups, skipped = {}, 0
    for gca, hmt in gca_pairs:
        key = hmt_rank.get(hmt)
        if key is None:                       # HMT filtered out by body site
            skipped += 1
            continue
        groups.setdefault(key, []).append((gca, hmt))
    if skipped:
        info(f"{skipped} assemblies skipped (HMT not in filtered taxonomy)")

    reps, no_stats = [], 0
    for key, members in sorted(groups.items()):
        # fewest contigs, then largest size, then GCA name
        members.sort(key=lambda gh: (*stats.get(gh[0], (BIG, -1))[:1],
                                     -stats.get(gh[0], (BIG, -1))[1], gh[0]))
        best_gca, best_hmt = members[0]
        c, s = stats.get(best_gca, (None, None))
        if c is None:
            no_stats += 1
            warn(f"group '{key}': no stats for any member; kept '{best_gca}' by name")
        reps.append({
            "gca": best_gca, "hmt_id": best_hmt, "rank_value": key,
            "contigs": "" if c is None else c,
            "size": "" if s is None else s,
            "n_candidates": len(members),
        })
    info(f"groups: {len(groups)}  representatives: {len(reps)}")
    if no_stats:
        warn(f"{no_stats} group(s) had no usable stats")
    return reps


# ── CLI ───────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tax-file", required=True, type=Path,
                    help="(filtered) HOMD taxonomy TSV from mod_IV_homd_filter.py")
    ap.add_argument("--gca-info", required=True, type=Path,
                    help="GCA_ID_info.txt (GCA -> HMT)")
    ap.add_argument("--stats", required=True, type=Path,
                    help="TSV with columns: gca, contigs, size")
    ap.add_argument("--rank", required=True, choices=["species", "genus", "family"],
                    help="taxonomic level to dereplicate at")
    ap.add_argument("--output", required=True, type=Path,
                    help="output TSV report of selected representatives")
    ap.add_argument("--plan", type=Path,
                    help="optional '<GCA> <HMT>' plan file for the downloader")
    args = ap.parse_args()

    for p in (args.tax_file, args.gca_info, args.stats):
        if not p.exists():
            die(f"file not found: {p}")

    hmt_rank = build_hmt_rank_map(args.tax_file, args.rank)
    info(f"HMT IDs carrying a '{args.rank}' value: {len(hmt_rank)}")
    gca_pairs = build_gca_hmt(args.gca_info)
    info(f"candidate GCA->HMT pairs: {len(gca_pairs)}")
    stats = build_stats(args.stats)
    info(f"assemblies with stats: {len(stats)}")

    reps = select(hmt_rank, gca_pairs, stats)
    if not reps:
        die("no representatives selected")

    cols = ["gca", "hmt_id", "rank_value", "contigs", "size", "n_candidates"]
    with args.output.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(reps)
    info(f"wrote representative report -> {args.output}")

    if args.plan:
        with args.plan.open("w", encoding="utf-8") as fh:
            for r in reps:
                fh.write(f"{r['gca']} {r['hmt_id']}\n")
        info(f"wrote download plan -> {args.plan}")


if __name__ == "__main__":
    main()
