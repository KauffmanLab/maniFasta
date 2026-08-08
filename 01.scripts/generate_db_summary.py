#!/usr/bin/env python3
"""
generate_db_summary.py
──────────────────────
Generate a methods-style summary of a finished maniFasta build by READING the
run directory's existing artifacts (it computes nothing new and fetches nothing).

It addresses the standard search-database reporting requirements:
  • number of sequences in the database
  • download/build date for publicly sourced sequences
  • origin of the sequences (per source)
  • how taxa were selected or excluded
  • how redundancy was handled
  • how functional and taxonomic annotations were assigned
  • how contaminant sequences were incorporated

Output, written into the run directory alongside build_summary.tsv:
  <LABEL>.database_methods.md    methods paragraph + per-source table

Artifacts read (all optional; missing ones degrade gracefully):
  <LABEL>.fasta                                  final DB FASTA (authoritative count)
  build_summary.tsv                              build-phase summary
  <LABEL>.manifest[.lineage].tsv                 per-protein source + taxonomy
  run_logs/source_plan.normalized.tsv            source registry (normalized)
  run_logs/source_plan.supported_collections.env homd_sites/homd_rank/crap_set/...
  run_data/homd/HOMD_representatives.tsv          dereplication provenance
  run_data/homd/HOMD_assembly_stats.tsv           candidate assembly stats
  run_data/homd/HOMD_prokka_selection.tsv         body-site candidate list

Taxonomy reporting rule: build_db.py writes a populated lineage_taxid column
into EVERY manifest, so a taxid count can never be used as evidence that
PHASE 2.5 ran. Ranked-lineage claims are made only when the manifest carries
tax_<rank> columns, which only add_lineage_to_metadata.py appends; otherwise
the methods text states plainly that enrichment was not performed. Placeholder
taxids (blank, non-numeric, or NCBI 'unidentified' 32644) are excluded from
represented-taxon counts, matching add_lineage_to_metadata.py's is_resolvable().

stdlib only.
"""

from __future__ import annotations

import argparse
import csv
import re
import shlex
import sys
from datetime import datetime
from pathlib import Path

csv.field_size_limit(min(sys.maxsize, 2**31 - 1))

# The NCBI "unidentified" taxid build_db.py assigns when nothing resolved.
# Keep in sync with add_lineage_to_metadata.py's UNIDENTIFIED_TAXID.
UNIDENTIFIED_TAXID = "32644"

# Default of add_lineage_to_metadata.py's --unclassified-label.
DEFAULT_UNCLASSIFIED_LABEL = "Unclassified"


def is_real_taxid(tid: str) -> bool:
    """Same rule as add_lineage_to_metadata.py's is_resolvable(): a taxid only
    counts as a represented taxon if it is numeric and is not the NCBI
    'unidentified' placeholder. Blank, 'NA' and 'Unclassified' all fail here."""
    tid = (tid or "").strip()
    return tid.isdigit() and tid != UNIDENTIFIED_TAXID


def info(m): print(f"[INFO] {m}", file=sys.stderr)
def warn(m): print(f"[WARN] {m}", file=sys.stderr)
def die(m):
    print(f"ERROR: {m}", file=sys.stderr)
    raise SystemExit(1)


# ── artifact resolution ───────────────────────────────────────────────────────

def first_existing(*candidates: Path):
    for c in candidates:
        if c and c.is_file():
            return c
    return None


def infer_label(run_dir: Path, explicit: str | None) -> str:
    if explicit:
        return explicit
    # A <LABEL>.fasta directly in the run dir is the strongest signal.
    fastas = [p for p in run_dir.glob("*.fasta") if p.is_file()]
    fastas = [p for p in fastas if not p.name.startswith("crap")]
    if len(fastas) == 1:
        return fastas[0].stem
    # Fall back to the run-dir name: <YYYYMMDD_HHMMSS>_<LABEL>
    m = re.match(r"^\d{8}_\d{6}_(.+)$", run_dir.name)
    if m:
        return m.group(1)
    return run_dir.name


def parse_build_date(run_dir: Path) -> str | None:
    """Build/download date from the run-dir timestamp (YYYYMMDD_HHMMSS_...)."""
    m = re.match(r"^(\d{8})_(\d{6})_", run_dir.name)
    if not m:
        return None
    try:
        dt = datetime.strptime(m.group(1) + m.group(2), "%Y%m%d%H%M%S")
        return dt.strftime("%Y-%m-%d")
    except ValueError:
        return None


# ── small readers ─────────────────────────────────────────────────────────────

def count_fasta(path: Path) -> int | None:
    if not path or not path.is_file():
        return None
    n = 0
    with path.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                n += 1
    return n


def read_supported_env(path: Path) -> dict:
    """Parse 'export key=value' lines from the supported-collections env."""
    env = {}
    if not path or not path.is_file():
        return env
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        line = line[len("export "):].strip() if line.startswith("export ") else line
        if "=" not in line:
            continue
        k, v = line.split("=", 1)
        try:
            parts = shlex.split(v)
            v = parts[0] if parts else ""
        except ValueError:
            v = v.strip().strip("'\"")
        env[k.strip()] = v
    return env


def parse_options(raw: str) -> dict:
    opts = {}
    for item in (raw or "").split(";"):
        item = item.strip()
        if "=" in item:
            k, v = item.split("=", 1)
            opts[k.strip()] = v.strip()
    return opts


def read_normalized_plan(path: Path) -> list[dict]:
    """Enabled, non-disabled source rows from the normalized plan."""
    rows = []
    if not path or not path.is_file():
        return rows
    with path.open(newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if (r.get("enabled") or "").strip().lower() != "true":
                continue
            if (r.get("status") or "").strip().lower() == "disabled":
                continue
            rows.append(r)
    return rows


def aggregate_manifest(path: Path,
                       unclassified_label: str = DEFAULT_UNCLASSIFIED_LABEL) -> dict:
    """
    Stream the manifest. Returns per-source counts, a placeholder-free
    distinct-taxid count, distinct counts for the tax_<rank> facets, and — only
    when the file actually carries lineage-enrichment columns — how many rows
    received a ranked lineage.

    Two different things are counted here and they must not be conflated:

      * A taxid column (lineage_taxid / ncbi_taxid) is written by build_db.py
        into EVERY manifest, enriched or not: write_manifest_row() always
        populates lineage_taxid, defaulting to the row's ncbi_taxid. Its
        presence therefore says nothing about whether PHASE 2.5 ran.
      * tax_<rank> columns are appended only by add_lineage_to_metadata.py,
        so their presence is the authoritative signal that ranked lineage
        enrichment happened. has_lineage_cols carries that signal.
    """
    out = {
        "total": 0, "per_source": {}, "taxids": 0, "ranks": {},
        "placeholder_rows": 0, "has_lineage_cols": False,
        "lineage_resolved": 0, "lineage_unclassified": 0,
    }
    if not path or not path.is_file():
        return out
    unclass_lc = (unclassified_label or DEFAULT_UNCLASSIFIED_LABEL).strip().lower()
    with path.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        cols = {c.lower(): c for c in (reader.fieldnames or [])}
        src_col = next((cols[c] for c in ("source", "source_label") if c in cols), None)
        tax_col = next((cols[c] for c in ("lineage_taxid", "ncbi_taxid", "taxid") if c in cols), None)
        lineage_cols = {c: cols[c] for c in cols if c.startswith("tax_")}
        # own_rank / own_name / lineage_path describe the row itself; only the
        # remaining tax_ columns are ranked facets worth counting distinctly.
        facet_cols = {c: v for c, v in lineage_cols.items()
                      if c not in ("tax_own_rank", "tax_own_name", "tax_lineage_path")}
        own_rank_col = lineage_cols.get("tax_own_rank")
        path_col = lineage_cols.get("tax_lineage_path")
        has_lineage = bool(lineage_cols)
        out["has_lineage_cols"] = has_lineage

        def row_resolved(row) -> bool:
            """add_lineage_to_metadata.py's unclassified_row() writes an empty
            own_rank (a resolved taxon always has at least 'no rank'), so a
            non-empty own_rank marks a genuinely resolved row regardless of what
            --unclassified-label was set to. lineage_path is the fallback for
            enriched files written with a different column set."""
            if own_rank_col:
                return bool((row.get(own_rank_col) or "").strip())
            if path_col:
                v = (row.get(path_col) or "").strip()
                return bool(v) and v.lower() != unclass_lc
            return False

        taxids: set[str] = set()
        rank_sets: dict[str, set] = {r: set() for r in facet_cols}
        per_source: dict[str, int] = {}
        total = 0
        placeholders = 0
        resolved = 0
        for row in reader:
            total += 1
            if src_col:
                s = (row.get(src_col) or "").strip() or "NA"
                per_source[s] = per_source.get(s, 0) + 1
            if tax_col:
                t = (row.get(tax_col) or "").strip()
                if is_real_taxid(t):
                    taxids.add(t)
                else:
                    placeholders += 1
            for rk, col in facet_cols.items():
                v = (row.get(col) or "").strip()
                if v and v.lower() not in ("na", unclass_lc):
                    rank_sets[rk].add(v)
            if has_lineage and row_resolved(row):
                resolved += 1
        out["total"] = total
        out["per_source"] = per_source
        out["taxids"] = len(taxids)
        out["ranks"] = {rk.replace("tax_", "", 1): len(s) for rk, s in rank_sets.items()}
        out["placeholder_rows"] = placeholders
        out["lineage_resolved"] = resolved if has_lineage else 0
        out["lineage_unclassified"] = (total - resolved) if has_lineage else 0
    return out


def read_homd_reps(path: Path) -> dict | None:
    if not path or not path.is_file():
        return None
    n_reps = 0
    cand_sum = 0
    contigs_one = 0
    with path.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        cols = {c.lower(): c for c in (reader.fieldnames or [])}
        ncand = cols.get("n_candidates")
        ccol = cols.get("contigs")
        for row in reader:
            n_reps += 1
            if ncand:
                try:
                    cand_sum += int(row[ncand])
                except (ValueError, KeyError):
                    pass
            if ccol:
                try:
                    if int(row[ccol]) == 1:
                        contigs_one += 1
                except (ValueError, KeyError):
                    pass
    return {"n_representatives": n_reps, "n_candidates": cand_sum, "n_single_contig": contigs_one}


def count_tsv_rows(path: Path) -> int | None:
    if not path or not path.is_file():
        return None
    with path.open(encoding="utf-8", errors="replace") as fh:
        n = sum(1 for _ in fh)
    return max(n - 1, 0)  # minus header


# ── assemble facts ────────────────────────────────────────────────────────────

def gather(run_dir: Path, label: str, tmp_dir: Path | None = None,
           unclassified_label: str = DEFAULT_UNCLASSIFIED_LABEL) -> dict:
    logs = run_dir / "run_logs"
    homd = run_dir / "run_data" / "homd"
    tmp = tmp_dir if tmp_dir else (run_dir / "tmp")

    fasta = first_existing(run_dir / f"{label}.fasta")
    manifest = first_existing(
        run_dir / f"{label}.manifest.lineage.tsv",
        run_dir / f"{label}.manifest.tsv",
    )
    plan = first_existing(
        logs / "source_plan.normalized.tsv",
        tmp / "source_plan.normalized.tsv",
    )
    supp_env = first_existing(
        logs / "source_plan.supported_collections.env",
        tmp / "source_plan.supported_collections.env",
    )
    reps = first_existing(homd / "HOMD_representatives.tsv")
    stats = first_existing(homd / "HOMD_assembly_stats.tsv")
    selection = first_existing(homd / "HOMD_prokka_selection.tsv")

    facts = {
        "pipeline": "maniFasta",
        "run_label": label,
        "run_dir": str(run_dir),
        "build_date": parse_build_date(run_dir),
        "total_sequences": count_fasta(fasta),
        "artifacts": {
            k: (str(v) if v else None) for k, v in {
                "final_fasta": fasta, "manifest": manifest, "plan": plan,
                "supported_env": supp_env, "homd_representatives": reps,
                "homd_assembly_stats": stats, "homd_selection": selection,
            }.items()
        },
    }

    env = read_supported_env(supp_env)
    manifest_agg = aggregate_manifest(manifest, unclassified_label)
    if manifest is not None and manifest.name.endswith(".manifest.lineage.tsv") \
            and not manifest_agg["has_lineage_cols"]:
        warn(f"{manifest.name} carries no tax_* columns; treating this build as "
             f"not lineage-enriched (PHASE 2.5 may have failed part-way).")
    plan_rows = read_normalized_plan(plan)

    # ── sources ──
    sources = []
    for r in plan_rows:
        sl = (r.get("source_label") or "").strip()
        opts = parse_options(r.get("options") or "")
        count = manifest_agg["per_source"].get(sl)
        sources.append({
            "source_label": sl,
            "source_group": (r.get("source_group") or "").strip(),
            "module": (r.get("module") or "").strip(),
            "collection": (r.get("collection") or "").strip(),
            "version": (r.get("version") or "").strip(),
            "citation": (r.get("citation") or "").strip(),
            "notes": (r.get("notes") or "").strip(),
            "options": opts,
            "sequence_count": count,
        })
    facts["sources"] = sources
    facts["manifest_total"] = manifest_agg["total"] or None
    facts["distinct_taxids"] = manifest_agg["taxids"] or None
    facts["taxonomic_ranks"] = {k: v for k, v in manifest_agg["ranks"].items() if v} or None
    facts["placeholder_taxid_rows"] = manifest_agg["placeholder_rows"] or None

    # Whether ranked lineage enrichment actually ran, and how far it got. This
    # is the ONLY thing the methods text may branch on when claiming lineage
    # assignment -- never the taxid count, which every manifest has.
    facts["lineage"] = ({
        "manifest": str(manifest) if manifest else None,
        "resolved": manifest_agg["lineage_resolved"],
        "unclassified": manifest_agg["lineage_unclassified"],
        "total": manifest_agg["total"],
        "label": unclassified_label,
    } if manifest_agg["has_lineage_cols"] else None)

    # ── HOMD dereplication ──
    rep_info = read_homd_reps(reps)
    if env.get("use_homd", "").lower() == "true" or rep_info:
        candidates = count_tsv_rows(selection)
        with_stats = count_tsv_rows(stats)
        facts["homd"] = {
            "version": env.get("homd_genomic_refseq_version") or None,
            "sites": env.get("homd_sites") or None,
            "rank": env.get("homd_rank") or None,
            "candidates_after_site_filter": candidates,
            "candidates_with_stats": with_stats,
            "representatives": (rep_info or {}).get("n_representatives"),
            "candidates_collapsed": (rep_info or {}).get("n_candidates"),
            "single_contig_reps": (rep_info or {}).get("n_single_contig"),
        }
    else:
        facts["homd"] = None

    # ── cRAP contaminants ──
    crap_set = env.get("crap_set", "none")
    if crap_set and crap_set != "none":
        crap_doi = None
        crap_label = None
        for s in sources:
            if (s["collection"] or "").lower() == "crap":
                crap_doi = s["citation"] or None
                crap_label = s["source_label"]
        facts["crap"] = {
            "set": crap_set,
            "doi": crap_doi,
            "source_label": crap_label,
            "sequence_count": (manifest_agg["per_source"].get(crap_label) if crap_label else None),
        }
    else:
        facts["crap"] = None

    return facts


# ── render markdown ───────────────────────────────────────────────────────────

def _fmt_int(n):
    return f"{n:,}" if isinstance(n, int) else "[REVIEW: not available]"


_RANK_PLURAL = {
    "domain": "domains", "kingdom": "kingdoms", "phylum": "phyla",
    "class": "classes", "order": "orders", "family": "families",
    "genus": "genera", "species": "species",
}


def _plural_rank(rank: str) -> str:
    return _RANK_PLURAL.get(rank, rank + "s")


def render_markdown(f: dict) -> str:
    L = []
    label = f["run_label"]
    date = f["build_date"] or "[REVIEW: build date unavailable]"
    total = f["total_sequences"]

    L.append(f"# Search database composition — {label}\n")
    L.append("@@SUBTITLE@@\n")

    homd = f.get("homd")
    crap = f.get("crap")
    grp = sorted({s["source_group"] for s in f["sources"] if s["source_group"]})
    n_src = len(f["sources"])

    # ── One flowing methods block (composition -> taxon selection -> redundancy
    #    -> annotation -> contaminants). No sub-headers; reads as methods prose. ──
    L.append("## Methods summary\n")

    # Paragraph 1: composition, taxonomic breadth, and how taxa were selected.
    p1 = []
    p1.append(
        f"The {label} search database comprised {_fmt_int(total)} protein sequence(s) "
        f"assembled with the maniFasta pipeline"
        + (f" on {date}" if f["build_date"] else "")
        + f". Sequences were drawn from {n_src} curated source"
        + ("s" if n_src != 1 else "")
        + (f" spanning {len(grp)} categories ({', '.join(grp)})" if grp else "")
        + "."
    )
    if f["build_date"]:
        p1.append("All publicly sourced sequences were retrieved on that build date.")
    if f.get("distinct_taxids"):
        ranks = f.get("taxonomic_ranks") or {}
        bits = [f"{v} {_plural_rank(k)}" for k, v in (("phylum", ranks.get("phylum")),
                                                      ("family", ranks.get("family")),
                                                      ("genus", ranks.get("genus")),
                                                      ("species", ranks.get("species"))) if v]
        extra = (" spanning " + ", ".join(bits)) if bits else ""
        p1.append(
            f"Across all sources the database represents {_fmt_int(f['distinct_taxids'])} "
            f"distinct NCBI taxa{extra}."
        )
        if f.get("placeholder_taxid_rows"):
            p1.append(
                f"A further {_fmt_int(f['placeholder_taxid_rows'])} sequence(s) carry no "
                f"resolved taxon ID (NCBI 'unidentified', taxid {UNIDENTIFIED_TAXID}) and are "
                f"counted as unclassified rather than as distinct represented taxa."
            )
    if homd and homd.get("sites"):
        p1.append(
            f"HOMD genomes were restricted to the {homd['sites']} body site(s) using the HOMD "
            f"taxonomy table prior to any redundancy reduction, while non-HOMD sources were curated "
            f"inclusion lists whose taxa were selected by the criteria noted per source below "
            f"(e.g., cultivation or pathogenicity evidence, tissue tropism, or literature curation) "
            f"rather than by automated taxonomic sweeps."
        )
    elif homd:
        p1.append(
            "HOMD genomes were included across all body sites with no site filter applied, while "
            "non-HOMD sources were curated inclusion lists whose taxa were selected by the criteria "
            "noted per source below rather than by automated taxonomic sweeps."
        )
    else:
        p1.append(
            "Taxa were selected from curated inclusion lists by the criteria noted per source below "
            "rather than by automated taxonomic sweeps."
        )

    # Paragraph 2: redundancy handling.
    p2 = []
    if homd and homd.get("rank"):
        sentence = (
            f"Within the site-filtered HOMD set, genome assemblies were dereplicated to one "
            f"representative proteome per {homd['rank']} "
            f"({_fmt_int(homd.get('representatives'))} representatives selected from "
            f"{_fmt_int(homd.get('candidates_collapsed') or homd.get('candidates_after_site_filter'))} "
            f"candidate assemblies), choosing for each group the assembly with the fewest contigs and "
            f"breaking ties by the largest genome size"
        )
        if homd.get("single_contig_reps") is not None and homd.get("representatives"):
            sentence += (
                f"; {homd['single_contig_reps']} of {homd['representatives']} representatives were "
                f"single-contig (closed or near-closed) assemblies."
            )
        else:
            sentence += "."
        p2.append(sentence)
    elif homd:
        p2.append(
            "HOMD genomes were included without taxonomic dereplication (one proteome per assembly "
            "passing the site filter)."
        )
    else:
        p2.append("No taxonomic dereplication was applied; all sequences from each source were retained.")

    # Paragraph 3: functional/taxonomic annotation and contaminants.
    p3 = []
    if homd:
        p3.append(
            "HOMD proteins carried PROKKA functional annotations as distributed by HOMD "
            f"(release {homd.get('version') or '[REVIEW: version]'}), with product descriptions and "
            "locus tags retained in the FASTA headers, while sequences contributed as curated FASTAs "
            "or accession lists retained their source-provided functional annotations."
        )
    else:
        p3.append(
            "Sequences contributed as curated FASTAs or accession lists retained their source-provided "
            "functional annotations."
        )
    lin = f.get("lineage")
    if lin:
        sentence = (
            "Taxonomic annotation was assigned by resolving each sequence's manifest taxon ID to a "
            "full ranked lineage (NCBI Taxonomy) during the lineage-enrichment step"
        )
        if lin["total"]:
            sentence += (
                f"; {_fmt_int(lin['resolved'])} of {_fmt_int(lin['total'])} sequence(s) received "
                f"ranked lineage information and {_fmt_int(lin['unclassified'])} remained "
                f"{lin['label']}"
            )
        p3.append(sentence + ".")
    elif f["artifacts"].get("manifest"):
        p3.append(
            "NCBI taxon IDs were retained in the manifest; ranked lineage enrichment was not "
            "performed for this build."
        )
    else:
        p3.append("[REVIEW: describe how taxonomic lineage was assigned (lineage-enrichment step).]")
    if crap:
        cnt = _fmt_int(crap["sequence_count"]) if crap.get("sequence_count") is not None else "[REVIEW: count]"
        p3.append(
            f"Common contaminant proteins were incorporated from the cRAP collection (set: {crap['set']}"
            + (f"; {crap['doi']}" if crap.get("doi") else "")
            + f"), contributing {cnt} sequence(s) merged into the final database alongside the "
            "biological sources."
        )
    else:
        p3.append("No contaminant (cRAP) sequences were included in this build.")

    L.append(" ".join(p1) + "\n")
    L.append(" ".join(p2) + "\n")
    L.append(" ".join(p3) + "\n")

    # ── Origin of sequences (table + per-source provenance), kept at the bottom ──
    L.append("## Origin of sequences\n")
    L.append("| Source | Category | Module / collection | Version | Sequences | Citation |")
    L.append("|---|---|---|---|---|---|")
    for s in f["sources"]:
        modcol = s["collection"] or s["module"]
        cnt = _fmt_int(s["sequence_count"]) if s["sequence_count"] is not None else "—"
        L.append(
            f"| {s['source_label'] or '—'} | {s['source_group'] or '—'} | {modcol or '—'} "
            f"| {s['version'] or '—'} | {cnt} | {s['citation'] or '—'} |"
        )
    L.append("")
    noted = [s for s in f["sources"] if s["notes"]]
    if noted:
        L.append("Source-specific provenance:\n")
        for s in noted:
            L.append(f"- **{s['source_label']}** — {s['notes']}")
        L.append("")

    md = "\n".join(L)
    if "[REVIEW:" in md:
        subtitle = "*Auto-generated from build artifacts. Resolve the [REVIEW: …] items below before submission.*\n"
    else:
        subtitle = "*Auto-generated from build artifacts.*\n"
    return md.replace("@@SUBTITLE@@\n", subtitle)


# ── CLI ───────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-dir", required=True, type=Path, help="finished build/run directory")
    ap.add_argument("--run-label", default=None, help="run label (inferred if omitted)")
    ap.add_argument("--out", default=None, type=Path, help="markdown output (default <run-dir>/<LABEL>.database_methods.md)")
    ap.add_argument("--tmp-dir", default=None, type=Path,
                    help="run TMPDIR holding the in-progress plan/env (default <run-dir>/tmp); "
                         "used when invoked from the pipeline before artifacts are packaged into run_logs/")
    ap.add_argument("--unclassified-label", default=DEFAULT_UNCLASSIFIED_LABEL,
                    help="Label add_lineage_to_metadata.py wrote for unresolved rows "
                         f"(default {DEFAULT_UNCLASSIFIED_LABEL}); only change it if that "
                         "step was run with a non-default --unclassified-label.")
    ap.add_argument("--print", dest="to_stdout", action="store_true", help="also echo the markdown to stdout")
    args = ap.parse_args()

    run_dir = args.run_dir
    if not run_dir.is_dir():
        die(f"run directory not found: {run_dir}")

    label = infer_label(run_dir, args.run_label)
    info(f"run label: {label}")

    facts = gather(run_dir, label, args.tmp_dir, args.unclassified_label)
    lin = facts.get("lineage")
    if lin:
        info(f"lineage enrichment detected: {lin['resolved']} resolved / "
             f"{lin['unclassified']} {lin['label']} of {lin['total']} rows")
    else:
        info("no lineage-enriched manifest; methods text will state that ranked "
             "lineage enrichment was not performed")

    md_path = args.out or (run_dir / f"{label}.database_methods.md")

    md = render_markdown(facts)
    md_path.write_text(md, encoding="utf-8")

    info(f"methods summary -> {md_path}")
    if facts["total_sequences"] is None:
        warn("final FASTA not found; total sequence count is unavailable (see [REVIEW] notes)")
    if args.to_stdout:
        print(md)


if __name__ == "__main__":
    main()
