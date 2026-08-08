#!/usr/bin/env python3
"""
source_planner.py
─────────────────
Source planner for the maniFasta source-registry model.

Reads 000.maniFasta.source_registry.tsv, validates enabled source rows,
resolves registry file fields, and emits one run-specific normalized source plan:

  source_plan.normalized.tsv
     Fully interpreted source plan consumed directly by runBuild_db.sh.

Path model
──────────
The registry may include a source_dir column. When source_dir is present for a
row, file fields such as input_list, fasta_file, and metadata_file are resolved
as:

  <mainDIR>/<source_dir>/<file field>

Examples:

  source_dir              input_list
  99.prepackaged_inputs   AllOralsDB_v1.0.microeuks.tsv
  00.setup                my_user_uploaded_species.tsv

Rows without source_dir keep the older behavior: relative paths are resolved
against --setup-dir, which is usually ${registry_path_base:-${setupDIR}}.
Absolute paths and shell-style paths beginning with $ or ~ are preserved.
"""

from __future__ import annotations

import argparse
import csv
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List

TRUE = {"true", "t", "yes", "y", "1", "on"}
FALSE = {"false", "f", "no", "n", "0", "off", ""}

# Canonical module identifiers (case-sensitive). These name the input LANE only;
# the data backend is the orthogonal fetch_source option (ncbi/uniprot/uniparc/
# pdb), which is why e.g. mod_II is not tied to NCBI.
#   mod_I   = proteome downloads        (was ncbi_proteome)
#   mod_II  = protein-accession fetch   (was ncbi_protein)
#   mod_III = user-supplied FASTA       (was user_fasta)
#   mod_IV  = supported collections     (was supported_collection)
VALID_MODULES = {"mod_I", "mod_II", "mod_III", "mod_IV"}

# Legacy code-names are HARD-REJECTED (no longer aliased). This map exists only
# to emit a helpful "rename to X" error so old registries fail loudly and
# pointably rather than with a generic "unsupported module".
LEGACY_MODULE_RENAMES = {
    "ncbi_proteome":        "mod_I",
    "ncbi_protein":         "mod_II",
    "user_fasta":           "mod_III",
    "supported_collection": "mod_IV",
    # old convenience shorthands also retired:
    "proteome":             "mod_I",
    "protein":              "mod_II",
}
SUPPORTED_COLLECTIONS = {"human_uniprot", "HOMD", "cRAP", "crap", "homd"}
# Fetch backends understood by the per-row/per-source fetch_source option.
VALID_FETCH_SOURCES = {"ncbi", "uniprot", "uniparc", "pdb"}
# Representative-selection rank for the HOMD collection (dereplication level).
# Absent/empty disables dereplication; otherwise ONE representative proteome is
# kept per group at this rank (fewest contigs, ties -> largest size), enforced
# downstream by mod_IV_select_representatives.py. Body-site filtering ('sites=') is
# always applied first, then dereplication runs on the survivors.
VALID_HOMD_RANKS = {"species", "genus", "family"}


def die(msg: str, code: int = 1) -> None:
    print(f"ERROR: {msg}", file=sys.stderr)
    raise SystemExit(code)


def warn(msg: str) -> None:
    print(f"[WARN] {msg}", file=sys.stderr)


def info(msg: str) -> None:
    print(f"[INFO] {msg}", file=sys.stderr)


def parse_bool(value: str, *, field: str) -> bool:
    v = (value or "").strip().lower()
    if v in TRUE:
        return True
    if v in FALSE:
        return False
    die(f"invalid boolean in {field}: {value!r}")


def parse_options(value: str) -> Dict[str, str]:
    """Parse key=value;key=value options."""
    opts: Dict[str, str] = {}
    value = (value or "").strip()
    if not value:
        return opts
    for item in value.split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" not in item:
            die(f"bad option {item!r}; options must be key=value separated by semicolons")
        k, v = item.split("=", 1)
        k = k.strip()
        v = v.strip()
        if not k:
            die(f"bad option {item!r}; empty option key")
        opts[k] = v
    return opts


def resolve_path(
    path_value: str,
    *,
    main_dir: Path,
    setup_dir: Path,
    source_dir: str = "",
) -> str:
    """
    Resolve a registry file field to an absolute path.

    New preferred behavior:
      source_dir=99.prepackaged_inputs + input_list=file.tsv
        -> <mainDIR>/99.prepackaged_inputs/file.tsv

      source_dir=00.setup + input_list=file.tsv
        -> <mainDIR>/00.setup/file.tsv

    Legacy behavior:
      If source_dir is blank, relative paths are resolved against setup_dir.
      This keeps older registries with paths like ../99.prepackaged_inputs/file.tsv
      working while registries are migrated.

    Empty fields stay empty. Absolute paths and shell-style paths beginning with
    $ or ~ are preserved for advanced/manual use.
    """
    p = (path_value or "").strip()
    sd = (source_dir or "").strip()

    if not p:
        return ""

    # Permit shell-style variables to survive untouched for advanced users.
    if p.startswith("$") or p.startswith("~"):
        return p

    pp = Path(p)
    if pp.is_absolute():
        return str(pp)

    if sd:
        source_root = Path(sd)
        if source_root.is_absolute() or ".." in source_root.parts:
            die(
                f"invalid source_dir={sd!r}; use a project-relative folder "
                "like 00.setup or 99.prepackaged_inputs"
            )
        if ".." in pp.parts:
            die(
                f"invalid file path {p!r} with source_dir={sd!r}; "
                "when source_dir is set, file fields should be filenames or "
                "subpaths within that source_dir, not ../ paths"
            )
        return str(main_dir / source_root / pp)

    return str(setup_dir / pp)


@dataclass
class SourceRow:
    source_id: str
    enabled: bool
    source_label: str
    source_group: str
    module: str
    collection: str
    source_dir: str
    input_list: str
    fasta_file: str
    metadata_file: str
    options_raw: str
    options: Dict[str, str]
    version: str
    citation: str
    notes: str

    @property
    def collection_norm(self) -> str:
        c = self.collection.strip()
        if c.lower() == "homd":
            return "HOMD"
        if c.lower() == "crap":
            return "cRAP"
        return c


def read_registry(path: Path, main_dir: Path, setup_dir: Path) -> List[SourceRow]:
    if not path.is_file():
        die(f"source registry not found: {path}")

    rows: List[SourceRow] = []
    with path.open(newline="", encoding="utf-8") as fh:
        data_lines = [line for line in fh if line.strip() and not line.startswith("#")]

    if not data_lines:
        die(f"source registry has no non-comment rows: {path}")

    reader = csv.DictReader(data_lines, delimiter="\t")
    required = [
        "source_id", "enabled", "source_label", "source_group", "module",
        "collection", "input_list", "fasta_file", "metadata_file", "options",
    ]
    missing = [c for c in required if c not in (reader.fieldnames or [])]
    if missing:
        die(f"source registry missing required column(s): {', '.join(missing)}")

    has_source_dir = "source_dir" in (reader.fieldnames or [])
    if not has_source_dir:
        warn(
            "source registry has no source_dir column; using legacy relative-path "
            "behavior against --setup-dir"
        )

    seen_ids = set()
    for n, r in enumerate(reader, start=2):
        source_id = (r.get("source_id") or "").strip()
        if not source_id:
            die(f"row {n}: source_id is required")
        if source_id in seen_ids:
            die(f"duplicate source_id: {source_id}")
        seen_ids.add(source_id)

        enabled = parse_bool(r.get("enabled", ""), field=f"enabled for {source_id}")
        module = (r.get("module") or "").strip()
        # HARD CUT: legacy code-names (ncbi_proteome/ncbi_protein/user_fasta/
        # supported_collection and the old protein/proteome shorthands) are no
        # longer accepted. They map cleanly onto mod_I..mod_IV, so fail loudly
        # with the exact replacement instead of a generic error.
        if module in LEGACY_MODULE_RENAMES:
            die(f"row {source_id}: module={module!r} is a retired code-name; "
                f"rename it to {LEGACY_MODULE_RENAMES[module]!r} in the registry "
                f"(modules are now mod_I..mod_IV; the data backend lives in "
                f"fetch_source, not the module name).")
        if module not in VALID_MODULES:
            die(f"row {source_id}: unsupported module={module!r}; "
                f"choose one of {sorted(VALID_MODULES)}")

        source_dir = (r.get("source_dir") or "").strip()
        opts_raw = r.get("options") or ""
        row = SourceRow(
            source_id=source_id,
            enabled=enabled,
            source_label=(r.get("source_label") or "").strip(),
            source_group=(r.get("source_group") or "").strip(),
            module=module,
            collection=(r.get("collection") or "").strip(),
            source_dir=source_dir,
            input_list=resolve_path(
                r.get("input_list") or "",
                main_dir=main_dir,
                setup_dir=setup_dir,
                source_dir=source_dir,
            ),
            fasta_file=resolve_path(
                r.get("fasta_file") or "",
                main_dir=main_dir,
                setup_dir=setup_dir,
                source_dir=source_dir,
            ),
            metadata_file=resolve_path(
                r.get("metadata_file") or "",
                main_dir=main_dir,
                setup_dir=setup_dir,
                source_dir=source_dir,
            ),
            options_raw=opts_raw,
            options=parse_options(opts_raw),
            version=(r.get("version") or "").strip(),
            citation=(r.get("citation") or "").strip(),
            notes=(r.get("notes") or "").strip(),
        )
        if enabled and not row.source_label:
            die(f"row {source_id}: source_label is required for enabled rows")
        rows.append(row)

    return rows


def validate_rows(rows: Iterable[SourceRow], *, strict_files: bool) -> None:
    rows = list(rows)
    enabled_rows = [row for row in rows if row.enabled]

    # source_label is the downstream provenance key written into manifests and
    # FASTA headers. Duplicate enabled labels would merge distinct sources under
    # one identity even when their source_id values are unique.
    enabled_labels = [row.source_label for row in enabled_rows]
    duplicate_labels = sorted({
        label for label in enabled_labels if enabled_labels.count(label) > 1
    })
    if duplicate_labels:
        die(f"enabled source_label values must be unique: {duplicate_labels}")

    # These supported collections are singleton choices. Multiple enabled rows
    # currently collapse into shared shell variables (last row wins), so reject
    # the registry instead of silently building with only the final selection.
    human_rows = [
        row for row in enabled_rows
        if row.module == "mod_IV" and row.collection_norm == "human_uniprot"
    ]
    homd_rows = [
        row for row in enabled_rows
        if row.module == "mod_IV" and row.collection_norm == "HOMD"
    ]
    if len(human_rows) > 1:
        die(
            "no more than one enabled human_uniprot row is allowed; found "
            + ", ".join(row.source_id for row in human_rows)
        )
    if len(homd_rows) > 1:
        die(
            "no more than one enabled HOMD row is allowed; found "
            + ", ".join(row.source_id for row in homd_rows)
        )

    for row in rows:
        if not row.enabled:
            continue

        # The 'rank' option drives HOMD representative selection only. Flag it
        # early if it strays onto a non-HOMD row so a misplaced option is caught
        # at plan time rather than silently ignored downstream.
        if "rank" in row.options and row.collection_norm != "HOMD":
            warn(f"{row.source_id}: 'rank' option applies only to the HOMD "
                 f"collection; it will be ignored here")

        if row.module == "mod_IV":
            if not row.collection:
                die(f"{row.source_id}: mod_IV (supported collection) rows require collection")
            if row.collection_norm not in SUPPORTED_COLLECTIONS:
                die(f"{row.source_id}: unsupported collection={row.collection!r}")
            # HOMD-only: validate the representative-selection rank if provided.
            # Empty/absent means "no dereplication" (keep every body-site match),
            # parallel to omitting 'sites='. A typo dies here instead of at fetch.
            if row.collection_norm == "HOMD":
                rank = row.options.get("rank", "").strip().lower()
                if rank and rank not in VALID_HOMD_RANKS:
                    die(f"{row.source_id}: invalid rank={row.options.get('rank')!r}; "
                        f"choose {'|'.join(sorted(VALID_HOMD_RANKS))} "
                        f"(or omit to keep every body-site match)")

        elif row.module == "mod_I":
            if not row.input_list:
                die(f"{row.source_id}: mod_I (proteome) rows require input_list")
            if row.fasta_file:
                die(f"{row.source_id}: mod_I (proteome) rows should not have fasta_file")

        elif row.module == "mod_II":
            if not row.input_list:
                die(f"{row.source_id}: mod_II (protein-accession) rows require input_list")
            if row.fasta_file:
                die(f"{row.source_id}: mod_II (protein-accession) rows should not have fasta_file")

        elif row.module == "mod_III":
            if not row.fasta_file:
                die(f"{row.source_id}: mod_III (user FASTA) rows require fasta_file")
            if row.input_list:
                die(f"{row.source_id}: mod_III (user FASTA) rows should not have input_list")

        # Fetch backend: ncbi is default; uniprot/uniparc/pdb select alternates.
        if row.module in ("mod_I", "mod_II"):
            fs = row.options.get("fetch_source", "ncbi").lower()
            if fs not in VALID_FETCH_SOURCES:
                die(f"{row.source_id}: invalid fetch_source={fs!r}; "
                    f"choose {'|'.join(sorted(VALID_FETCH_SOURCES))}")
            # UniParc and PDB are protein-accession backends only; there is no
            # proteome-level concept for them, so they require mod_II.
            if fs in {"uniparc", "pdb"} and row.module != "mod_II":
                die(f"{row.source_id}: fetch_source={fs} is only valid for "
                    f"mod_II (protein-accession) rows, not {row.module}")

        if strict_files:
            for field, p in [
                ("input_list", row.input_list),
                ("fasta_file", row.fasta_file),
                ("metadata_file", row.metadata_file),
            ]:
                if p and not (p.startswith("$") or p.startswith("~")) and not Path(p).is_file():
                    die(f"{row.source_id}: {field} not found: {p}")


def write_normalized_plan(rows: List[SourceRow], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    cols = [
        "source_id", "enabled", "source_label", "source_group", "module", "collection",
        "input_list", "fasta_file", "metadata_file", "options", "version", "citation", "notes", "status",
    ]
    with out_path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for row in rows:
            status = "enabled" if row.enabled else "disabled"
            w.writerow({
                "source_id": row.source_id,
                "enabled": "true" if row.enabled else "false",
                "source_label": row.source_label,
                "source_group": row.source_group,
                "module": row.module,
                "collection": row.collection_norm,
                "input_list": row.input_list,
                "fasta_file": row.fasta_file,
                "metadata_file": row.metadata_file,
                "options": row.options_raw,
                "version": row.version,
                "citation": row.citation,
                "notes": row.notes,
                "status": status,
            })


def main() -> None:
    ap = argparse.ArgumentParser(description="Prepare normalized maniFasta source plan from source registry")
    ap.add_argument("--registry", required=True, help="000.maniFasta.source_registry.tsv")
    ap.add_argument("--main-dir", required=True, help="mainDIR")
    ap.add_argument("--setup-dir", required=True, help="setupDIR / legacy relative path base")
    ap.add_argument("--ref-dir", required=True, help="run_data directory")
    ap.add_argument("--tmp-dir", required=True, help="TMPDIR")
    ap.add_argument("--out-normalized-plan", required=True, help="Output normalized source plan TSV")
    ap.add_argument("--no-strict-files", action="store_true", help="Do not fail when listed input files are missing")
    args = ap.parse_args()

    registry = Path(args.registry)
    main_dir = Path(args.main_dir)
    setup_dir = Path(args.setup_dir)
    ref_dir = Path(args.ref_dir)
    tmp_dir = Path(args.tmp_dir)
    normalized_plan = Path(args.out_normalized_plan)

    rows = read_registry(registry, main_dir, setup_dir)
    validate_rows(rows, strict_files=not args.no_strict_files)
    write_normalized_plan(rows, normalized_plan)

    enabled = [r for r in rows if r.enabled]
    info(f"Source registry rows: {len(rows)} total, {len(enabled)} enabled")
    info(f"Main directory: {main_dir}")
    info(f"Legacy path base: {setup_dir}")
    info(f"Normalized source plan: {normalized_plan}")

    # Detailed source plan summary for Slurm logs.
    print("[INFO] Enabled source plan:", file=sys.stderr)
    for r in enabled:
        detail = r.collection_norm if r.module == "mod_IV" else (
            r.input_list or r.fasta_file
        )
        pieces = [
            f"label={r.source_label}",
            f"group={r.source_group or 'NA'}",
            f"module={r.module}",
            f"detail={detail or 'NA'}",
        ]
        if r.source_dir:
            pieces.append(f"source_dir={r.source_dir}")
        if r.options_raw:
            pieces.append(f"options={r.options_raw}")
        if r.version:
            pieces.append(f"version={r.version}")
        if r.input_list:
            pieces.append(f"input_list={r.input_list}")
        if r.fasta_file:
            pieces.append(f"fasta_file={r.fasta_file}")
        if r.metadata_file:
            pieces.append(f"metadata_file={r.metadata_file}")
        if r.notes:
            pieces.append(f"notes={r.notes}")
        print(f"[INFO]   {r.source_id}: " + " | ".join(pieces), file=sys.stderr)

    # Keep these variables referenced so linters do not complain if checks are added later.
    _ = (ref_dir, tmp_dir)

    print("[INFO] Source planner complete.", file=sys.stderr)


if __name__ == "__main__":
    main()
