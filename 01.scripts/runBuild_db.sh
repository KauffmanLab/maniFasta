#!/bin/bash

## Originated by Christopher Handelmann

set -euo pipefail

ts() { date +%H:%M:%S; }

phase_banner() {
    local msg="$1"
    echo "" >&2
    echo "============================================================" >&2
    echo "[$(ts)] $msg" >&2
    echo "============================================================" >&2
}

module_banner() {
    local msg="$1"
    echo "" >&2
    echo "---- [$(ts)] $msg ----" >&2
}

source "${Config_file}"
source "${setupDIR}/000.maniFasta.local.env"


# Backfill directory variables for older configs that do not yet define them.
: "${setupDIR:=${mainDIR}/00.setup}"
: "${scriptsDIR:=${mainDIR}/01.scripts}"
: "${dbRootDIR:=${mainDIR}}"
: "${logDIR:=${mainDIR}}"
: "${run_reference_data_subdir:=run_data}"
: "${run_tmp_subdir:=tmp}"

# ── Build output directory ───────────────────────────────────────────────────
# Define the timestamp once at the very start of the run. All run-scoped paths
# derive from this single value so logs, downloaded reference data, outputs, and
# run documentation are internally consistent.
# Named: <timestamp>_<run_label>, e.g. 20260608_090000_AllOralsDB_v1.0
RUN_LABEL="${run_label:-test}"
RUN_TIMESTAMP="$(date +%Y%m%d_%H%M%S)"
RUN_STAMP="${RUN_TIMESTAMP}_${RUN_LABEL}"
RUN_DIR="${dbRootDIR}/${RUN_STAMP}"
RUN_REFDIR="${RUN_DIR}/${run_reference_data_subdir}"
RUN_TMPDIR="${RUN_DIR}/${run_tmp_subdir}"

export RUN_LABEL RUN_TIMESTAMP RUN_STAMP RUN_DIR RUN_REFDIR RUN_TMPDIR

# Re-source config after exporting RUN_REFDIR/RUN_TMPDIR so derived config paths
# such as refDIR, TMPDIR, user_fasta_file, proteomes_manifest, etc. resolve into
# this run directory.
source "${Config_file}"

: "${setupDIR:=${mainDIR}/00.setup}"
: "${scriptsDIR:=${mainDIR}/01.scripts}"
: "${refDIR:=${RUN_REFDIR}}"
: "${TMPDIR:=${RUN_TMPDIR}}"

mkdir -p "${RUN_DIR}" "${refDIR}" "${TMPDIR}"

# Record the exact maniFasta version used for this build.
if command -v git >/dev/null 2>&1 &&
   git -C "${mainDIR}" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
    mf_git_commit="$(git -C "${mainDIR}" rev-parse HEAD)"
    mf_git_tag="$(git -C "${mainDIR}" describe --tags --exact-match HEAD 2>/dev/null || true)"
    [[ -n "${mf_git_tag}" ]] || mf_git_tag="Untagged"
else
    mf_git_commit="Unavailable"
    mf_git_tag="Unavailable"
fi

printf 'field\tvalue\nmanifasta_git_tag\t%s\nmanifasta_git_commit\t%s\n' \
    "${mf_git_tag}" "${mf_git_commit}" \
    > "${RUN_DIR}/maniFasta_software_provenance.tsv"

echo "[INFO] [$(ts)] Run timestamp          : ${RUN_TIMESTAMP}" >&2
echo "[INFO] [$(ts)] Build output directory : ${RUN_DIR}" >&2
echo "[INFO] [$(ts)] Run reference data dir : ${refDIR}" >&2
echo "[INFO] [$(ts)] Run temp directory     : ${TMPDIR}" >&2

# ── allow-missing accession handling ─────────────────────────────────────────
# When allow_missing_accessions=true (config), unreturned accessions are skipped
# by the fetchers and recorded in this shared file, which is then handed to
# build_db.py so they appear in the build summary's failed/skipped section.
FETCH_MISSING="${TMPDIR}/fetch_missing_accessions.tsv"
ALLOW_MISSING_FLAG=()
if [[ "${allow_missing_accessions:-false}" == "true" ]]; then
    ALLOW_MISSING_FLAG=( --allow-missing --missing-out "$FETCH_MISSING" )
    echo "[INFO] [$(ts)] allow_missing_accessions=true: unreturned accessions will be skipped and reported." >&2
else
    echo "[INFO] [$(ts)] allow_missing_accessions=false: any unreturned accession aborts the build (strict)." >&2
fi

# ── UniParc representative-organism policy ───────────────────────────────────
# organism_pick (config) controls how mod_II_fetchUniParcAccessions.py chooses ONE
# representative organism when a UPI maps to multiple organisms across its
# UniParc cross-references. Applied globally to every UniParc fetch in this
# build (both the --plan whole-source path and the partitioned --input path).
#   best   : restrict to the highest-authority (active, db-rank) cross-refs,
#            then take the most frequent organism within that tier (default)
#   common : take the most frequent organism across ALL cross-references
#   first  : take the first cross-reference's organism as listed by UniProt
organism_pick="${organism_pick:-best}"
case "$organism_pick" in
    best|common|first) ;;
    *)
        echo "ERROR: [$(ts)] invalid organism_pick=${organism_pick}; choose best | common | first" >&2
        exit 1
        ;;
esac
echo "[INFO] [$(ts)] UniParc organism-pick policy: ${organism_pick}" >&2

# ── UniProt proteome-type policy ─────────────────────────────────────────────
# proteome_type (config) is the global default proteome type used by
# mod_I_fetchUniProtProteomes.py for UniProt-backed proteome (Module I) sources that
# do not set their own uniprot_proteome_type option in the source registry.
# A per-source uniprot_proteome_type=... option still overrides this.
#   best           : reference > representative > other > redundant by rank (default)
#   reference      : require a Reference proteome (falls back to best if none)
#   representative : require a Representative proteome (falls back to best if none)
#   any            : take the first proteome the search returns
#   pan            : resolve to the linked pan-proteome (falls back to the proteome)
proteome_type="${proteome_type:-best}"
case "$proteome_type" in
    best|reference|representative|any|pan) ;;
    *)
        echo "ERROR: [$(ts)] invalid proteome_type=${proteome_type}; choose best | reference | representative | any | pan" >&2
        exit 1
        ;;
esac
echo "[INFO] [$(ts)] UniProt proteome-type policy: ${proteome_type}" >&2

# ── UniProt taxon-scope inclusion policy ─────────────────────────────────────
# Global defaults for mod_II uniprot_scope=taxon sources: whether a taxonomy_id
# bulk pull includes TrEMBL (else reviewed Swiss-Prot only) and isoforms. Passed
# to mod_II_fetchUniProtAccessions.py; a per-source uniprot_include_trembl= /
# uniprot_include_isoform= option in the source registry still overrides these.
#   uniprot_include_trembl   true : include TrEMBL + Swiss-Prot | false : reviewed only (default)
#   uniprot_include_isoform  true : include isoforms            | false : canonical only (default)
uniprot_include_trembl="${uniprot_include_trembl:-false}"
uniprot_include_isoform="${uniprot_include_isoform:-false}"
case "$uniprot_include_trembl" in
    true|false) ;;
    *)
        echo "ERROR: [$(ts)] invalid uniprot_include_trembl=${uniprot_include_trembl}; choose true | false" >&2
        exit 1
        ;;
esac
case "$uniprot_include_isoform" in
    true|false) ;;
    *)
        echo "ERROR: [$(ts)] invalid uniprot_include_isoform=${uniprot_include_isoform}; choose true | false" >&2
        exit 1
        ;;
esac
echo "[INFO] [$(ts)] UniProt taxon inclusion: trembl=${uniprot_include_trembl} isoform=${uniprot_include_isoform}" >&2

# Override placeholder config values with real timestamped paths
db_fasta_output="${RUN_DIR}/${RUN_LABEL}.fasta"
db_map_output="${RUN_DIR}/${RUN_LABEL}.manifest.tsv"

# ── Source registry planner ──────────────────────────────────────────────────
phase_banner "PHASE 0: Initialize run and parse source registry"

if [[ -n "${source_registry:-}" ]]; then
    echo "[INFO] [$(ts)] Source registry detected: ${source_registry}" >&2
    SOURCE_PLAN_NORMALIZED="${TMPDIR}/source_plan.normalized.tsv"

    STRICT_FILES_FLAG=()
    if [[ "${fail_on_missing_source_file:-true}" != "true" ]]; then
        STRICT_FILES_FLAG=( --no-strict-files )
    fi

    python3 "${scriptsDIR}/source_planner.py" \
        --registry "${source_registry}" \
        --main-dir "${mainDIR}" \
        --setup-dir "${registry_path_base:-${setupDIR}}" \
        --ref-dir "${refDIR}" \
        --tmp-dir "${TMPDIR}" \
        --out-normalized-plan "${SOURCE_PLAN_NORMALIZED}" \
        "${STRICT_FILES_FLAG[@]}"

    # Directly derive mod_IV (supported collection) runtime variables from the
    # normalized source plan. Helper scripts still expect these variables, but
    # they now come from source_plan.normalized.tsv.
    SOURCE_PLAN_SUPPORTED_ENV="${TMPDIR}/source_plan.supported_collections.env"

    python3 - "${SOURCE_PLAN_NORMALIZED}" "${refDIR}" "${SOURCE_PLAN_SUPPORTED_ENV}" "${homd_prokka_base_url:-}" <<'PY'
import csv
import datetime
import re
import shlex
import sys
import urllib.request
from pathlib import Path

def _ts(): return datetime.datetime.now().strftime("%H:%M:%S")

source_plan = Path(sys.argv[1])
ref_dir = Path(sys.argv[2])
out_env = Path(sys.argv[3])
# Unversioned HOMD PROKKA root from the config (.../genomes/PROKKA). Used only to
# auto-detect the newest release when the registry omits genomic_refseq_version.
homd_base_url = (sys.argv[4] if len(sys.argv) > 4 else "").strip()

TRUE = {"true", "t", "yes", "y", "1", "on"}
FALSE = {"false", "f", "no", "n", "0", "off", ""}

def parse_options(raw: str) -> dict[str, str]:
    opts = {}
    for item in (raw or "").split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" not in item:
            raise SystemExit(f"ERROR: bad option {item!r}; expected key=value")
        k, v = item.split("=", 1)
        opts[k.strip()] = v.strip()
    return opts

def parse_bool(value: str, default: bool = True) -> bool:
    if value is None:
        return default
    v = str(value).strip().lower()
    if v in TRUE:
        return True
    if v in FALSE:
        return False
    raise SystemExit(f"ERROR: invalid boolean option value: {value!r}")

def q(value) -> str:
    if isinstance(value, bool):
        value = "true" if value else "false"
    return shlex.quote(str(value))

def normalize_homd_version(v: str) -> str:
    """Canonicalize a HOMD PROKKA release token: 'v11.03' -> 'V11.03'.

    HOMD's directory names use an uppercase 'V'; a lowercase 'v' 404s. Tokens that
    are not V<major>.<minor> (e.g. 'current') and empty strings pass through
    untouched.
    """
    v = (v or "").strip()
    m = re.fullmatch(r"[Vv]?(\d+\.\d+)", v)
    return "V" + m.group(1) if m else v

def resolve_latest_homd_version(base_url: str) -> str:
    """Return the newest V<major>.<minor> release directory under base_url.

    Parses the FTP autoindex listing and picks the highest version numerically
    (so V11.10 sorts above V11.9). The literal directory name is preserved.
    """
    index_url = base_url.rstrip("/") + "/"
    try:
        req = urllib.request.Request(index_url, headers={"User-Agent": "maniFasta/1.0"})
        with urllib.request.urlopen(req, timeout=60) as r:
            html = r.read().decode("utf-8", "replace")
    except Exception as exc:
        raise SystemExit(
            f"ERROR: could not list HOMD PROKKA releases at {index_url} ({exc}).\n"
            f"       Set genomic_refseq_version= on the HOMD row of the source registry "
            f"to pin a release explicitly."
        )
    found = re.findall(r'href="[^"]*?(V\d+\.\d+)/?"', html, flags=re.IGNORECASE)
    if not found:
        raise SystemExit(
            f"ERROR: no V<major>.<minor> release directories found at {index_url}.\n"
            f"       Set genomic_refseq_version= on the HOMD row of the source registry "
            f"to pin a release explicitly."
        )
    def _key(v):
        mj, mn = re.match(r"[Vv](\d+)\.(\d+)", v).groups()
        return (int(mj), int(mn))
    return max({normalize_homd_version(f) for f in found}, key=_key)

# Defaults when supported collections are absent/disabled.
env = {
    "use_human": False,
    "use_human_download": True,
    "human_protein_set": "canonical",
    "human_fasta_file": ref_dir / "human" / "human_proteins.canonical.fasta",

    "use_homd": False,
    "use_homd_download": True,
    "homd_genomic_refseq_version": "",
    "homd_sites": "",
    "homd_rank": "",
    "homd_tax_raw_file": ref_dir / "homd" / "homd_taxonomy_raw.tsv",
    "homd_tax_file": ref_dir / "homd" / "homd_taxonomy.tsv",
    "homd_filtered_tax": ref_dir / "homd" / "homd_taxonomy_filtered.tsv",
    "gca_info_file": ref_dir / "homd" / "GCA_ID_info_.txt",
    "homd_fasta_file": ref_dir / "homd" / "homd_proteins.fasta",
    "homd_proteomes_dir": ref_dir / "homd" / "proteomes",

    "crap_set": "none",
    # cRAP is run-scoped: the exact combined FASTA used for this build is
    # written flat into run_data/ alongside other staged FASTAs.
    "crap_dir": ref_dir,
    "crap_fasta": "",
}

with source_plan.open(newline="", encoding="utf-8") as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    for row in reader:
        if (row.get("enabled") or "").strip().lower() != "true":
            continue
        if (row.get("status") or "").strip().lower() == "disabled":
            continue
        if (row.get("module") or "").strip() != "mod_IV":
            continue

        source_id = (row.get("source_id") or "").strip()
        coll = (row.get("collection") or "").strip()
        if coll.lower() == "homd":
            coll = "HOMD"
        elif coll.lower() == "crap":
            coll = "cRAP"
        opts = parse_options(row.get("options") or "")

        if coll == "human_uniprot":
            protein_set = opts.get("protein_set", "canonical")
            env["use_human"] = True
            env["use_human_download"] = parse_bool(opts.get("download", "true"), True)
            env["human_protein_set"] = protein_set
            env["human_fasta_file"] = ref_dir / "human" / f"human_proteins.{protein_set}.fasta"
            print(f"[INFO] [{_ts()}] Supported collection parsed: {source_id} -> HUMAN protein_set={protein_set} download={env['use_human_download']}", file=sys.stderr)

        elif coll == "HOMD":
            # Version resolution — the single source of truth for the HOMD release.
            # Precedence: genomic_refseq_version= option, then the row's version
            # column, then auto-detect the newest release under homd_prokka_base_url.
            # The resolved value names GCA_ID_info_<version>.txt, is recorded in the
            # methods manifest, and builds every PROKKA URL in the mod_IV fetchers.
            download = parse_bool(opts.get("download", "true"), True)
            version = normalize_homd_version(
                opts.get("genomic_refseq_version", (row.get("version") or "").strip())
            )
            if not version:
                if not homd_base_url:
                    raise SystemExit(
                        "ERROR: HOMD genomic_refseq_version is unset and homd_prokka_base_url "
                        "is empty in the config, so the latest release cannot be auto-detected."
                    )
                if not download:
                    raise SystemExit(
                        "ERROR: HOMD genomic_refseq_version is unset and download=false, so the "
                        "latest release cannot be auto-detected. Pin genomic_refseq_version= on "
                        "the HOMD row to match the pre-staged files."
                    )
                version = resolve_latest_homd_version(homd_base_url)
                print(f"[INFO] [{_ts()}] HOMD genomic_refseq_version unset -> auto-detected "
                      f"latest release: {version}", file=sys.stderr)
            sites = opts.get("sites", "")
            sites = " ".join(s.strip() for s in sites.replace(",", " ").split() if s.strip())
            # Representative-selection rank (dereplication level). Empty = keep
            # every body-site match. Lowercased so it matches the fetch script's
            # species|genus|family guard; the planner already validated it.
            rank = opts.get("rank", "").strip().lower()
            env["use_homd"] = True
            env["use_homd_download"] = download
            env["homd_genomic_refseq_version"] = version
            env["homd_sites"] = sites
            env["homd_rank"] = rank
            env["gca_info_file"] = ref_dir / "homd" / f"GCA_ID_info_{version}.txt"
            print(f"[INFO] [{_ts()}] Supported collection parsed: {source_id} -> HOMD version={version} sites={sites or 'all'} rank={rank or 'none'} download={env['use_homd_download']}", file=sys.stderr)

        elif coll == "cRAP":
            env["crap_set"] = opts.get("set", "ccp")
            print(f"[INFO] [{_ts()}] Supported collection parsed: {source_id} -> cRAP set={env['crap_set']}", file=sys.stderr)

        else:
            raise SystemExit(f"ERROR: unsupported collection={coll!r} in {source_id}")

out_env.parent.mkdir(parents=True, exist_ok=True)
with out_env.open("w", encoding="utf-8") as out:
    out.write("# Auto-generated by runBuild_db.sh from source_plan.normalized.tsv\n")
    out.write("# Direct mod_IV (supported collection) runtime variables.\n")
    for key in sorted(env):
        out.write(f"export {key}={q(env[key])}\n")

print(f"[INFO] [{_ts()}] Supported-collection env written: {out_env}", file=sys.stderr)
PY

    # shellcheck disable=SC1090
    source "${SOURCE_PLAN_SUPPORTED_ENV}"
    echo "[INFO] [$(ts)] Supported-collection variables loaded from normalized source plan: ${SOURCE_PLAN_SUPPORTED_ENV}" >&2
else
    echo "ERROR: [$(ts)] source_registry is not set. maniFasta now requires a source registry." >&2
    exit 1
fi

# ── Run artifact packaging and cleanup ───────────────────────────────────────
# On every exit, collect the files needed to understand/reproduce the run.
# On success, remove the run-scoped TMPDIR after packaging.
# On failure, preserve TMPDIR for debugging.
_package_run_artifacts_and_cleanup() {
    local exit_code=$?

    if [[ -n "${RUN_DIR:-}" && -d "$RUN_DIR" ]]; then
        phase_banner "PHASE 3: Package run provenance, scripts, setup files, and logs"

        local run_setup="${RUN_DIR}/run_setup"
        local run_scripts="${RUN_DIR}/run_scripts"
        local run_logs="${RUN_DIR}/run_logs"

        mkdir -p "$run_setup" "$run_scripts" "$run_logs"

        # ── run_setup: user/run-definition files and source inputs ───────────
        [[ -n "${Config_file:-}" && -f "${Config_file}" ]] && \
            cp -f "$Config_file" "$run_setup/" && \
            echo "[INFO] [$(ts)] Copied config to run_setup/" >&2

        [[ -n "${source_registry:-}" && -f "${source_registry}" ]] && \
            cp -f "$source_registry" "$run_setup/" && \
            echo "[INFO] [$(ts)] Copied source registry to run_setup/" >&2

        # All input files referenced in source_plan.normalized.tsv. These are
        # the setup/source files supplied by the user or curated for this build.
        if [[ -n "${SOURCE_PLAN_NORMALIZED:-}" && -f "${SOURCE_PLAN_NORMALIZED}" ]]; then
            python3 - "$SOURCE_PLAN_NORMALIZED" "$run_setup" << 'PYSETUP'
import csv
import datetime
import sys
from pathlib import Path
import shutil

def _ts(): return datetime.datetime.now().strftime("%H:%M:%S")

plan = Path(sys.argv[1])
out_dir = Path(sys.argv[2])
copied = set()

with plan.open(newline="", encoding="utf-8") as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    for row in reader:
        if (row.get("enabled") or "").strip().lower() != "true":
            continue
        if (row.get("status") or "").strip().lower() == "disabled":
            continue
        for field in ["input_list", "fasta_file", "metadata_file"]:
            p = (row.get(field) or "").strip()
            if not p or p in copied or p.startswith("$") or p.startswith("~"):
                continue
            src = Path(p)
            if src.is_file():
                shutil.copy2(src, out_dir / src.name)
                copied.add(p)

print(f"[INFO] [{_ts()}] Copied {len(copied)} source input files to run_setup/", file=sys.stderr)
PYSETUP
        fi

        # ── run_scripts: copy the script set used to launch/build this run ───
        # Exclude caches, old/backups, and bytecode so the directory is compact.
        if [[ -n "${scriptsDIR:-}" && -d "${scriptsDIR}" ]]; then
            python3 - "$scriptsDIR" "$run_scripts" << 'PYSCRIPTS'
import datetime
import sys
from pathlib import Path
import shutil

def _ts(): return datetime.datetime.now().strftime("%H:%M:%S")

src_dir = Path(sys.argv[1])
dst_dir = Path(sys.argv[2])

ignore_dirs = {".git", "__pycache__", "oldies"}
ignore_suffixes = {".pyc", ".pyo"}
ignore_contains = [".bak", "~"]

copied = 0

for src in src_dir.rglob("*"):
    rel = src.relative_to(src_dir)
    if any(part in ignore_dirs for part in rel.parts):
        continue
    if src.is_dir():
        continue
    name = src.name
    if any(token in name for token in ignore_contains):
        continue
    if src.suffix in ignore_suffixes:
        continue

    dst = dst_dir / rel
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)
    copied += 1

print(f"[INFO] [{_ts()}] Copied {copied} script files to run_scripts/", file=sys.stderr)
PYSCRIPTS
        fi

        # ── run_logs: machine-readable run plans/checks and execution logs ───
        [[ -n "${SOURCE_PLAN_NORMALIZED:-}" && -f "${SOURCE_PLAN_NORMALIZED}" ]] && \
            cp -f "$SOURCE_PLAN_NORMALIZED" "$run_logs/"
        [[ -n "${SOURCE_PLAN_SUPPORTED_ENV:-}" && -f "${SOURCE_PLAN_SUPPORTED_ENV}" ]] && \
            cp -f "$SOURCE_PLAN_SUPPORTED_ENV" "$run_logs/"
        [[ -n "${SOURCE_STAGE_FLAGS:-}" && -f "${SOURCE_STAGE_FLAGS}" ]] && \
            cp -f "$SOURCE_STAGE_FLAGS" "$run_logs/"
        [[ -n "${PROT_FETCH_REPORT:-}" && -f "${PROT_FETCH_REPORT}" ]] && \
            cp -f "$PROT_FETCH_REPORT" "$run_logs/"
        [[ -n "${FETCH_MISSING:-}" && -f "${FETCH_MISSING}" ]] && \
            cp -f "$FETCH_MISSING" "$run_logs/" && \
            echo "[INFO] [$(ts)] Skipped/missing accessions report copied to run_logs/" >&2
        [[ -n "${Species_resolved:-}" && -f "${Species_resolved}" ]] && \
            cp -f "$Species_resolved" "$run_logs/"
        [[ -n "${proteomes_manifest:-}" && -f "${proteomes_manifest}" ]] && \
            cp -f "$proteomes_manifest" "$run_logs/"

        # Slurm log. Move it into run_logs when possible; on non-Slurm/local
        # runs this simply does nothing.
        if [[ -n "${SLURM_JOB_ID:-}" ]]; then
            local log="${logDIR:-${mainDIR}}/slurm-${SLURM_JOB_ID}.out"
            if [[ -f "$log" ]]; then
                mv -f "$log" "${run_logs}/slurm-${SLURM_JOB_ID}.out" || true
                echo "[INFO] [$(ts)] Slurm log moved to: run_logs/" >&2
            fi
        fi

        echo "[INFO] [$(ts)] Run setup files packaged in: ${run_setup}" >&2
        echo "[INFO] [$(ts)] Run scripts packaged in: ${run_scripts}" >&2
        echo "[INFO] [$(ts)] Run logs/checks packaged in: ${run_logs}" >&2
    fi

    if [[ $exit_code -eq 0 ]]; then
        echo "[INFO] [$(ts)] Run completed successfully. Cleaning up TMPDIR: ${TMPDIR}" >&2
        rm -rf "${TMPDIR:?}"
    else
        echo "[INFO] [$(ts)] Run failed. TMPDIR preserved for debugging: ${TMPDIR}" >&2
    fi
}
trap '_package_run_artifacts_and_cleanup' EXIT

phase_banner "PHASE 1: Collect and stage source data"
module_banner "Modules I–III: stage NCBI proteome lists, NCBI protein accessions, and user FASTAs"

# ═════════════════════════════════════════════════════════════════════════════
# SOURCE PLAN MODULE STAGING
#
# Enabled normalized-plan rows are routed by module:
#   mod_I   -> combined species/genome fetch list
#   mod_II  -> combined protein accession list + metadata
#   mod_III -> staged user FASTA + metadata
#   mod_IV  -> handled by the supported-collection sections
# ═════════════════════════════════════════════════════════════════════════════

# Staging directories for source-plan outputs
FETCH_STAGE="${TMPDIR}/fetch_lists_combined"
FASTA_STAGE="${TMPDIR}/fasta_override_combined"
PROT_STAGE="${TMPDIR}/prot_fetch_combined"
mkdir -p "$FETCH_STAGE" "$FASTA_STAGE" "$PROT_STAGE"

# These accumulate across source-plan rows
COMBINED_FETCH_LIST="${FETCH_STAGE}/combined_fetch_list.tsv"
COMBINED_FETCH_LIST_NCBI="${FETCH_STAGE}/combined_fetch_list.ncbi.tsv"
COMBINED_FETCH_LIST_UNIPROT="${FETCH_STAGE}/combined_fetch_list.uniprot.tsv"
COMBINED_FASTA="${FASTA_STAGE}/user_supplied_combined.fasta"
COMBINED_META="${FASTA_STAGE}/user_supplied_metadata.tsv"
COMBINED_PROT_LIST="${PROT_STAGE}/combined_protein_accessions.tsv"
COMBINED_PROT_LIST_NCBI="${PROT_STAGE}/combined_protein_accessions.ncbi.tsv"
COMBINED_UNIPROT_LIST="${PROT_STAGE}/combined_protein_accessions.uniprot.tsv"
COMBINED_UNIPARC_LIST="${PROT_STAGE}/combined_protein_accessions.uniparc.tsv"
COMBINED_PDB_LIST="${PROT_STAGE}/combined_protein_accessions.pdb.tsv"
COMBINED_PENDING_LIST="${PROT_STAGE}/combined_protein_accessions.pending.tsv"
COMBINED_PROT_META="${PROT_STAGE}/combined_protein_metadata.tsv"
PROT_FETCH_REPORT="${PROT_STAGE}/combined_proteins.fetch_report.tsv"
SOURCE_STAGE_FLAGS="${TMPDIR}/source_plan_stage_flags.env"

HAS_FETCH_ROWS=false
HAS_FETCH_UNIPROT_ROWS=false
HAS_FASTA_ROWS=false
HAS_PROT_ROWS=false

if [[ -n "${SOURCE_PLAN_NORMALIZED:-}" && -f "${SOURCE_PLAN_NORMALIZED}" ]]; then
    echo "[INFO] [$(ts)] Reading normalized source plan: ${SOURCE_PLAN_NORMALIZED}" >&2

    python3 - \
        "${SOURCE_PLAN_NORMALIZED}" \
        "${COMBINED_FETCH_LIST}" \
        "${COMBINED_FASTA}" \
        "${COMBINED_META}" \
        "${COMBINED_PROT_LIST}" \
        "${COMBINED_PROT_META}" \
        "${SOURCE_STAGE_FLAGS}" <<'PY'
import csv
import shlex
import sys
from pathlib import Path

(
    source_plan,
    combined_fetch_list,
    combined_fasta,
    combined_meta,
    combined_prot_list,
    combined_prot_meta,
    source_stage_flags,
) = [Path(x) for x in sys.argv[1:8]]

for p in [
    combined_fetch_list,
    combined_fasta,
    combined_meta,
    combined_prot_list,
    combined_prot_meta,
    source_stage_flags,
]:
    p.parent.mkdir(parents=True, exist_ok=True)
    if p.exists():
        p.unlink()

fetch_header_written = False
fasta_meta_header_written = False
prot_list_header_written = False
prot_meta_header_written = False

has_fetch_rows = False
has_fasta_rows = False
has_prot_rows = False

counts = {
    "mod_I_rows": 0,
    "mod_II_rows": 0,
    "mod_III_rows": 0,
    "mod_IV_rows": 0,
    "fetch_records": 0,
    "protein_accessions": 0,
    "user_fasta_records": 0,
    "metadata_rows": 0,
}


def info(msg: str) -> None:
    import datetime
    ts = datetime.datetime.now().strftime("%H:%M:%S")
    print(f"[INFO] [{ts}] {msg}", file=sys.stderr)


def die(msg: str) -> None:
    import datetime
    ts = datetime.datetime.now().strftime("%H:%M:%S")
    print(f"ERROR: [{ts}] {msg}", file=sys.stderr)
    raise SystemExit(1)


def clean_lines(path: Path) -> list[str]:
    return [
        line.rstrip("\n").replace("\r", "")
        for line in path.open(encoding="utf-8", errors="replace")
        if line.strip() and not line.startswith("#")
    ]


def read_table(path: Path) -> tuple[list[str], list[list[str]]]:
    lines = clean_lines(path)
    if not lines:
        die(f"no data rows found in {path}")
    header = lines[0].split("\t")
    rows = [line.split("\t") for line in lines[1:]]
    return header, rows


_FETCH_SPP_NAMES = ["species", "species name"]
_FETCH_ACC_NAMES = ["accession"]
_FETCH_TAX_NAMES = ["taxonid", "taxid"]
_FETCH_NOTE_NAMES = ["note"]
_FETCH_SOURCE_NAMES = ["fetch_source"]
_FETCH_CANON_HEADER = ["species", "accession", "taxonID", "note", "fetch_source"]


def _fetch_col_index(header_lower: list, names: list) -> int:
    return next((header_lower.index(n) for n in names if n in header_lower), -1)


def append_source_fetch_list(label: str, path_value: str) -> None:
    """Append one source's species-info rows to the combined fetch list.

    Columns are mapped by NAME (case-insensitive), not by position. Different
    mod_I (proteome) sources can have genuinely different native column orders
    (legacy 'species/accession/taxonID/note' vs newer
    'species/taxid/accession/fetch_source/category/subcategory/note'); reusing
    one source's header positionally for another's rows silently swaps
    accession <-> taxid whenever the orders differ, producing wrong-organism
    results (coincidental nuccore-UID matches) without ever raising an error.

    A per-row `fetch_source` column (ncbi|uniprot) is read the same way
    mod_II reads it for protein accessions, so a single mod_I source can mix
    NCBI and UniProt rows; mod_I_partition_proteome_sources.py splits them
    downstream. Absent column/value defaults to ncbi, so existing lists with
    no fetch_source column behave exactly as before.
    """
    global fetch_header_written, has_fetch_rows

    src = Path(path_value)
    if not src.is_file():
        die(f"fetch_list not found for '{label}': {src}")

    header, rows = read_table(src)
    header_lower = [h.strip().lower() for h in header]
    i_spp  = _fetch_col_index(header_lower, _FETCH_SPP_NAMES)
    i_acc  = _fetch_col_index(header_lower, _FETCH_ACC_NAMES)
    i_tax  = _fetch_col_index(header_lower, _FETCH_TAX_NAMES)
    i_note = _fetch_col_index(header_lower, _FETCH_NOTE_NAMES)
    i_fs   = _fetch_col_index(header_lower, _FETCH_SOURCE_NAMES)
    if i_spp == -1 and i_acc == -1 and i_tax == -1:
        die(
            f"fetch_list for '{label}' ({src}) has no recognizable "
            f"species/accession/taxonID column. Found columns: {header}"
        )

    def get(row: list, idx: int) -> str:
        return row[idx].strip() if 0 <= idx < len(row) else ""

    mode = "a" if fetch_header_written else "w"
    with combined_fetch_list.open(mode, encoding="utf-8", newline="") as out:
        if not fetch_header_written:
            out.write("source_label\t" + "\t".join(_FETCH_CANON_HEADER) + "\n")
        for row in rows:
            spp, acc, tax, note = get(row, i_spp), get(row, i_acc), get(row, i_tax), get(row, i_note)
            if not (spp or acc or tax):
                continue
            fetch_source = (get(row, i_fs) or "ncbi").lower()
            out.write(label + "\t" + "\t".join([spp, acc, tax, note, fetch_source]) + "\n")
            counts["fetch_records"] += 1

    fetch_header_written = True
    has_fetch_rows = True


def append_protein_accessions(label: str, path_value: str) -> None:
    """Append accession-fetch rows with explicit protein_id/fetch_id provenance."""
    global prot_list_header_written, has_prot_rows

    src = Path(path_value)
    if not src.is_file():
        die(f"protein accession input_list not found for '{label}': {src}")

    header, rows = read_table(src)
    lower = [h.strip().lower() for h in header]

    def idx(name: str) -> int:
        return lower.index(name) if name in lower else -1

    acc_idx = idx("accession")
    protein_id_idx = idx("protein_id")
    fetch_id_idx = idx("fetch_id")
    fetch_type_idx = idx("fetch_type")
    source_accession_idx = idx("source_accession")
    fetch_source_idx = idx("fetch_source")
    note_idx = idx("note")

    if acc_idx < 0 and fetch_id_idx < 0:
        die(f"no 'accession' or 'fetch_id' column in {src}. Found: {header}")

    out_header = [
        "source_label", "protein_id", "fetch_id", "fetch_type",
        "source_accession", "accession", "fetch_source", "note",
    ]

    mode = "a" if prot_list_header_written else "w"
    with combined_prot_list.open(mode, encoding="utf-8", newline="") as out:
        if not prot_list_header_written:
            out.write("\t".join(out_header) + "\n")
        for row in rows:
            while len(row) < len(header):
                row.append("")

            accession = row[acc_idx].strip() if acc_idx >= 0 else ""
            fetch_id = row[fetch_id_idx].strip() if fetch_id_idx >= 0 else accession
            protein_id = row[protein_id_idx].strip() if protein_id_idx >= 0 else (fetch_id or accession)
            fetch_type = row[fetch_type_idx].strip() if fetch_type_idx >= 0 else "accession"
            source_accession = row[source_accession_idx].strip() if source_accession_idx >= 0 else (accession or fetch_id)
            note = row[note_idx].strip() if note_idx >= 0 else ""
            fetch_source = row[fetch_source_idx].strip() if fetch_source_idx >= 0 else ""

            if not fetch_id:
                die(f"empty fetch_id/accession in protein input for {label}: {src}")
            if not protein_id:
                die(f"empty protein_id in protein input for {label}: {src}")

            out.write("\t".join([
                label,
                protein_id,
                fetch_id,
                fetch_type or "accession",
                source_accession or "NA",
                accession or fetch_id,
                fetch_source,
                note,
            ]) + "\n")
            counts["protein_accessions"] += 1

    prot_list_header_written = True
    has_prot_rows = True

METADATA_CANONICAL_COLUMNS = [
    "protein_id", "source", "ncbi_taxid", "genus", "species",
    "strain_or_gene", "gene", "protein_name", "assembly", "description",
    "fetch_id", "fetch_type", "source_accession", "returned_fasta_id",
]


def append_metadata_file(label: str, path_value: str, out_path: Path, header_written_attr: str) -> bool:
    """Append metadata by column name into the canonical build_db.py schema.

    Source metadata tables are allowed to contain different optional columns and
    column orders. Every input is projected into one canonical superset before it
    is appended, so a short table cannot define the positional meaning of later,
    wider tables. The registry source_label deliberately overrides any source
    value carried by the input table.
    """
    global fasta_meta_header_written, prot_meta_header_written

    src = Path(path_value)
    if not src.is_file():
        die(f"metadata file not found: {src}")

    lines = clean_lines(src)
    if not lines:
        die(f"no data in metadata file: {src}")

    header = [h.strip() for h in lines[0].split("\t")]
    lower = [h.lower() for h in header]
    if "protein_id" not in lower:
        die(f"metadata file is missing required column 'protein_id': {src}; found {header}")
    if len(lower) != len(set(lower)):
        duplicates = sorted({name for name in lower if lower.count(name) > 1})
        die(f"metadata file has duplicate column names: {src}; duplicates={duplicates}")

    input_index = {name: idx for idx, name in enumerate(lower)}
    unknown = [name for name in lower if name not in METADATA_CANONICAL_COLUMNS]
    if unknown:
        info(f"Ignoring unsupported metadata columns in {src.name}: {unknown}")

    already_written = (
        fasta_meta_header_written if header_written_attr == "fasta"
        else prot_meta_header_written
    )
    mode = "a" if already_written else "w"

    with out_path.open(mode, encoding="utf-8", newline="") as out:
        if not already_written:
            out.write("\t".join(METADATA_CANONICAL_COLUMNS) + "\n")
        for line in lines[1:]:
            raw = line.split("\t")
            values = {
                column: (raw[idx].strip() if idx < len(raw) else "")
                for column, idx in input_index.items()
            }
            values["source"] = label
            out.write("\t".join(values.get(column, "") for column in METADATA_CANONICAL_COLUMNS) + "\n")
            counts["metadata_rows"] += 1

    if header_written_attr == "fasta":
        fasta_meta_header_written = True
    else:
        prot_meta_header_written = True

    return True

def append_stub_metadata_from_fasta(label: str, fasta_path: Path) -> int:
    global fasta_meta_header_written

    header = METADATA_CANONICAL_COLUMNS

    rows = []
    with fasta_path.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                pid = line[1:].split()[0].strip()
                values = {column: "" for column in METADATA_CANONICAL_COLUMNS}
                values.update({
                    "protein_id": pid,
                    "source": label,
                    "ncbi_taxid": "32644",
                })
                rows.append([values[column] for column in METADATA_CANONICAL_COLUMNS])

    mode = "a" if fasta_meta_header_written else "w"
    with combined_meta.open(mode, encoding="utf-8", newline="") as out:
        if not fasta_meta_header_written:
            out.write("\t".join(header) + "\n")
        for row in rows:
            out.write("\t".join(row) + "\n")

    fasta_meta_header_written = True
    counts["metadata_rows"] += len(rows)
    return len(rows)


def append_user_fasta(label: str, fasta_value: str, metadata_value: str) -> None:
    global has_fasta_rows

    src = Path(fasta_value)
    if not src.is_file():
        die(f"fasta_file not found for '{label}': {src}")

    first = src.open("rb").read(1)
    if first != b">":
        die(
            f"fasta_file for '{label}' does not appear to be FASTA. "
            f"Expected first character '>', found {first!r}. Path: {src}"
        )

    info(f"Appending user_fasta for {label}: {src}")
    with src.open("r", encoding="utf-8", errors="replace") as inp, \
         combined_fasta.open("a", encoding="utf-8", newline="") as out:
        copied_any = False
        for line in inp:
            out.write(line)
            copied_any = True
            if line.startswith(">"):
                counts["user_fasta_records"] += 1
        if copied_any and not str(line).endswith("\n"):
            out.write("\n")

    if metadata_value:
        info(f"Using supplied metadata for {label}: {metadata_value}")
        append_metadata_file(label, metadata_value, combined_meta, "fasta")
    else:
        n = append_stub_metadata_from_fasta(label, src)
        info(f"  -> {n} stub metadata rows written for {label}")

    has_fasta_rows = True


with source_plan.open(newline="", encoding="utf-8") as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    required = {
        "source_id", "enabled", "source_label", "module",
        "collection", "input_list", "fasta_file", "metadata_file", "status",
    }
    missing = required - set(reader.fieldnames or [])
    if missing:
        die(f"normalized source plan missing columns: {sorted(missing)}")

    for row in reader:
        enabled = (row.get("enabled") or "").strip().lower() == "true"
        status = (row.get("status") or "").strip().lower()
        if not enabled or status == "disabled":
            continue

        source_id = (row.get("source_id") or "").strip()
        label = (row.get("source_label") or "").strip()
        module = (row.get("module") or "").strip()
        collection = (row.get("collection") or "").strip()
        input_list = (row.get("input_list") or "").strip()
        fasta_file = (row.get("fasta_file") or "").strip()
        metadata_file = (row.get("metadata_file") or "").strip()


        # UniProt/UniParc/PDB-backed whole sources are fetched after staging by
        # their dedicated fetchers via the --plan path. Defer them here so they
        # are NOT also routed through the NCBI partition path, which would
        # double-fetch and duplicate every protein in the build. (Per-ROW
        # uniprot/uniparc/pdb declared via a fetch_source column inside an
        # otherwise-NCBI list still flow through staging -> partition ->
        # the *_fromlist fetchers; only whole-source declarations are deferred.)
        _opts = {}
        for _it in (row.get("options") or "").split(";"):
            if "=" in _it:
                _k, _v = _it.split("=", 1)
                _opts[_k.strip().lower()] = _v.strip()
        _fs = _opts.get("fetch_source", "ncbi").lower()
        if _fs in ("uniprot", "uniparc", "pdb"):
            info(f"Deferring {_fs}-backed source: {source_id}")
            continue


        if module == "mod_IV":
            counts["mod_IV_rows"] += 1
            info(f"Supported collection will be handled later: {source_id} ({collection})")
        elif module == "mod_I":
            counts["mod_I_rows"] += 1
            info(f"Queuing mod_I (proteome) for {label}: {input_list}")
            append_source_fetch_list(label, input_list)
        elif module == "mod_II":
            counts["mod_II_rows"] += 1
            info(f"Queuing mod_II (protein accessions) for {label}: {input_list}")
            append_protein_accessions(label, input_list)
            if metadata_file:
                info(f"Queuing mod_II metadata for {label}: {metadata_file}")
                append_metadata_file(label, metadata_file, combined_prot_meta, "protein")
        elif module == "mod_III":
            counts["mod_III_rows"] += 1
            append_user_fasta(label, fasta_file, metadata_file)
        else:
            die(f"unsupported module in normalized source plan: {module!r} for {source_id}")

with source_stage_flags.open("w", encoding="utf-8") as out:
    def export_bool(name: str, value: bool) -> None:
        out.write(f"export {name}={'true' if value else 'false'}\n")
    export_bool("HAS_FETCH_ROWS", has_fetch_rows)
    export_bool("HAS_FASTA_ROWS", has_fasta_rows)
    export_bool("HAS_PROT_ROWS", has_prot_rows)

info(
    "Source-plan staging complete: "
    + ", ".join(f"{k}={v}" for k, v in counts.items())
)
info(f"Stage flags: {source_stage_flags}")
PY

    # shellcheck disable=SC1090
    source "${SOURCE_STAGE_FLAGS}"
    echo "[INFO] [$(ts)] Source-plan stage flags loaded from: ${SOURCE_STAGE_FLAGS}" >&2

else
    echo "[INFO] [$(ts)] No normalized source plan found; skipping registry-driven source staging." >&2
fi

# Append any later fetcher-generated metadata through the same canonical,
# name-aligned schema used during initial source staging. This prevents a table
# with one optional-column layout from changing the positional interpretation of
# a later table with another layout.
append_metadata_normalized() {
    local input_file="$1"
    local output_file="$2"

    python3 - "$input_file" "$output_file" <<'PYMETA'
import csv
import sys
from pathlib import Path

source = Path(sys.argv[1])
target = Path(sys.argv[2])
canonical = [
    "protein_id", "source", "ncbi_taxid", "genus", "species",
    "strain_or_gene", "gene", "protein_name", "assembly", "description",
    "fetch_id", "fetch_type", "source_accession", "returned_fasta_id",
]

if not source.is_file():
    raise SystemExit(f"ERROR: metadata file not found: {source}")

with source.open("r", encoding="utf-8", errors="replace", newline="") as fh:
    reader = csv.DictReader(
        (line for line in fh if line.strip() and not line.startswith("#")),
        delimiter="\t",
        quoting=csv.QUOTE_NONE,
    )
    if reader.fieldnames is None:
        raise SystemExit(f"ERROR: metadata file has no header: {source}")
    fieldnames = [name.strip() for name in reader.fieldnames]
    lowered = [name.lower() for name in fieldnames]
    if "protein_id" not in lowered:
        raise SystemExit(
            f"ERROR: metadata file is missing required column 'protein_id': "
            f"{source}; found {fieldnames}"
        )
    if len(lowered) != len(set(lowered)):
        duplicates = sorted({name for name in lowered if lowered.count(name) > 1})
        raise SystemExit(
            f"ERROR: metadata file has duplicate column names: "
            f"{source}; duplicates={duplicates}"
        )
    rows = reader

    target.parent.mkdir(parents=True, exist_ok=True)
    target_exists = target.is_file() and target.stat().st_size > 0
    if target_exists:
        with target.open("r", encoding="utf-8", errors="replace") as current:
            target_header = current.readline().rstrip("\r\n").split("\t")
        if target_header != canonical:
            raise SystemExit(
                f"ERROR: existing combined metadata has a non-canonical header: "
                f"{target}; found {target_header}; expected {canonical}"
            )

    mode = "a" if target_exists else "w"
    with target.open(mode, encoding="utf-8", newline="") as out:
        writer = csv.DictWriter(
            out,
            fieldnames=canonical,
            delimiter="\t",
            lineterminator="\n",
            quoting=csv.QUOTE_NONE,
            extrasaction="ignore",
        )
        if not target_exists:
            writer.writeheader()
        for row in rows:
            normalized = {
                (key or "").strip().lower(): (value or "").strip()
                for key, value in row.items()
                if key is not None
            }
            writer.writerow({column: normalized.get(column, "") for column in canonical})
PYMETA
}

# ── Promote staged outputs into config-expected paths ─────────────────────────

if [[ "$HAS_FETCH_ROWS" == "true" ]]; then
    echo "[INFO] [$(ts)] Combined fetch list written: ${COMBINED_FETCH_LIST}" >&2

    # Route rows by their fetch_source column (per-row column wins, else
    # ncbi), mirroring how the mod_II protein-accession list is partitioned.
    # Previously mod_I ignored this column entirely, so a uniprot-tagged row
    # inside an otherwise-NCBI list was silently fetched from NCBI instead of
    # UniProt. A per-registry-row options flag (fetch_source=uniprot) still
    # overrides this at the whole-source level -- those rows are deferred
    # before they ever reach this combined list (see the staging step above).
    python3 mod_I_partition_proteome_sources.py \
        -i "$COMBINED_FETCH_LIST" \
        --default-source ncbi \
        --ncbi-out "$COMBINED_FETCH_LIST_NCBI" \
        --uniprot-out "$COMBINED_FETCH_LIST_UNIPROT"

    # Override the Species_fetch_list config var so mod_I_runFetch_Species_INFO.sh
    # and mod_I_runDownload_UserSpecified_protein_FASTA.sh pick it up automatically.
    # Must be exported so subshells invoked via 'bash script.sh' inherit it.
    # Only the NCBI partition flows through that path now; any uniprot-tagged
    # rows are routed to mod_I_fetchUniProtProteomes.py -i further down.
    Species_fetch_list="$COMBINED_FETCH_LIST_NCBI"
    export Species_fetch_list
    user_specified_proteins=true

    if [[ -s "$COMBINED_FETCH_LIST_UNIPROT" ]] && [[ "$(grep -c . "$COMBINED_FETCH_LIST_UNIPROT")" -gt 1 ]]; then
        HAS_FETCH_UNIPROT_ROWS=true
    fi
    export HAS_FETCH_UNIPROT_ROWS
fi

if [[ "$HAS_FASTA_ROWS" == "true" ]]; then
    echo "[INFO] [$(ts)] Combined user-supplied FASTA written: ${COMBINED_FASTA}" >&2
    mkdir -p "$(dirname "$user_fasta_file")"
    cp -f "$COMBINED_FASTA" "$user_fasta_file"
    cp -f "$COMBINED_META"  "$user_metadata_file"
    user_supplied_proteins=true
fi

if [[ "$HAS_PROT_ROWS" == "true" ]]; then
    echo "[INFO] [$(ts)] Protein accession list ready: ${COMBINED_PROT_LIST}" >&2

    PROT_FASTA="${PROT_STAGE}/combined_proteins.fasta"
    PROT_FETCH_REPORT="${PROT_STAGE}/combined_proteins.fetch_report.tsv"

    # Route the combined protein list by per-row fetch_source (column wins,
    # else default ncbi): NCBI -> .ncbi, UniProt -> .uniprot, UniParc -> .uniparc,
    # PDB -> .pdb, anything still unimplemented -> .pending.
    python3 mod_II_partition_protein_sources.py \
        -i "$COMBINED_PROT_LIST" \
        --default-source ncbi \
        --ncbi-out "$COMBINED_PROT_LIST_NCBI" \
        --uniprot-out "$COMBINED_UNIPROT_LIST" \
        --uniparc-out "$COMBINED_UNIPARC_LIST" \
        --pdb-out "$COMBINED_PDB_LIST" \
        --pending-out "$COMBINED_PENDING_LIST"

    # Only run the NCBI exact-ID efetch if the NCBI partition actually has
    # rows. A mod_II accession list can be entirely uniprot/uniparc/pdb
    # (per-row fetch_source column) with zero NCBI accessions -- previously
    # mod_II_fetchProteinAccessions_batch.py aborted the whole build in that
    # case ("no accessions to fetch"), even though the non-NCBI accessions
    # were valid and get fetched below via the *_fromlist blocks. Any source
    # for module II should be enough to proceed; NCBI is not required.
    if [[ -s "$COMBINED_PROT_LIST_NCBI" ]] && [[ "$(grep -c . "$COMBINED_PROT_LIST_NCBI")" -gt 1 ]]; then
        echo "[INFO] [$(ts)] Fetching protein sequences via exact-ID efetch..." >&2

        # Determine sleep rate. NCBI allows 10 req/s with a key, 3 req/s
        # without; ncbi_sleep_api_key / ncbi_sleep_no_key in the config
        # override the shipped 0.11 / 0.4 defaults.
        if [[ -n "${ncbi_api_key:-}" ]]; then
            PROT_SLEEP="${ncbi_sleep_api_key:-0.11}"
        else
            PROT_SLEEP="${ncbi_sleep_no_key:-0.4}"
        fi

        # IDs per efetch POST. 200 is NCBI's documented ceiling for the
        # protein db and remains the default; lower it via ncbi_batch_size
        # in the config if NCBI starts rejecting large POSTs.
        PROT_BATCH="${ncbi_batch_size:-200}"

        EMAIL_FLAG=(); KEY_FLAG=()
        [[ -n "${ncbi_email:-}"   ]] && EMAIL_FLAG=( --email   "$ncbi_email"   )
        [[ -n "${ncbi_api_key:-}" ]] && KEY_FLAG=(   --api_key "$ncbi_api_key" )

        python3 mod_II_fetchProteinAccessions_batch.py \
            -i  "$COMBINED_PROT_LIST_NCBI" \
            -o  "$PROT_FASTA" \
            --report "$PROT_FETCH_REPORT" \
            --sleep "$PROT_SLEEP" \
            --batch-size "$PROT_BATCH" \
            "${EMAIL_FLAG[@]}" \
            "${KEY_FLAG[@]}" \
            "${ALLOW_MISSING_FLAG[@]}"

        # Append exact-ID rewritten FASTA to user-supplied build inputs.
        mkdir -p "$(dirname "$user_fasta_file")"
        cat "$PROT_FASTA" >> "$user_fasta_file"

        user_supplied_proteins=true
        echo "[INFO] [$(ts)] Protein accession exact fetch complete." >&2
        echo "[INFO] [$(ts)] Protein fetch report: ${PROT_FETCH_REPORT}" >&2
    else
        echo "[INFO] [$(ts)] No NCBI-type accessions in combined protein list; skipping exact-ID efetch (any uniprot/uniparc/pdb per-row accessions are fetched below)." >&2
    fi

    # Metadata for ALL non-deferred mod_II rows (any backend) lives in
    # COMBINED_PROT_META regardless of which partition their accessions
    # ended up in, so merge it independently of whether the NCBI leg above ran.
    if [[ -f "$COMBINED_PROT_META" ]]; then
        mkdir -p "$(dirname "$user_metadata_file")"
        append_metadata_normalized "$COMBINED_PROT_META" "$user_metadata_file"
    fi
fi


# UniProt protein accessions                                            ← new
UP_PROT_FASTA="${PROT_STAGE}/uniprot_proteins.fasta"
UP_PROT_META="${PROT_STAGE}/uniprot_proteins.metadata.tsv"
python3 mod_II_fetchUniProtAccessions.py \
    --plan     "$SOURCE_PLAN_NORMALIZED"        \
    -o         "$UP_PROT_FASTA"                 \
    --metadata "$UP_PROT_META"                  \
    --report   "${PROT_STAGE}/uniprot_proteins.fetch_report.tsv" \
    --sleep    "${UNIPROT_SLEEP:-0.4}"          \
    --include-trembl  "${uniprot_include_trembl:-false}"  \
    --include-isoform "${uniprot_include_isoform:-false}" \
    "${ALLOW_MISSING_FLAG[@]}"

if [[ -s "$UP_PROT_FASTA" ]]; then
    mkdir -p "$(dirname "$user_fasta_file")"
    cat "$UP_PROT_FASTA" >> "$user_fasta_file"
    append_metadata_normalized "$UP_PROT_META" "$user_metadata_file"
    user_supplied_proteins=true
fi

# Per-row UniProt accessions partitioned out of mixed protein files (fetch_source column)
if [[ -s "$COMBINED_UNIPROT_LIST" ]] && [[ "$(grep -c . "$COMBINED_UNIPROT_LIST")" -gt 1 ]]; then
    UPL_FASTA="${PROT_STAGE}/uniprot_proteins_fromlist.fasta"
    UPL_META="${PROT_STAGE}/uniprot_proteins_fromlist.metadata.tsv"
    python3 mod_II_fetchUniProtAccessions.py \
        -i         "$COMBINED_UNIPROT_LIST" \
        -o         "$UPL_FASTA" \
        --metadata "$UPL_META" \
        --report   "${PROT_STAGE}/uniprot_proteins_fromlist.fetch_report.tsv" \
        --sleep    "${UNIPROT_SLEEP:-0.4}" \
        "${ALLOW_MISSING_FLAG[@]}"
    if [[ -s "$UPL_FASTA" ]]; then
        mkdir -p "$(dirname "$user_fasta_file")"
        cat "$UPL_FASTA" >> "$user_fasta_file"
        append_metadata_normalized "$UPL_META" "$user_metadata_file"
        user_supplied_proteins=true
    fi
fi

# UniParc protein accessions (plan-declared sources)                    ← new
UPARC_PROT_FASTA="${PROT_STAGE}/uniparc_proteins.fasta"
UPARC_PROT_META="${PROT_STAGE}/uniparc_proteins.metadata.tsv"
python3 mod_II_fetchUniParcAccessions.py \
    --plan          "$SOURCE_PLAN_NORMALIZED"        \
    -o              "$UPARC_PROT_FASTA"              \
    --metadata      "$UPARC_PROT_META"               \
    --report        "${PROT_STAGE}/uniparc_proteins.fetch_report.tsv" \
    --organism-pick "$organism_pick"                 \
    --sleep         "${UNIPROT_SLEEP:-0.4}"          \
    "${ALLOW_MISSING_FLAG[@]}"

if [[ -s "$UPARC_PROT_FASTA" ]]; then
    mkdir -p "$(dirname "$user_fasta_file")"
    cat "$UPARC_PROT_FASTA" >> "$user_fasta_file"
    append_metadata_normalized "$UPARC_PROT_META" "$user_metadata_file"
    user_supplied_proteins=true
fi

# Per-row UniParc accessions partitioned out of mixed protein files (fetch_source column)
if [[ -s "$COMBINED_UNIPARC_LIST" ]] && [[ "$(grep -c . "$COMBINED_UNIPARC_LIST")" -gt 1 ]]; then
    UPARCL_FASTA="${PROT_STAGE}/uniparc_proteins_fromlist.fasta"
    UPARCL_META="${PROT_STAGE}/uniparc_proteins_fromlist.metadata.tsv"
    python3 mod_II_fetchUniParcAccessions.py \
        -i              "$COMBINED_UNIPARC_LIST" \
        -o              "$UPARCL_FASTA" \
        --metadata      "$UPARCL_META" \
        --report        "${PROT_STAGE}/uniparc_proteins_fromlist.fetch_report.tsv" \
        --organism-pick "$organism_pick" \
        --sleep         "${UNIPROT_SLEEP:-0.4}" \
        "${ALLOW_MISSING_FLAG[@]}"
    if [[ -s "$UPARCL_FASTA" ]]; then
        mkdir -p "$(dirname "$user_fasta_file")"
        cat "$UPARCL_FASTA" >> "$user_fasta_file"
        append_metadata_normalized "$UPARCL_META" "$user_metadata_file"
        user_supplied_proteins=true
    fi
fi

# PDB protein entities (plan-declared sources)                          ← new
PDB_PROT_FASTA="${PROT_STAGE}/pdb_proteins.fasta"
PDB_PROT_META="${PROT_STAGE}/pdb_proteins.metadata.tsv"
python3 mod_II_fetchPDBAccessions.py \
    --plan     "$SOURCE_PLAN_NORMALIZED"        \
    -o         "$PDB_PROT_FASTA"                \
    --metadata "$PDB_PROT_META"                 \
    --report   "${PROT_STAGE}/pdb_proteins.fetch_report.tsv" \
    --sleep    "${PDB_SLEEP:-0.4}"              \
    "${ALLOW_MISSING_FLAG[@]}"

if [[ -s "$PDB_PROT_FASTA" ]]; then
    mkdir -p "$(dirname "$user_fasta_file")"
    cat "$PDB_PROT_FASTA" >> "$user_fasta_file"
    append_metadata_normalized "$PDB_PROT_META" "$user_metadata_file"
    user_supplied_proteins=true
fi

# Per-row PDB accessions partitioned out of mixed protein files (fetch_source column)
if [[ -s "$COMBINED_PDB_LIST" ]] && [[ "$(grep -c . "$COMBINED_PDB_LIST")" -gt 1 ]]; then
    PDBL_FASTA="${PROT_STAGE}/pdb_proteins_fromlist.fasta"
    PDBL_META="${PROT_STAGE}/pdb_proteins_fromlist.metadata.tsv"
    python3 mod_II_fetchPDBAccessions.py \
        -i         "$COMBINED_PDB_LIST" \
        -o         "$PDBL_FASTA" \
        --metadata "$PDBL_META" \
        --report   "${PROT_STAGE}/pdb_proteins_fromlist.fetch_report.tsv" \
        --sleep    "${PDB_SLEEP:-0.4}" \
        "${ALLOW_MISSING_FLAG[@]}"
    if [[ -s "$PDBL_FASTA" ]]; then
        mkdir -p "$(dirname "$user_fasta_file")"
        cat "$PDBL_FASTA" >> "$user_fasta_file"
        append_metadata_normalized "$PDBL_META" "$user_metadata_file"
        user_supplied_proteins=true
    fi
fi

module_banner "Module IV-A: supported collection — human UniProt proteins"

# ═════════════════════════════════════════════════════════════════════════════
# HUMAN PROTEINS
#
# use_human        = true | false  — include human proteins at all
# use_human_download = true | false  — download fresh vs. use cached file
#
# use_human=false  : human proteins excluded entirely from the build
# use_human=true, use_human_download=true  : download fresh from UniProt
# use_human=true, use_human_download=false : use existing file at $human_fasta_file
# ═════════════════════════════════════════════════════════════════════════════

HUMAN_FLAG=()

if [[ "${use_human:-true}" == "true" ]]; then
    if [[ "${use_human_download:-true}" == "true" ]]; then
        echo "[INFO] [$(ts)] Fetching human proteins from UniProt (${human_protein_set})..." >&2
        bash mod_IV_runFetch_Human_Proteins.sh
    else
        echo "[INFO] [$(ts)] Skipping human protein download; using existing file." >&2
        [[ -f "$human_fasta_file" ]] || {
            echo "ERROR: [$(ts)] use_human_download=false but file not found: ${human_fasta_file}" >&2
            exit 1
        }
    fi
    HUMAN_FLAG=( --human_fasta "$human_fasta_file" )
else
    echo "[INFO] [$(ts)] Human proteins excluded from this build (use_human=false)." >&2
fi

module_banner "Module IV-B: supported collection — HOMD proteomes"

# ═════════════════════════════════════════════════════════════════════════════
# HOMD
#
# use_homd = true | false
#
# If true:
#   1. Download fresh HOMD taxonomy table
#   2. Filter by body site and download matching proteomes
#      - homd_sites set   -> filtered taxonomy + filtered proteomes
#      - homd_sites empty -> full taxonomy    + all proteomes
#      - homd_rank set    -> dereplicate survivors to one representative per
#                            species|genus|family (fewest contigs, ties -> size),
#                            applied AFTER the body-site filter
#   3. Pass HOMD args to build
#
# If false: all HOMD steps are skipped entirely
# ═════════════════════════════════════════════════════════════════════════════

HOMD_ARGS=()

if [[ "${use_homd:-true}" == "true" ]]; then

    if [[ "${use_homd_download:-true}" == "true" ]]; then
        echo "[INFO] [$(ts)] Fetching HOMD taxonomy table..." >&2
        bash mod_IV_runFetch_HOMD_Taxonomy.sh

        echo "[INFO] [$(ts)] Fetching HOMD proteomes..." >&2
        bash mod_IV_runFetch_HOMD_Proteomes.sh
    else
        echo "[INFO] [$(ts)] Skipping HOMD download; using cached files." >&2
        [[ -f "$homd_fasta_file" ]] || {
            echo "ERROR: [$(ts)] use_homd_download=false but homd_fasta_file not found: ${homd_fasta_file}" >&2; exit 1
        }
        [[ -f "$gca_info_file" ]] || {
            echo "ERROR: [$(ts)] use_homd_download=false but gca_info_file not found: ${gca_info_file}" >&2; exit 1
        }
        if [[ -n "${homd_sites:-}" && ! -f "$homd_filtered_tax" ]]; then
            echo "[INFO] [$(ts)] Filtered taxonomy not cached; running mod_IV_homd_filter.py..." >&2
            mkdir -p "$(dirname "$homd_filtered_tax")"
            python3 mod_IV_homd_filter.py \
                --file   "$homd_tax_file" \
                --sites  $homd_sites \
                --output "$homd_filtered_tax"
        fi
    fi

    if [[ -n "${homd_sites:-}" ]]; then
        HOMD_TAX_TO_USE="$homd_filtered_tax"
    else
        HOMD_TAX_TO_USE="$homd_tax_file"
    fi

    HOMD_ARGS=(
        --homd_fasta "$homd_fasta_file"
        --gca_info   "$gca_info_file"
        --homd_tax   "$HOMD_TAX_TO_USE"
    )

else
    echo "[INFO] [$(ts)] HOMD excluded from this build (use_homd=false)." >&2
fi

module_banner "Module I: fetch user-specified proteomes from NCBI"

# ═════════════════════════════════════════════════════════════════════════════
# USER-SPECIFIED SPECIES
# (mod_I proteome rows from the normalized source plan)
#
# If true:
#   1. Fetch species taxonomy info via mod_I_Fetch_UserSpecifiedSpecies_INFO.sh
#   2. Download proteomes via mod_I_download_UserSpecified_protein_FASTA.sh
#   3. Pass --proteomes_manifest to build
# ═════════════════════════════════════════════════════════════════════════════

MANIFEST_FLAG=()

if [[ "${user_specified_proteins:-false}" == "true" ]]; then
    echo "[INFO] [$(ts)] Fetching species taxonomy info..." >&2
    bash mod_I_runFetch_Species_INFO.sh
    echo "[INFO] [$(ts)] Downloading species proteomes..." >&2
    bash mod_I_runDownload_UserSpecified_protein_FASTA.sh
    [[ -f "$proteomes_manifest" ]] || {
        echo "ERROR: [$(ts)] proteomes_manifest not found: ${proteomes_manifest}" >&2; exit 1; }
else
    echo "[INFO] [$(ts)] No NCBI user-specified proteomes in this build." >&2
fi

# UniProt proteomes / pan-proteomes — appends to the SAME manifest      ← new
mkdir -p "$proteomes_dir"
python3 mod_I_fetchUniProtProteomes.py \
    --plan          "$SOURCE_PLAN_NORMALIZED" \
    --proteomes-dir "$proteomes_dir"          \
    --manifest      "$proteomes_manifest"     \
    --proteome-type "$proteome_type"          \
    --sleep         "${UNIPROT_SLEEP:-0.4}"

# Per-row UniProt proteome rows partitioned out of an otherwise-NCBI mod_I
# input_list (fetch_source column) — the flat-list counterpart to the
# --plan call above, mirroring mod_II's *_fromlist pattern.
if [[ "${HAS_FETCH_UNIPROT_ROWS:-false}" == "true" ]]; then
    echo "[INFO] [$(ts)] Fetching per-row UniProt-tagged proteomes..." >&2
    python3 mod_I_fetchUniProtProteomes.py \
        -i              "$COMBINED_FETCH_LIST_UNIPROT" \
        --proteomes-dir "$proteomes_dir"                \
        --manifest      "$proteomes_manifest"           \
        --proteome-type "$proteome_type"                \
        --sleep         "${UNIPROT_SLEEP:-0.4}"
fi

[[ -f "$proteomes_manifest" ]] && MANIFEST_FLAG=( --proteomes_manifest "$proteomes_manifest" )







module_banner "Modules II–III: finalize user-supplied and accession-fetched proteins"

# ═════════════════════════════════════════════════════════════════════════════
# USER-SUPPLIED PROTEINS
# (mod_III user-FASTA and mod_II protein rows from the normalized source plan)
# ═════════════════════════════════════════════════════════════════════════════

USER_FASTA_FLAG=()
USER_META_FLAG=()

if [[ "${user_supplied_proteins:-false}" == "true" ]]; then
    [[ -n "${user_fasta_file:-}" ]]    && USER_FASTA_FLAG=( --user_fasta     "$user_fasta_file"    )
    [[ -n "${user_metadata_file:-}" ]] && USER_META_FLAG=(  --user_metadata  "$user_metadata_file" )
else
    echo "[INFO] [$(ts)] User-supplied proteins excluded from this build." >&2
fi

module_banner "Module IV-C: supported collection — cRAP contaminants"

# cRAP is prepared immediately before database construction so the exact
# run-used contaminant FASTA is captured in run_data.

# ═════════════════════════════════════════════════════════════════════════════
# BUILD DATABASE
# ═════════════════════════════════════════════════════════════════════════════

NCBI_EMAIL_FLAG=()
NCBI_KEY_FLAG=()
# ── cRAP contaminants — download if needed ───────────────────────────────────
CRAP_FLAG=()

if [[ -n "${crap_set:-}" && "${crap_set}" != "none" ]]; then
    # cRAP is run-scoped. Raw downloads are temporary; the exact combined FASTA
    # used by this build is written flat into run_data/.
    crap_dir="${crap_dir:-${refDIR}}"
    CRAP_TMP="${TMPDIR}/crap_downloads"
    mkdir -p "$crap_dir" "$CRAP_TMP"

    ZENODO="https://zenodo.org/records/15115102/files"
    CRAP_COMBINED="${crap_dir}/crap_combined_${crap_set}.fasta"

    _dl_crap() {
        local name="$1"
        local url="$2"
        local dest="${CRAP_TMP}/${name}"

        echo "[INFO] [$(ts)] Downloading ${name} from Zenodo..." >&2
        curl --http1.1 -gSL --retry 3 --retry-delay 2 \
            --connect-timeout 30 --max-time 300 \
            "${url}?download=1" -o "$dest"
        [[ -s "$dest" ]] || { echo "ERROR: [$(ts)] failed to download ${name}" >&2; exit 1; }
    }

    > "$CRAP_COMBINED"
    case "${crap_set}" in
        ccp)
            _dl_crap "crap_ccp.fasta"         "${ZENODO}/crap_ccp.fasta"
            cat "${CRAP_TMP}/crap_ccp.fasta" >> "$CRAP_COMBINED" ;;
        gpm)
            _dl_crap "crap_gpm.fasta"         "${ZENODO}/crap_gpm.fasta"
            cat "${CRAP_TMP}/crap_gpm.fasta" >> "$CRAP_COMBINED" ;;
        maxquant)
            _dl_crap "crap_maxquant.fasta.gz" "${ZENODO}/crap_maxquant.fasta.gz"
            gunzip -c "${CRAP_TMP}/crap_maxquant.fasta.gz" >> "$CRAP_COMBINED" ;;
        all)
            _dl_crap "crap_ccp.fasta"         "${ZENODO}/crap_ccp.fasta"
            _dl_crap "crap_gpm.fasta"         "${ZENODO}/crap_gpm.fasta"
            _dl_crap "crap_maxquant.fasta.gz" "${ZENODO}/crap_maxquant.fasta.gz"
            cat "${CRAP_TMP}/crap_ccp.fasta"               >> "$CRAP_COMBINED"
            cat "${CRAP_TMP}/crap_gpm.fasta"               >> "$CRAP_COMBINED"
            gunzip -c "${CRAP_TMP}/crap_maxquant.fasta.gz" >> "$CRAP_COMBINED" ;;
        *)
            echo "ERROR: [$(ts)] unknown crap_set='${crap_set}'. Choose: ccp|gpm|maxquant|all" >&2
            exit 1 ;;
    esac

    crap_fasta="$CRAP_COMBINED"
    CRAP_FLAG=( --crap_fasta "$crap_fasta" )
    CRAP_COUNT="$(grep -c '^>' "$crap_fasta" || true)"
    echo "[INFO] [$(ts)] cRAP set='${crap_set}': ${CRAP_COUNT} sequences" >&2
    echo "[INFO] [$(ts)] cRAP combined FASTA for this run: ${crap_fasta}" >&2
fi

NCBI_EMAIL_FLAG=()
NCBI_KEY_FLAG=()
[[ -n "${ncbi_email:-}"     ]] && NCBI_EMAIL_FLAG=( --ncbi_email   "$ncbi_email"     )
[[ -n "${ncbi_api_key:-}"   ]] && NCBI_KEY_FLAG=(   --ncbi_api_key "$ncbi_api_key"   )

# Hand any skipped/missing accessions (from --allow-missing fetchers) to build_db
# so they are surfaced in the build summary's failed/skipped section.
FAILED_ACC_FLAG=()
[[ -n "${FETCH_MISSING:-}" && -f "${FETCH_MISSING}" ]] && FAILED_ACC_FLAG=( --failed_accessions "$FETCH_MISSING" )

phase_banner "PHASE 2: Build harmonized maniFasta database"

SOURCE_PLAN_FLAG=()
[[ -n "${SOURCE_PLAN_NORMALIZED:-}" && -f "${SOURCE_PLAN_NORMALIZED}" ]] && \
    SOURCE_PLAN_FLAG=( --source_plan "${SOURCE_PLAN_NORMALIZED}" )

python3 build_db.py \
    "${HOMD_ARGS[@]}"                         \
    "${HUMAN_FLAG[@]}"                        \
    "${CRAP_FLAG[@]}"                         \
    "${MANIFEST_FLAG[@]}"                     \
    "${USER_FASTA_FLAG[@]}"                   \
    "${USER_META_FLAG[@]}"                    \
    "${NCBI_EMAIL_FLAG[@]}"                   \
    "${NCBI_KEY_FLAG[@]}"                     \
    "${FAILED_ACC_FLAG[@]}"                   \
    "${SOURCE_PLAN_FLAG[@]}"                  \
    --out_fasta          "$db_fasta_output"  \
    --out_map            "$db_map_output"    \
    --run_label          "${RUN_LABEL}"      \
    --run_dir            "${RUN_DIR}"        \
    --db_prefix          "${db_prefix:-maniFastaDB}"  \
    --out_summary        "${RUN_DIR}/build_summary.tsv"

# ═════════════════════════════════════════════════════════════════════════════
# PHASE 2.5: Taxonomic lineage enrichment (post-build annotation)
#
# Resolve each protein's lineage_taxid (the manifest's rollup taxid: defaults to
# ncbi_taxid, or a registry lineage_taxid_override) to a full ranked NCBI
# lineage and append tax_* columns to a copy of the manifest. The enriched file
# feeds taxonomy summaries / the sunburst graphic. The original manifest is left
# untouched.
#
# Config knobs (00.setup/000.maniFasta.config):
#   add_lineage         true|false  enable this step (default true)
#   lineage_cache_file  path         persistent lineage cache reused across runs
#                                    (default ${mainDIR}/taxid_lineage_cache.tsv)
# Reuses ncbi_email / ncbi_api_key for the NCBI Taxonomy rate limit. This step is
# NON-FATAL: a failure warns but does not fail an otherwise-successful build.
# ═════════════════════════════════════════════════════════════════════════════
phase_banner "PHASE 2.5: Taxonomic lineage enrichment"

if [[ "${add_lineage:-true}" == "true" ]]; then
    LINEAGE_MANIFEST="${RUN_DIR}/${RUN_LABEL}.manifest.lineage.tsv"
    LINEAGE_CACHE="${lineage_cache_file:-${mainDIR}/taxid_lineage_cache.tsv}"

    LINEAGE_EMAIL_FLAG=(); LINEAGE_KEY_FLAG=(); LINEAGE_CACHE_IN=()
    [[ -n "${ncbi_email:-}"   ]] && LINEAGE_EMAIL_FLAG=( --email   "$ncbi_email" )
    [[ -n "${ncbi_api_key:-}" ]] && LINEAGE_KEY_FLAG=(   --api-key "$ncbi_api_key" )
    [[ -f "$LINEAGE_CACHE"    ]] && LINEAGE_CACHE_IN=(   --lineage-cache "$LINEAGE_CACHE" )

    # If no email is configured and live fetch isn't explicitly allowed, stay
    # offline: enrich only from the cache, label the rest Unclassified.
    NO_FETCH_FLAG=()
    if [[ -z "${ncbi_email:-}" && "${lineage_allow_fetch:-true}" != "true" ]]; then
        NO_FETCH_FLAG=( --no-fetch )
        echo "[INFO] [$(ts)] No ncbi_email and lineage_allow_fetch!=true: lineage from cache only." >&2
    fi

    echo "[INFO] [$(ts)] Enriching manifest taxonomy (lineage cache: ${LINEAGE_CACHE})" >&2
    if python3 add_lineage_to_metadata.py \
            --metadata    "$db_map_output"   \
            --out         "$LINEAGE_MANIFEST" \
            --taxid-col   lineage_taxid       \
            --write-cache "$LINEAGE_CACHE"    \
            "${LINEAGE_CACHE_IN[@]}"          \
            "${LINEAGE_EMAIL_FLAG[@]}"        \
            "${LINEAGE_KEY_FLAG[@]}"          \
            "${NO_FETCH_FLAG[@]}"; then
        echo "[INFO] [$(ts)] Lineage-enriched manifest: ${LINEAGE_MANIFEST}" >&2
    else
        echo "[WARN] [$(ts)] Lineage enrichment failed; the build manifest is unaffected." >&2
        echo "[WARN] [$(ts)] Re-run add_lineage_to_metadata.py on ${db_map_output} to retry." >&2
    fi
else
    echo "[INFO] [$(ts)] add_lineage=false: skipping taxonomic lineage enrichment." >&2
fi

# ═════════════════════════════════════════════════════════════════════════════
# PHASE 2.6: Static taxonomy sunburst (optional, headless matplotlib)
#
# Renders a publication-ready sunburst (PNG + SVG) from the lineage-enriched
# manifest produced in PHASE 2.5. Runs on a compute node with no display.
# Config knobs:
#   plot_taxonomy              true|false   render the figure (default true)
#   taxonomy_plot_color_by     phylum|domain_or_realm   (default phylum)
#   taxonomy_plot_label_ranks  comma list   ranks to label (default phylum,species)
#   taxonomy_plot_html         true|false   also emit interactive HTML (default true)
# NON-FATAL and only runs if PHASE 2.5 produced an enriched manifest. Requires
# matplotlib in the active python3 environment; if absent, this step warns and
# the build still succeeds.
# ═════════════════════════════════════════════════════════════════════════════
if [[ "${plot_taxonomy:-true}" == "true" && -n "${LINEAGE_MANIFEST:-}" && -f "${LINEAGE_MANIFEST:-/nonexistent}" ]]; then
    phase_banner "PHASE 2.6: Static taxonomy sunburst"
    TAX_PLOT_SVG="${RUN_DIR}/${RUN_LABEL}.taxonomy_sunburst.svg"
    echo "[INFO] [$(ts)] Rendering taxonomy sunburst from ${LINEAGE_MANIFEST}" >&2

    # SVG is the stdlib backend (no matplotlib/numpy) and always renders.
    python3 plot_taxonomy.py \
        --metadata    "$LINEAGE_MANIFEST" \
        --out         "$TAX_PLOT_SVG" \
        --color-by    "${taxonomy_plot_color_by:-phylum}" \
        --label-ranks "${taxonomy_plot_label_ranks:-phylum,species}" \
        --label-style "${taxonomy_plot_label_style:-callout}" \
        --subtitle    "${RUN_LABEL}" \
        && echo "[INFO] [$(ts)] Taxonomy sunburst (SVG): ${TAX_PLOT_SVG}" >&2 \
        || echo "[WARN] [$(ts)] taxonomy SVG render failed (non-fatal)." >&2

    # PNG is optional and needs matplotlib in the environment. Off by default so
    # clusters without matplotlib don't log an error every build.
    if [[ "${taxonomy_plot_png:-false}" == "true" ]]; then
        TAX_PLOT_PNG="${RUN_DIR}/${RUN_LABEL}.taxonomy_sunburst.png"
        python3 plot_taxonomy.py \
            --metadata    "$LINEAGE_MANIFEST" \
            --out         "$TAX_PLOT_PNG" \
            --color-by    "${taxonomy_plot_color_by:-phylum}" \
            --label-ranks "${taxonomy_plot_label_ranks:-phylum,species}" \
            --label-style "${taxonomy_plot_label_style:-callout}" \
            --subtitle    "${RUN_LABEL}" \
            && echo "[INFO] [$(ts)] Taxonomy sunburst (PNG): ${TAX_PLOT_PNG}" >&2 \
            || echo "[WARN] [$(ts)] taxonomy PNG render failed (non-fatal; matplotlib required for PNG)." >&2
    fi

    # Interactive HTML sunburst. stdlib only (imports plot_taxonomy.py for the
    # tree/colour logic) so it always renders; hover tooltips carry the full
    # taxon labels the static figure has to truncate.
    if [[ "${taxonomy_plot_html:-true}" == "true" ]]; then
        TAX_PLOT_HTML="${RUN_DIR}/${RUN_LABEL}.taxonomy_sunburst.html"
        python3 plot_taxonomy_html.py \
            --metadata      "$LINEAGE_MANIFEST" \
            --out           "$TAX_PLOT_HTML" \
            --plot-taxonomy "${scriptsDIR}/plot_taxonomy.py" \
            --color-by      "${taxonomy_plot_color_by:-phylum}" \
            --title         "${RUN_LABEL}" \
            --subtitle      "Taxonomic composition of the search database" \
            && echo "[INFO] [$(ts)] Taxonomy sunburst (HTML): ${TAX_PLOT_HTML}" >&2 \
            || echo "[WARN] [$(ts)] taxonomy HTML render failed (non-fatal)." >&2
    fi
elif [[ "${plot_taxonomy:-true}" == "true" ]]; then
    echo "[INFO] [$(ts)] plot_taxonomy=true but no enriched manifest available; skipping sunburst." >&2
fi

# ═════════════════════════════════════════════════════════════════════════════
# PHASE 2.7: Database methods summary (auto-generated reporting)
#
# Reads this run's build artifacts (final FASTA, manifest / lineage manifest,
# normalized source plan, HOMD provenance, supported-collections env) and writes
# a methods-style summary next to build_summary.tsv:
#   ${RUN_LABEL}.database_methods.md
# It addresses standard search-database reporting requirements: sequence count,
# build/download date, per-source origin, taxon selection, redundancy handling
# (incl. HOMD representative selection), functional/taxonomic annotation, and
# contaminant incorporation. It only READS artifacts; it computes/fetches nothing.
#
# The normalized plan and supported-collections env still live in $TMPDIR at this
# point (they are packaged into run_logs/ later by the EXIT trap), so $TMPDIR is
# passed explicitly via --tmp-dir.
#
# Config knob (00.setup/000.maniFasta.config):
#   write_db_summary  true|false  enable this step (default true)
# NON-FATAL: a failure warns but never fails an otherwise-successful build.
# ═════════════════════════════════════════════════════════════════════════════
if [[ "${write_db_summary:-true}" == "true" ]]; then
    phase_banner "PHASE 2.7: Database methods summary"
    if python3 generate_db_summary.py \
            --run-dir   "$RUN_DIR" \
            --run-label "$RUN_LABEL" \
            --tmp-dir   "$TMPDIR"; then
        echo "[INFO] [$(ts)] Methods summary: ${RUN_DIR}/${RUN_LABEL}.database_methods.md" >&2
    else
        echo "[WARN] [$(ts)] Database methods summary failed (non-fatal)." >&2
        echo "[WARN] [$(ts)] Re-run manually: python3 generate_db_summary.py --run-dir ${RUN_DIR}" >&2
    fi
else
    echo "[INFO] [$(ts)] write_db_summary=false: skipping database methods summary." >&2
fi
