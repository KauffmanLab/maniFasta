#!/usr/bin/env bash

## Originated by Christopher Handelmann

## HOMD proteome fetcher — PROKKA source (GCA-named faa files).
##
## CHANGED: instead of extracting GCA accessions and re-downloading proteins
## from NCBI via ncbi-datasets-cli, this version downloads the PROKKA-annotated
## proteomes that HOMD itself distributes, one file per GenBank assembly:
##
##     https://www.homd.org/ftp/genomes/PROKKA/<VERSION>/faa/<GCA>.faa
##
## The <VERSION> component (e.g. V11.03) comes from the HOMD row of
## source_registry.tsv; homd_prokka_base_url in the config stops at .../PROKKA.
##
## Per the HOMD PROKKA README (v11.02), the faa files are named with the GCA
## assembly ID, so selection maps directly onto GCA_ID_info.txt with no SEQF
## indirection. Workflow:
##
##   1. Filter taxonomy by body site (mod_IV_homd_filter.py)   -> HMT IDs
##   2. GCA_ID_info.txt                                  -> GCA -> HMT
##   3. Select GCAs whose HMT passes the site filter
##   2.5 (optional) Harvest contigs+size per candidate from PROKKA summary/
##   2.6 (optional) Dereplicate to ONE representative per taxonomic rank
##   4. Download each <GCA>.faa from HOMD/PROKKA (winners only if dereplicating)
##   5. Prepend the GCA accession to every header so build_db.py's assembly
##      lookup resolves:   >GCA_000006605.1_<PROKKA locus tag> ...
##   6. Concatenate into a single FASTA at $homd_fasta_file
##
## Representative selection (steps 2.5/2.6) is OPT-IN via $homd_rank. When
## $homd_rank is empty the script behaves exactly as before (download every
## body-site match), so this change is non-breaking.
##
## NOTE: PROKKA locus tags contain underscores, so build_db.py must extract the
## assembly via a regex anchored on the leading GC[AF]_<digits>.<ver> rather than
## rsplit('_', 1). See the accompanying build_db.py patch.

set -euo pipefail

source "${Config_file}"
source "${setupDIR}/000.maniFasta.local.env"

# ── PROKKA source configuration ───────────────────────────────────────────────
# homd_prokka_base_url is the UNVERSIONED release root (.../genomes/PROKKA); the
# release version is chosen on the HOMD row of source_registry.tsv via
# 'genomic_refseq_version=' (e.g. V11.03), or auto-detected as the newest release
# when that option is omitted. The planner resolves it to a concrete version and
# exports homd_genomic_refseq_version, so this script and
# mod_IV_runFetch_HOMD_Taxonomy.sh always read the same genome set.
: "${homd_prokka_base_url:=https://www.homd.org/ftp/genomes/PROKKA}"

# Resolved release root: <base>/<version>. Every PROKKA URL below hangs off this,
# so the version lives in exactly one place.
if [[ -z "${homd_genomic_refseq_version:-}" ]]; then
    echo "ERROR: 'homd_genomic_refseq_version' is not set." >&2
    echo "       It is resolved by the source planner from the HOMD row's" >&2
    echo "       genomic_refseq_version= option (or auto-detected). Ensure the" >&2
    echo "       supported-collections env was sourced before running this script." >&2
    exit 1
fi
homd_prokka_release_url="${homd_prokka_base_url%/}/${homd_genomic_refseq_version}"

: "${homd_prokka_faa_url:=${homd_prokka_release_url}/faa}"
# PROKKA per-genome annotation tables (<GCA>.tsv: locus_tag, ftype, length_bp,
# gene, EC_number, COG, product). Used to attach gene symbols to the .faa, which
# carries only the product. Set homd_fetch_gene_names=false to skip these
# downloads and leave gene=NA (original behaviour).
: "${homd_prokka_tsv_url:=${homd_prokka_release_url}/tsv}"
: "${homd_fetch_gene_names:=true}"
# PROKKA per-genome stats live here as <GCA>.txt (contigs:/bases: lines).
: "${homd_prokka_summary_url:=${homd_prokka_release_url}/summary}"
# Parallel faa/summary downloads. HOMD is a single modest server; keep this polite.
: "${homd_prokka_jobs:=4}"

# ── Representative-selection configuration (opt-in) ────────────────────────────
# homd_rank: "" (disabled) | "species" | "genus" | "family".
# When set, the pipeline keeps ONE assembly per group at that rank, choosing the
# most contiguous (fewest contigs; ties -> largest genome size).
: "${homd_rank:=}"

# ── Validate required config vars ─────────────────────────────────────────────

for var in homd_tax_file gca_info_file homd_fasta_file homd_proteomes_dir; do
    [[ -n "${!var:-}" ]] || { echo "ERROR: '$var' is not set in config." >&2; exit 1; }
done

[[ -f "$homd_tax_file"  ]] || { echo "ERROR: homd_tax_file not found: $homd_tax_file"  >&2; exit 1; }
[[ -f "$gca_info_file"  ]] || { echo "ERROR: gca_info_file not found: $gca_info_file"  >&2; exit 1; }

command -v curl    >/dev/null || { echo "ERROR: curl not found"    >&2; exit 1; }
command -v python3 >/dev/null || { echo "ERROR: python3 not found" >&2; exit 1; }

# Validate rank value early (before any downloads) if dereplication is requested.
if [[ -n "${homd_rank:-}" ]]; then
    case "$homd_rank" in
        species|genus|family) : ;;
        *) echo "ERROR: homd_rank must be one of: species, genus, family (got '$homd_rank')" >&2; exit 1 ;;
    esac
fi

: "${TMPDIR:=/tmp}"
mkdir -p "$TMPDIR"
mkdir -p "$(dirname "$homd_fasta_file")"
mkdir -p "$homd_proteomes_dir"

# Run-scoped HOMD directory (…/{run}/run_data/homd). All HOMD mid-step
# provenance — body-site selection, assembly stats, representative report —
# lands here, beside the taxonomy tables and the final homd_fasta_file.
HOMD_DIR="$(dirname "$homd_fasta_file")"

# ── Step 1: Filter taxonomy by body site ──────────────────────────────────────

if [[ -n "${homd_sites:-}" ]]; then
    echo "[INFO] Filtering HOMD taxonomy to site(s): ${homd_sites}" >&2
    mkdir -p "$(dirname "$homd_filtered_tax")"
    python3 mod_IV_homd_filter.py \
        --file   "$homd_tax_file" \
        --sites  $homd_sites \
        --output "$homd_filtered_tax"
    TAX_TO_USE="$homd_filtered_tax"
else
    echo "[INFO] No site filter specified; using full HOMD taxonomy." >&2
    TAX_TO_USE="$homd_tax_file"
fi

# ── Step 2: Select GCA assemblies whose HMT passes the body-site filter ───────
#
# Output: SELECTED_PLAN  (whitespace: "<GCA> <HMT>")  — these are exactly the
# genomes the old NCBI path would have downloaded, now sourced from PROKKA.
SELECTED_PLAN="${TMPDIR}/HOMD_prokka_selected.txt"
PLAN_REPORT="${HOMD_DIR}/HOMD_prokka_selection.tsv"

python3 - "$TAX_TO_USE" "$gca_info_file" "$SELECTED_PLAN" "$PLAN_REPORT" << 'PY'
import sys
import re
from pathlib import Path

tax_file  = Path(sys.argv[1])
gca_file  = Path(sys.argv[2])
plan_out  = Path(sys.argv[3])
report_out = Path(sys.argv[4])

def info(m): print(f"[INFO] {m}", file=sys.stderr)
def die(m):
    print(f"ERROR: {m}", file=sys.stderr); raise SystemExit(1)

# eHOMD changed HMT-ID zero-padding between releases (HMT-047 -> HMT-0047) while
# PROKKA/<version>/GCA_ID_info.txt kept the 3-digit form. Canonicalise both sides
# of the taxonomy <-> GCA_ID_info join so either padding works. Must match
# norm_hmt() in build_db.py and mod_IV_select_representatives.py.
_HMT_RE = re.compile(r"^HMT-?0*([0-9]+)$", re.IGNORECASE)

def norm_hmt(hmt):
    hmt = (hmt or "").strip()
    m = _HMT_RE.match(hmt)
    return f"HMT-{int(m.group(1)):04d}" if m else hmt

# ── Filtered HMT set (homd_filter output: banner line, then header w/ HMT-ID) ──
hmt_set = set()
with tax_file.open(encoding="utf-8") as fh:
    header = None
    idx_hmt = None
    for line in fh:
        line = line.rstrip("\n")
        if not line or line.startswith("HOMD.org Taxon Data::"):
            continue
        cols = line.split("\t")
        if header is None:
            header = [c.strip() for c in cols]
            idx_hmt = next((i for i, h in enumerate(header) if "HMT" in h.upper()), None)
            if idx_hmt is None:
                die(f"no HMT column in taxonomy file: {tax_file} (header={header})")
            continue
        if len(cols) > idx_hmt and cols[idx_hmt].strip():
            hmt_set.add(norm_hmt(cols[idx_hmt].strip()))
if not hmt_set:
    die(f"no HMT IDs parsed from taxonomy: {tax_file}")
info(f"Filtered HMT IDs: {len(hmt_set)}")

# ── GCA -> HMT from GCA_ID_info.txt. Only the two leading single-token columns
# (GCA-ID, HMT-ID) are read, so trailing 'Status'/multi-word 'Strain' columns
# are irrelevant. ──────────────────────────────────────────────────────────────
selected = []
seen = set()
with gca_file.open(encoding="utf-8") as fh:
    for line in fh:
        line = line.strip()
        if not line or line.startswith("GCA-ID"):
            continue
        cols = line.split()
        if len(cols) < 2:
            continue
        gca, hmt = cols[0].strip(), norm_hmt(cols[1].strip())
        if gca and hmt in hmt_set and gca not in seen:
            selected.append((gca, hmt))
            seen.add(gca)

if not selected:
    die("no GCA assemblies matched the requested body site(s).")
info(f"GCA assemblies selected by body-site filter: {len(selected)}")

with plan_out.open("w", encoding="utf-8") as fh:
    for gca, hmt in selected:
        fh.write(f"{gca} {hmt}\n")
with report_out.open("w", encoding="utf-8") as fh:
    fh.write("gca\thmt_id\n")
    for gca, hmt in selected:
        fh.write(f"{gca}\t{hmt}\n")
info(f"Selection report: {report_out}")
PY

PLAN_COUNT="$(wc -l < "$SELECTED_PLAN" | tr -d ' ')"
echo "[INFO] Candidate proteomes after body-site filter: ${PLAN_COUNT}" >&2
[[ "$PLAN_COUNT" -gt 0 ]] || { echo "ERROR: empty selection." >&2; exit 1; }

# ── Steps 2.5 / 2.6: optional dereplication to one representative per rank ─────
#
# DOWNLOAD_PLAN is what Step 3 actually downloads. By default it is the full
# body-site selection (old behavior). If $homd_rank is set, we first harvest
# per-assembly stats, then collapse to one representative per taxonomic group.
DOWNLOAD_PLAN="$SELECTED_PLAN"

if [[ -n "${homd_rank:-}" ]]; then
    echo "[INFO] Representative selection ENABLED at rank: ${homd_rank}" >&2

    # ── Step 2.5: Harvest assembly stats (contigs, size) from PROKKA summary/ ──
    # Each summary/<GCA>.txt is ~250 bytes and carries 'contigs:' and 'bases:'.
    # Counts are byte-identical to deriving them from fna/<GCA>.fna.
    STATS_TSV="${HOMD_DIR}/HOMD_assembly_stats.tsv"
    STATS_STAGE="$(mktemp -d -p "$TMPDIR" homd_stats.XXXXXX)"
    export SUMMARY_BASE_URL="$homd_prokka_summary_url"
    export STATS_STAGE

    stat_one() {
        local gca="$1"
        local txt
        txt="$(curl --http1.1 -gsSL --retry 5 --retry-all-errors --retry-delay 2 \
                    --connect-timeout 30 --max-time 120 \
                    "${SUMMARY_BASE_URL}/${gca}.txt")" || return 0
        # Parse 'contigs:' and 'bases:'; emit only when both are present.
        awk -v gca="$gca" '
            /^contigs:/ { c = $2 }
            /^bases:/   { b = $2 }
            END { if (c != "" && b != "") printf "%s\t%s\t%s\n", gca, c, b }
          ' <<< "$txt" > "${STATS_STAGE}/${gca}.stat"
    }
    export -f stat_one

    echo "[INFO] Harvesting assembly stats for ${PLAN_COUNT} candidates (${homd_prokka_jobs} job(s))..." >&2
    xargs -P "$homd_prokka_jobs" -n 2 bash -c 'stat_one "$1"' _ < "$SELECTED_PLAN"

    { printf 'gca\tcontigs\tsize\n'; cat "${STATS_STAGE}"/*.stat 2>/dev/null; } > "$STATS_TSV"
    rm -rf "$STATS_STAGE"
    STATS_COUNT="$(($(wc -l < "$STATS_TSV") - 1))"
    echo "[INFO] Stats harvested: ${STATS_COUNT} / ${PLAN_COUNT} assemblies -> ${STATS_TSV}" >&2
    [[ "$STATS_COUNT" -gt 0 ]] || { echo "ERROR: no assembly stats could be harvested." >&2; exit 1; }

    # ── FALLBACK (commented): derive stats from fna/<GCA>.fna instead of summary/.
    # Produces identical contig/base counts; use if summary/ ever disappears.
    #   export FNA_BASE_URL="${homd_prokka_release_url}/fna"; export STATS_STAGE
    #   stat_one() {
    #       local gca="$1"
    #       curl --http1.1 -gsSL --retry 5 --retry-all-errors --retry-delay 2 \
    #            --connect-timeout 30 --max-time 600 "${FNA_BASE_URL}/${gca}.fna" \
    #       | awk -v gca="$gca" '/^>/{c++;next}{b+=length($0)} END{if(c>0)printf "%s\t%d\t%d\n",gca,c,b}' \
    #       > "${STATS_STAGE}/${gca}.stat"
    #   }; export -f stat_one

    # ── Step 2.6: Pick one representative per taxonomic group ──────────────────
    homd_reps_report="${HOMD_DIR}/HOMD_representatives.tsv"
    REPS_PLAN="${TMPDIR}/HOMD_reps_plan.txt"
    python3 mod_IV_select_representatives.py \
        --tax-file "$TAX_TO_USE" \
        --gca-info "$gca_info_file" \
        --stats    "$STATS_TSV" \
        --rank     "$homd_rank" \
        --output   "$homd_reps_report" \
        --plan     "$REPS_PLAN"

    DOWNLOAD_PLAN="$REPS_PLAN"
    echo "[INFO] Representative report: ${homd_reps_report}" >&2
else
    echo "[INFO] Representative selection DISABLED; downloading all body-site matches." >&2
fi

DL_COUNT="$(wc -l < "$DOWNLOAD_PLAN" | tr -d ' ')"
echo "[INFO] Proteomes to download from HOMD PROKKA: ${DL_COUNT}" >&2
[[ "$DL_COUNT" -gt 0 ]] || { echo "ERROR: empty download plan." >&2; exit 1; }

# ── Step 3: Download <GCA>.faa (parallel), prepend GCA, stage per file ─────────

STAGE="$(mktemp -d -p "$TMPDIR" homd_prokka_stage.XXXXXX)"

# Prefix every .faa header with the GCA (so build_db.py resolves the assembly)
# and, when a PROKKA .tsv is present, append ' GN=<gene>' for loci that carry a
# gene symbol. The join is on the RAW PROKKA locus tag, which the .faa and .tsv
# share verbatim — so it is unaffected by the GCA prefix added here or by
# build_db.py's later prefix consolidation. FS='\t' on the .tsv keeps empty
# gene cells from collapsing and shifting columns. Falls back to a plain
# prefix when no usable .tsv is available, so gene simply stays NA.
prefix_and_annotate() {
    local gca="$1" faa="$2" tsv="$3" out="$4"
    if [[ -s "$tsv" ]] && head -n1 "$tsv" | grep -qi 'locus_tag'; then
        awk -v gca="$gca" '
            BEGIN { FS = "\t"; lc = 0; gc = 0 }
            FNR == NR {
                if (FNR == 1) {
                    for (i = 1; i <= NF; i++) {
                        h = tolower($i)
                        if (h == "locus_tag") lc = i
                        else if (h == "gene") gc = i
                    }
                    next
                }
                if (lc > 0 && gc > 0) {
                    g = $gc
                    if (g != "" && g != "-") gene[$lc] = g
                }
                next
            }
            /^>/ {
                line = substr($0, 2)
                sp = index(line, " ")
                locus = (sp > 0) ? substr(line, 1, sp - 1) : line
                g = gene[locus]
                if (g != "") print ">" gca "_" line " GN=" g
                else          print ">" gca "_" line
                next
            }
            { print }
        ' "$tsv" "$faa" > "$out"
    else
        sed "s/^>/>${gca}_/" "$faa" > "$out"
    fi
}
export -f prefix_and_annotate

dl_one() {
    # args: <GCA> <HMT(unused)>
    local gca="$1"
    local fname="${gca}.faa"
    local url="${FAA_BASE_URL}/${fname}"
    local raw="${STAGE}/${fname}.raw"
    local out="${STAGE}/${fname}.prefixed"
    local tsv="${STAGE}/${gca}.tsv.raw"

    for attempt in 1 2 3 4 5; do
        if curl --http1.1 -gsSL \
            --retry 5 --retry-all-errors --retry-delay 2 \
            --connect-timeout 30 --max-time 600 \
            "$url" -o "$raw" \
            && [[ -s "$raw" ]] \
            && head -n1 "$raw" | grep -q '^>'; then
            # Best-effort PROKKA .tsv fetch for gene symbols (never fatal).
            if [[ "${FETCH_GENE_NAMES}" == "true" ]]; then
                curl --http1.1 -gsSL \
                    --retry 3 --retry-all-errors --retry-delay 2 \
                    --connect-timeout 30 --max-time 300 \
                    "${TSV_BASE_URL}/${gca}.tsv" -o "$tsv" 2>/dev/null || true
            fi
            prefix_and_annotate "$gca" "$raw" "$tsv" "$out"
            rm -f "$raw" "$tsv"
            return 0
        fi
        sleep 3
    done
    echo "[WARN] failed to download ${fname} after retries" >&2
    rm -f "$raw" "$tsv"
    return 0   # don't abort the parallel run; shortfall detected below
}
export -f dl_one
export FAA_BASE_URL="$homd_prokka_faa_url"
export TSV_BASE_URL="$homd_prokka_tsv_url"
export FETCH_GENE_NAMES="$homd_fetch_gene_names"
export STAGE

echo "[INFO] Downloading ${DL_COUNT} proteomes with ${homd_prokka_jobs} parallel job(s)..." >&2
xargs -P "$homd_prokka_jobs" -n 2 bash -c 'dl_one "$1" "$2"' _ < "$DOWNLOAD_PLAN"

# ── Step 4: Concatenate and validate ──────────────────────────────────────────

COMBINED_FAA="$(mktemp -p "$TMPDIR" homd_combined.XXXXXX.faa)"
: > "$COMBINED_FAA"

DOWNLOADED=0
while IFS= read -r -d '' f; do
    cat "$f" >> "$COMBINED_FAA"
    DOWNLOADED=$(( DOWNLOADED + 1 ))
done < <(find "$STAGE" -type f -name '*.prefixed' -print0)

rm -rf "$STAGE"

echo "[INFO] Proteomes successfully downloaded: ${DOWNLOADED} / ${DL_COUNT}" >&2
if [[ "$DOWNLOADED" -lt "$DL_COUNT" ]]; then
    echo "[WARN] $(( DL_COUNT - DOWNLOADED )) proteome(s) failed to download; see WARN lines above." >&2
fi

[[ -s "$COMBINED_FAA" ]] || {
    echo "ERROR: No protein sequences were downloaded. Check network and GCA selection." >&2
    rm -f "$COMBINED_FAA"; exit 1
}

FIRST_LINE="$(grep -m1 '.' "$COMBINED_FAA" || true)"
if [[ "${FIRST_LINE:0:1}" != ">" ]]; then
    echo "ERROR: Combined output does not appear to be a FASTA." >&2
    echo "       First line: ${FIRST_LINE:0:120}" >&2
    rm -f "$COMBINED_FAA"; exit 1
fi

mv -f "$COMBINED_FAA" "$homd_fasta_file"

SEQ_COUNT="$(grep -c '^>' "$homd_fasta_file")"
echo "[INFO] Done. Total HOMD (PROKKA) proteins written: ${SEQ_COUNT}" >&2
echo "[INFO] Output: ${homd_fasta_file}" >&2
echo "[INFO] Selection provenance: ${PLAN_REPORT}" >&2
if [[ -n "${homd_rank:-}" ]]; then
    echo "[INFO] Representative provenance: ${homd_reps_report}" >&2
fi
