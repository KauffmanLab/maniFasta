#!/bin/bash

## Downloads HOMD taxonomy table and GCA_ID_info.txt.
## GCA filename now encodes the HOMD genomic version to prevent
## collisions when running multiple builds with different versions.

set -euo pipefail

source "${Config_file}"

# ── Validate required config vars ─────────────────────────────────────────────

for var in homd_tax_raw_file gca_info_file homd_genomic_refseq_version; do
    [[ -n "${!var:-}" ]] || {
        echo "ERROR: '$var' is not set in config." >&2; exit 1
    }
done

mkdir -p "$(dirname "$homd_tax_raw_file")"

# gca_info_file is set in config (encodes version, lives in 02.reference_data/)
mkdir -p "$(dirname "$gca_info_file")"

# ── Download 1: HOMD taxonomy table ───────────────────────────────────────────

HOMD_TAX_URL="https://www.homd.org/download/dld_taxtable_all/browser"

# Optional pinned local table. Set homd_tax_local_file in config to skip the
# download and use a version-controlled copy instead. All validation below
# applies to both paths.
if [[ -n "${homd_tax_local_file:-}" ]]; then
    echo "==> Using PINNED local HOMD taxonomy table (download skipped)" >&2
    echo "==> Source : ${homd_tax_local_file}" >&2
    echo "==> Output : ${homd_tax_raw_file}"   >&2

    [[ -s "$homd_tax_local_file" ]] || {
        echo "ERROR: homd_tax_local_file is set but empty or missing: ${homd_tax_local_file}" >&2
        exit 1
    }

    cp -f "$homd_tax_local_file" "$homd_tax_raw_file"
else
    echo "==> Downloading HOMD taxonomy table..."  >&2
    echo "==> URL    : ${HOMD_TAX_URL}"            >&2
    echo "==> Output : ${homd_tax_raw_file}"       >&2

    curl --http1.1 -gSL \
        --retry 5 --retry-all-errors --retry-delay 2 \
        --connect-timeout 30 --max-time 600 \
        "$HOMD_TAX_URL" -o "$homd_tax_raw_file"
fi

[[ -s "$homd_tax_raw_file" ]] || {
    echo "ERROR: Taxonomy file is empty or missing: ${homd_tax_raw_file}" >&2; exit 1
}

FIRST_LINE="$(grep -m1 '.' "$homd_tax_raw_file" || true)"
if [[ "$FIRST_LINE" != "HOMD.org Taxon Data::"* ]]; then
    echo "ERROR: Downloaded file does not appear to be a valid HOMD taxonomy table." >&2
    echo "       Expected first line starting with 'HOMD.org Taxon Data::'" >&2
    echo "       Got: ${FIRST_LINE:0:120}" >&2
    exit 1
fi

ROW_COUNT="$(tail -n +2 "$homd_tax_raw_file" | grep -c '.' || true)"
echo "==> Done. Rows in table (excluding header): ${ROW_COUNT}" >&2
echo "==> Output: ${homd_tax_raw_file}" >&2

# Copy raw taxonomy to homd_tax_file (the name mod_IV_runFetch_HOMD_Proteomes.sh expects).
# mod_IV_homd_filter.py handles the metadata line in the raw format, so no transformation needed.
cp -f "$homd_tax_raw_file" "$homd_tax_file"
echo "==> Copied to: ${homd_tax_file}" >&2

# ── Download 2: GCA_ID_info.txt ───────────────────────────────────────────────
# Pull from the PROKKA release so GCA_ID_info.txt lists exactly the genomes that
# have PROKKA faa files (same set mod_IV_runFetch_HOMD_Proteomes.sh downloads).
#
# homd_prokka_base_url is the UNVERSIONED root (.../genomes/PROKKA). The version
# is resolved by the source planner from the HOMD row's genomic_refseq_version=
# option (or auto-detected) and arrives here as homd_genomic_refseq_version,
# which is validated as non-empty above.
: "${homd_prokka_base_url:=https://www.homd.org/ftp/genomes/PROKKA}"
homd_prokka_release_url="${homd_prokka_base_url%/}/${homd_genomic_refseq_version}"
GCA_INFO_URL="${homd_prokka_release_url}/GCA_ID_info.txt"

echo "" >&2
echo "==> Downloading HOMD GCA_ID_info.txt..."            >&2
echo "==> Version: ${homd_genomic_refseq_version}"        >&2
echo "==> URL    : ${GCA_INFO_URL}"                       >&2
echo "==> Output : ${gca_info_file}"                      >&2

curl --http1.1 -gSL \
    --retry 5 --retry-all-errors --retry-delay 2 \
    --connect-timeout 30 --max-time 600 \
    "$GCA_INFO_URL" -o "$gca_info_file"

[[ -s "$gca_info_file" ]] || {
    echo "ERROR: GCA_ID_info.txt is empty or missing: ${gca_info_file}" >&2; exit 1
}

FIRST_DATA_LINE="$(grep -m1 '^GCA_' "$gca_info_file" || true)"
if [[ -z "$FIRST_DATA_LINE" ]]; then
    echo "ERROR: Downloaded GCA_ID_info.txt does not appear valid." >&2
    echo "       Expected lines starting with 'GCA_'" >&2
    exit 1
fi

GCA_COUNT="$(grep -c '^GCA_' "$gca_info_file" || true)"
echo "==> Done. GCA entries: ${GCA_COUNT}" >&2
echo "==> Output: ${gca_info_file}" >&2
