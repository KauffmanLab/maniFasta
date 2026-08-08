#!/bin/bash

## Originated by Christopher Handelmann

set -euo pipefail

source ${Config_file}

# ── Validate config ───────────────────────────────────────────────────────────

case "${human_protein_set:-}" in
    canonical | canonical_isoform | canonical_trembl | canonical_isoform_trembl) ;;
    "")
        echo "ERROR: human_protein_set is not set in config." >&2
        echo "Valid options: canonical | canonical_isoform | canonical_trembl | canonical_isoform_trembl" >&2
        exit 1
        ;;
    *)
        echo "ERROR: Unknown human_protein_set: '${human_protein_set}'" >&2
        echo "Valid options: canonical | canonical_isoform | canonical_trembl | canonical_isoform_trembl" >&2
        exit 1
        ;;
esac

[[ -n "${human_fasta_file:-}" ]] || {
    echo "ERROR: human_fasta_file is not set in config." >&2; exit 1
}

# ── URL selection ─────────────────────────────────────────────────────────────
#
#   canonical               SwissProt only, no isoforms
#                           ~20,000 sequences
#
#   canonical_isoform       SwissProt + manually reviewed isoforms
#                           ~42,000 sequences
#
#   canonical_trembl        SwissProt + TrEMBL (unreviewed), no isoforms
#                           ~100,000+ sequences
#
#   canonical_isoform_trembl  SwissProt + isoforms + TrEMBL (unreviewed)
#                           ~250,000+ sequences; downloaded compressed

COMPRESSED=false

case "${human_protein_set}" in

    canonical)
        URL="https://rest.uniprot.org/uniprotkb/stream?\
format=fasta&\
query=%28%28proteome%3AUP000005640%29+AND+reviewed%3Dtrue%29"
        ;;

    canonical_isoform)
        URL="https://rest.uniprot.org/uniprotkb/stream?\
format=fasta&\
includeIsoform=true&\
query=%28%28proteome%3AUP000005640%29+AND+reviewed%3Dtrue%29"
        ;;

    canonical_trembl)
        # Swiss-Prot + TrEMBL, no isoforms: drop the reviewed filter (keeps
        # both reviewed and unreviewed) and omit includeIsoform.
        URL="https://rest.uniprot.org/uniprotkb/stream?\
format=fasta&\
query=%28%28proteome%3AUP000005640%29%29"
        ;;

    canonical_isoform_trembl)
        URL="https://rest.uniprot.org/uniprotkb/stream?\
compressed=true&\
format=fasta&\
includeIsoform=true&\
query=%28%28proteome%3AUP000005640%29%29"
        COMPRESSED=true
        ;;
esac

echo "==> Human protein set : ${human_protein_set}" >&2
echo "==> Output file       : ${human_fasta_file}"  >&2
echo "==> Compressed        : ${COMPRESSED}"        >&2

# ── Prepare output directory ──────────────────────────────────────────────────

mkdir -p "$(dirname "$human_fasta_file")"

# ── Download ──────────────────────────────────────────────────────────────────
#
# Retry logic mirrors fetch_to_file() in Fetch_Species_INFO_2.0.sh:
#   --retry 5 --retry-all-errors --retry-delay 2

if [[ "$COMPRESSED" == "true" ]]; then

    TMP_GZ="${human_fasta_file}.gz"

    echo "==> Downloading compressed FASTA to ${TMP_GZ} ..." >&2

    curl --http1.1 -gSL \
        --retry 5 --retry-all-errors --retry-delay 2 \
        --connect-timeout 30 --max-time 3600 \
        "$URL" -o "$TMP_GZ"

    echo "==> Decompressing ..." >&2
    gunzip -c "$TMP_GZ" > "$human_fasta_file"
    rm -f "$TMP_GZ"

else

    echo "==> Downloading FASTA to ${human_fasta_file} ..." >&2

    curl --http1.1 -gSL \
        --retry 5 --retry-all-errors --retry-delay 2 \
        --connect-timeout 30 --max-time 3600 \
        "$URL" -o "$human_fasta_file"

fi

# ── Validate output ───────────────────────────────────────────────────────────

[[ -s "$human_fasta_file" ]] || {
    echo "ERROR: Output file is empty or missing: ${human_fasta_file}" >&2; exit 1
}

# Confirm the file looks like a FASTA (first non-empty line starts with >)
FIRST_LINE="$(grep -m1 '.' "$human_fasta_file" || true)"
if [[ "${FIRST_LINE:0:1}" != ">" ]]; then
    echo "ERROR: Downloaded file does not appear to be a FASTA." >&2
    echo "       First line: ${FIRST_LINE:0:120}" >&2
    exit 1
fi

# Count sequences
SEQ_COUNT="$(grep -c '^>' "$human_fasta_file")"
echo "==> Done. Sequences written: ${SEQ_COUNT}" >&2
echo "==> Output: ${human_fasta_file}" >&2

# ── Record downloaded accessions ──────────────────────────────────────────────

if [[ -n "${RUN_DIR:-}" ]]; then
    HUMAN_ACC_DIR="${RUN_DIR}/run_logs"
else
    # Fallback for standalone invocation outside the pipeline: sit next to the FASTA.
    HUMAN_ACC_DIR="$(dirname "$human_fasta_file")"
fi
mkdir -p "$HUMAN_ACC_DIR"
HUMAN_ACC_FILE="${HUMAN_ACC_DIR}/human_proteins.${human_protein_set}.accessions.tsv"

{
    printf 'accession\tfetch_source\treview_status\tsource_label\tnote\n'
    grep '^>' "$human_fasta_file" \
        | awk -F'|' -v set="$human_protein_set" 'BEGIN{OFS="\t"}
            {
                db  = $1; sub(/^>/, "", db)
                acc = $2
                if (acc == "") next
                status = (db == "sp" ? "reviewed" : (db == "tr" ? "unreviewed" : db))
                print acc, "uniprot", status, "HUMAN", "human_uniprot " set " UP000005640"
            }'
} > "$HUMAN_ACC_FILE"

HUMAN_ACC_COUNT="$(($(grep -c . "$HUMAN_ACC_FILE") - 1))"
echo "==> Accession record: ${HUMAN_ACC_FILE} (${HUMAN_ACC_COUNT} accessions)" >&2
