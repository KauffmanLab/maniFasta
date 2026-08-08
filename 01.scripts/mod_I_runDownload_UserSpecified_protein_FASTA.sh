#!/usr/bin/env bash
set -euo pipefail

source "${Config_file}"
source "${setupDIR}/000.maniFasta.local.env"

# Make datasets binary findable by bare name for inner scripts
export PATH="$(dirname "$DATASETS_BIN"):$PATH"

mkdir -p "$TMPDIR"
mkdir -p "$proteomes_dir"

# Set sleep rate based on whether an API key is present
if [[ -n "${ncbi_api_key:-}" ]]; then
    SLEEP_SEC="0.11"
else
    SLEEP_SEC="0.4"
fi

EMAIL_FLAG=()
KEY_FLAG=()
[[ -n "${ncbi_email:-}"   ]] && EMAIL_FLAG=( --email   "$ncbi_email"   )
[[ -n "${ncbi_api_key:-}" ]] && KEY_FLAG=(   --api_key "$ncbi_api_key" )

# efetch is called via HTTP (eutils), so no conda env wrapper is needed here
bash ./mod_I_download_UserSpecified_protein_FASTA.sh \
  -i "$Species_resolved" \
  -o "$proteomes_dir"    \
  --sleep "$SLEEP_SEC"   \
  "${EMAIL_FLAG[@]}"     \
  "${KEY_FLAG[@]}"
