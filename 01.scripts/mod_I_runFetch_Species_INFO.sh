#!/bin/bash

source ${Config_file}

mkdir -p "$(dirname "$Species_resolved")"

# Set sleep rate based on whether an API key is present.
# NCBI allows 10 req/s with a key and 3 req/s without. Both are overridable via
# ncbi_sleep_api_key / ncbi_sleep_no_key in 00.setup/000.maniFasta.config.
if [[ -n "${ncbi_api_key:-}" ]]; then
    SLEEP_SEC="${ncbi_sleep_api_key:-0.11}"
else
    SLEEP_SEC="${ncbi_sleep_no_key:-0.4}"
fi

EMAIL_FLAG=()
KEY_FLAG=()
[[ -n "${ncbi_email:-}"   ]] && EMAIL_FLAG=( --email      "$ncbi_email"   )
[[ -n "${ncbi_api_key:-}" ]] && KEY_FLAG=(   --api_key    "$ncbi_api_key" )

bash mod_I_Fetch_UserSpecifiedSpecies_INFO.sh \
  --input  "$Species_fetch_list" \
  --output "$Species_resolved"   \
  --sleep  "$SLEEP_SEC"          \
  "${EMAIL_FLAG[@]}"             \
  "${KEY_FLAG[@]}"
