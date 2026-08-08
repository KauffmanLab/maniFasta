#!/bin/bash

## Originated by Christopher Handelmann

# Define the config file here so the rest of the scripts can dynamically source it
HERE="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"   # = …/01.scripts
cd "$HERE"                                                  # submit from 01.scripts
Config_file="$(dirname "$HERE")/00.setup/000.maniFasta.config"

export Config_file
# Source the config to set variables in this shell
source "${Config_file}"

# Set run variables
runScript=runBuild_db.sh

# Slurm log goes into the project directory alongside the run output dirs.
# %j is replaced by the job ID so logs from different runs don't collide.
LOG_DIR="${mainDIR}"

# Submit run script to slurm scheduler [choose cluster in config file]

SBATCH_ARGS=(
  --job-name="$build_JobName"
  --nodes="$build_Nodes"
  --cpus-per-task="$build_cpuspertask"
  --mem="$build_mem"
  --time="$build_Time"
  --output="${LOG_DIR}/slurm-%j.out"
  --error="${LOG_DIR}/slurm-%j.out"
)

# Cluster, partition, and QOS may be left blank when not required.
if [[ -n "${cluster:-}" ]]; then
  SBATCH_ARGS+=(--cluster="$cluster")
fi

if [[ -n "${partition:-}" ]]; then
  SBATCH_ARGS+=(--partition="$partition")
fi

if [[ -n "${qos:-}" ]]; then
  SBATCH_ARGS+=(--qos="$qos")
fi

sbatch "${SBATCH_ARGS[@]}" "${runScript}"