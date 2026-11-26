#!/usr/bin/env bash
set -euo pipefail

# Get the absolute directory containing this script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Get the project root directory
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

source "${PROJECT_DIR}/configs/paths.sh"

SRA_DIR="${BASEDIR}/SRA"
FASTQ_DIR="${BASEDIR}/FASTQ"
# Override SRA_RUNS_FILE to convert a different module's run list.
RUNS_FILE="${SRA_RUNS_FILE:-${PROJECT_DIR}/course_materials/Module_II/metadata/selected_runs_m2.txt}"
if [[ "$RUNS_FILE" != /* ]]; then
    RUNS_FILE="${PROJECT_DIR}/${RUNS_FILE}"
fi

if [[ ! -s "$RUNS_FILE" ]]; then
    echo "ERROR: Run list not found or empty: $RUNS_FILE" >&2
    exit 1
fi

command -v fasterq-dump >/dev/null || {
    echo "ERROR: fasterq-dump not found. Activate the SRA Toolkit environment." >&2
    exit 1
}

mkdir -p "$FASTQ_DIR"

while IFS= read -r run || [[ -n "$run" ]]; do
    run="${run%$'\r'}"

    [[ -z "$run" || "$run" == \#* ]] && continue

    echo "Converting $run..."

    SRA_FILE="$SRA_DIR/$run/$run.sra"

    if [[ ! -s "$SRA_FILE" ]]; then
        echo "ERROR: SRA file not found or empty: $SRA_FILE" >&2
        exit 1
    fi

    fasterq-dump "$SRA_FILE" \
        --outdir "$FASTQ_DIR" \
        --outfile "$run.fastq" \
        --threads 4

    echo "Converted: $run"

done < "$RUNS_FILE"

echo "All conversions completed!"