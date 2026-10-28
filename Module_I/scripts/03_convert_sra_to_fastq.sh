#!/usr/bin/env bash
set -euo pipefail

# Get the absolute directory containing this script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Get the project root directory
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

# Define paths
SRA_DIR="$PROJECT_DIR/data/SRA"
FASTQ_DIR="$PROJECT_DIR/data/FASTQ"
RUN_LIST="$PROJECT_DIR/metadata/selected_runs.txt"

# Create output directory
mkdir -p "$FASTQ_DIR"

# Convert each SRA run to FASTQ
while IFS= read -r run || [[ -n "$run" ]]; do

    # Skip empty lines and comments
    [[ -z "$run" || "$run" == \#* ]] && continue

    echo "Converting $run..."

    # Locate the downloaded SRA file
    SRA_FILE="$SRA_DIR/$run/$run.sra"

    # Check if the SRA file exists
    if [[ ! -f "$SRA_FILE" ]]; then
        echo "ERROR: SRA file not found: $SRA_FILE" >&2
        exit 1
    fi

    # Convert SRA to FASTQ using fasterq-dump
    fasterq-dump "$SRA_FILE" \
        --outdir "$FASTQ_DIR" \
        --outfile "$run.fastq" \
        --threads 4

    echo "Converted: $run"

done < "$RUN_LIST"

echo "All conversions completed!"