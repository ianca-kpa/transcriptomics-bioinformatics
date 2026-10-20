
#!/usr/bin/env bash
set -euo pipefail

# Get the absolute directory containing this script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Get the project root directory
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

# Define paths
SRA_DIR="$PROJECT_DIR/data/SRA"
RUN_LIST="$PROJECT_DIR/metadata/selected_runs.txt"

# Create output directory
mkdir -p "$SRA_DIR"

# Download each SRA run
while IFS= read -r run || [[ -n "$run" ]]; do

    # Skip empty lines and comments
    [[ -z "$run" || "$run" == \#* ]] && continue

    echo "Downloading $run..."

    prefetch "$run" -O "$SRA_DIR"

    echo "Downloaded: $run"

done < "$RUN_LIST"

echo "All downloads completed!"
