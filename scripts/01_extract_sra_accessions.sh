
#!/usr/bin/env bash
set -euo pipefail

# Get the absolute directory containing this script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Get the project root directory
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

# Module I sample metadata
METADATA_DIR="$PROJECT_DIR/course_materials/Module_I/metadata"
INPUT_FILE="$METADATA_DIR/SraRunTable.csv"
OUTPUT_FILE="$METADATA_DIR/run_id.txt"

# Check metadata input
if [[ ! -s "$INPUT_FILE" ]]; then
    echo "ERROR: SRA metadata not found or empty: $INPUT_FILE" >&2
    exit 1
fi

# Extract run IDs from the first column, skipping the header
tail -n +2 "$INPUT_FILE" | cut -d',' -f1 > "$OUTPUT_FILE" 

echo "Processing completed!"
echo "Run IDs saved to: $OUTPUT_FILE"
