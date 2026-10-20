
#!/usr/bin/env bash
set -euo pipefail

# Get the absolute directory containing this script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Get the project root directory
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

# Define paths
METADATA_DIR="$PROJECT_DIR/metadata"
INPUT_FILE="$METADATA_DIR/SraRunTable.csv"
OUTPUT_FILE="$METADATA_DIR/run_id.txt"

# Extract run IDs from the first column, skipping the header
tail -n +2 "$INPUT_FILE" | cut -d',' -f1 > "$OUTPUT_FILE"

echo "Processing completed!"
echo "Run IDs saved to: $OUTPUT_FILE"
