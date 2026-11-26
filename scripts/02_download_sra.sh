
#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

source "${PROJECT_DIR}/configs/paths.sh"

SRA_DIR="${BASEDIR}/SRA"
# Override SRA_RUNS_FILE to download a different module's run list.
RUNS_FILE="${SRA_RUNS_FILE:-${PROJECT_DIR}/course_materials/Module_I/metadata/selected_runs.txt}"
if [[ "$RUNS_FILE" != /* ]]; then
    RUNS_FILE="${PROJECT_DIR}/${RUNS_FILE}"
fi

if [[ ! -s "$RUNS_FILE" ]]; then
    echo "ERROR: Run list not found or empty: $RUNS_FILE" >&2
    exit 1
fi

command -v prefetch >/dev/null || {
    echo "ERROR: prefetch not found. Activate the SRA Toolkit environment." >&2
    exit 1
}

mkdir -p "$SRA_DIR"

echo "Run list: $RUNS_FILE"
echo "SRA output: $SRA_DIR"

while IFS= read -r run || [[ -n "$run" ]]; do
    run="${run%$'\r'}"

    [[ -z "$run" || "$run" == \#* ]] && continue

    echo "Downloading $run..."

    prefetch "$run" -O "$SRA_DIR"

    echo "Downloaded: $run"

done < "$RUNS_FILE"

echo "All downloads completed!"
