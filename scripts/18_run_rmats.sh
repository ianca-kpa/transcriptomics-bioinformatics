
#!/usr/bin/env bash
set -euo pipefail

# Module II - Event-level alternative splicing with rMATS
# Adapted from Johns Hopkins University course

# Load local configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
source "${PROJECT_DIR}/configs/paths.sh"

# Directories
DATADIR="$STAR_BAM_DIR"
WORKDIR="${RMATS_DIR}"
ANNOT="$REFERENCE_GTF"

OUTDIR="${WORKDIR}/Ctrl-vs-Test"
TMPDIR="${WORKDIR}/tmp_Ctrl-vs-Test"

# Group 1: Month 1
CTRL_IDS=(
    SRR3734796
    SRR3734797
    SRR3734798
)

# Group 2: Month 4
TEST_IDS=(
    SRR3734820
    SRR3734824
    SRR3734825
)

# Create BAM lists
# mkdir -p "$WORKDIR"
mkdir -p "$OUTDIR" "$TMPDIR"

CTRLBAML="${WORKDIR}/ctrl_bams.rmats.list"
TESTBAML="${WORKDIR}/test_bams.rmats.list"

make_bam_list() {
    local -a paths=()
    local id bam

    for id in "$@"; do
        bam="${DATADIR}/${id}_Aligned.sortedByCoord.out.bam"

        if [[ ! -s "$bam" ]]; then
            echo "ERROR: Missing or empty BAM: $bam" >&2
            return 1
        fi

        paths+=("$bam")
    done

    local IFS=,
    printf '%s\n' "${paths[*]}"
}

make_bam_list "${CTRL_IDS[@]}" > "$CTRLBAML"
make_bam_list "${TEST_IDS[@]}" > "$TESTBAML"

# Check annotation and software
test -s "$ANNOT" || {
    echo "ERROR: Reference GTF not found: $ANNOT" >&2
    exit 1
}

command -v rmats.py >/dev/null || {
    echo "ERROR: rmats.py not found" >&2
    exit 1
}

# Run rMATS
rmats.py \
    --b1 "$CTRLBAML" \
    --b2 "$TESTBAML" \
    --gtf "$ANNOT" \
    --od "$OUTDIR" \
    --tmp "$TMPDIR" \
    -t single \
    --readLength 100 \
    --novelSS \
    --allow-clipping \
    --nthread 2 \
    > "${WORKDIR}/rmats.log" 2>&1

echo "rMATS analysis completed!"
