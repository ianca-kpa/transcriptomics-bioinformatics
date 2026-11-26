
#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Module II — Check read lengths before running rMATS
#
# Available modes:
#   one_sample  : first 10,000 alignments from one sample
#   six_samples : first 10,000 alignments from each sample
#   full_file   : all alignments from one sample
#
# Activate only ONE mode. Comment out the others using #.
# ============================================================

# 1. Choose analysis mode
MODE="one_sample"
# MODE="six_samples"
# MODE="full_file"

# Sample used in one_sample and full_file modes
SAMPLE_ID="SRR3734796"

# Number of alignments to inspect in limited modes
MAX_READS=10000


# 2. Load project paths
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

source "${PROJECT_DIR}/configs/paths.sh"


# 3. Check dependencies and parameters
if ! command -v samtools >/dev/null 2>&1; then
    echo "ERROR: samtools is not available." >&2
    exit 1
fi

if ! [[ "$MAX_READS" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: MAX_READS must be a positive integer." >&2
    exit 1
fi


# 4. Define experimental samples
SAMPLE_IDS=(
    SRR3734796
    SRR3734797
    SRR3734798
    SRR3734820
    SRR3734824
    SRR3734825
)


# 5. Function to check read lengths
check_bam() {
    local id="$1"
    local limit="$2"

    # limit = 0 means inspect the entire BAM
    local bam="${STAR_BAM_DIR}/${id}_Aligned.sortedByCoord.out.bam"
    local -a statuses

    # Check BAM existence
    if [[ ! -s "$bam" ]]; then
        echo "ERROR: Missing or empty BAM: $bam" >&2
        return 1
    fi

    # Check BAM header and EOF marker
    samtools quickcheck -v "$bam" || {
        echo "ERROR: BAM failed samtools quickcheck: $bam" >&2
        return 1
    }

    echo
    echo "Sample: $id"

    if (( limit == 0 )); then
        echo "Scope: entire BAM"
    else
        echo "Scope: first $limit alignments"
    fi

    # samtools -F 2308 excludes:
    #   4    = unmapped reads
    #   256  = secondary alignments
    #   2048 = supplementary alignments
    #
    # SAM column 10 ($10) contains the read sequence.
    # length($10) returns its length in nucleotides.
    #
    # When AWK stops after MAX_READS, samtools may receive
    # SIGPIPE (exit 141). This is expected in limited modes.

    set +e

    samtools view -F 2308 "$bam" | \
    awk -v limit="$limit" '
        $10 != "*" {
            len = length($10)
            counts[len]++
            total++
            if (limit > 0 && total >= limit) {
                exit
            }
        }

        END {
            if (total == 0) {
                print "ERROR: No alignments found." > "/dev/stderr"
                exit 1
            }
            for (len in counts) {
                printf "%d nt: %d reads (%.1f%%)\n",
                    len,
                    counts[len],
                    100 * counts[len] / total
            }
        }
    ' | sort -n

    statuses=("${PIPESTATUS[@]}")

    set -e

    # Check AWK and sort
    if (( statuses[1] != 0 || statuses[2] != 0 )); then
        echo "ERROR: Failed to summarize $id." >&2
        return 1
    fi

    # Check samtools, allowing expected SIGPIPE
    if (( statuses[0] != 0 )); then
        if (( limit > 0 && statuses[0] == 141 )); then
            : # Expected when AWK stops early
        else
            echo "ERROR: samtools failed for $id." >&2
            return 1
        fi
    fi
}


# 6. Execute selected analysis mode
# case "$VARIAVEL" in
#     valor1)
#         comando1
#         ;;
#     valor2)
#         comando2
#         ;;
#     *)
#         comando_padrao
#         ;;
# esac

case "$MODE" in
    one_sample)
        check_bam "$SAMPLE_ID" "$MAX_READS"
        ;;
    six_samples)
        for id in "${SAMPLE_IDS[@]}"; do
            check_bam "$id" "$MAX_READS"
        done
        ;;
    full_file)
        check_bam "$SAMPLE_ID" 0
        ;;
    *)
        echo "ERROR: Invalid MODE: $MODE" >&2
        exit 1
        ;;
esac
