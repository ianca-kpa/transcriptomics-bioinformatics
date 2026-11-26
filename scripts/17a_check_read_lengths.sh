
#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Module II — Validate RNA-seq reads before running rMATS
#
# Checks:
#   1. Read length distribution
#   2. Single-end or paired-end alignment flags
#   3. Basic BAM integrity
#
# Available modes:
#   one_sample  : first N alignments from one sample
#   six_samples : first N alignments from each sample
#   full_file   : all alignments from one sample
#
# Activate only ONE mode.
# ============================================================


# ------------------------------------------------------------
# 1. Analysis settings
# ------------------------------------------------------------

MODE="one_sample"
# MODE="one_sample"
# MODE="full_file"

# Sample used in one_sample and full_file modes
SAMPLE_ID="SRR3734796"

# Number of alignments inspected in limited modes
MAX_READS=10000


# ------------------------------------------------------------
# 2. Load project configuration
# ------------------------------------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

source "${PROJECT_DIR}/configs/paths.sh"

# Use the system SAMtools.
# This avoids older versions installed in some Conda environments.
SAMTOOLS="/usr/bin/samtools"


# ------------------------------------------------------------
# 3. Check dependencies and parameters
# ------------------------------------------------------------

if [[ ! -x "$SAMTOOLS" ]]; then
    echo "ERROR: SAMtools not found: $SAMTOOLS" >&2
    exit 1
fi

if ! [[ "$MAX_READS" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: MAX_READS must be a positive integer." >&2
    exit 1
fi


# ------------------------------------------------------------
# 4. Define experimental samples
# ------------------------------------------------------------

# Month 1: three biological replicates
# Month 4: three biological replicates

SAMPLE_IDS=(
    SRR3734796
    SRR3734797
    SRR3734798
    SRR3734820
    SRR3734824
    SRR3734825
)


# ------------------------------------------------------------
# 5. Function to validate BAM files
# ------------------------------------------------------------

check_bam() {
    local id="$1"
    local limit="$2"

    # limit = 0 means inspect the entire BAM.
    local bam="${STAR_BAM_DIR}/${id}_Aligned.sortedByCoord.out.bam"
    local -a statuses

    # Check that the BAM exists and is not empty.
    if [[ ! -s "$bam" ]]; then
        echo "ERROR: Missing or empty BAM: $bam" >&2
        return 1
    fi

    # Check BAM header and EOF marker.
    if ! "$SAMTOOLS" quickcheck -v "$bam"; then
        echo "ERROR: BAM failed quickcheck: $bam" >&2
        return 1
    fi

    echo
    echo "========================================"
    echo "Sample: $id"

    if (( limit == 0 )); then
        echo "Scope: entire BAM"
    else
        echo "Scope: first $limit primary mapped alignments"
    fi

    echo "========================================"

    # --------------------------------------------------------
    # 5.1. Extract primary mapped alignments
    # --------------------------------------------------------

    # samtools view -F 2308 excludes:
    #
    #   4    = unmapped reads
    #   256  = secondary alignments
    #   2048 = supplementary alignments
    #
    # SAM columns:
    #   $2  = SAM FLAG
    #   $10 = Read sequence
    #
    # SAM FLAG bit 0x1 indicates that a read originates
    # from a paired sequencing template.
    #
    # The FLAG is an integer. Its least significant bit
    # can be checked using flag % 2.
    #
    # When AWK stops after MAX_READS, samtools may receive
    # SIGPIPE (exit status 141). This is expected when
    # analyzing only part of a BAM file.

    set +e

    "$SAMTOOLS" view -F 2308 "$bam" |
    awk -v limit="$limit" '

        $10 != "*" {

            # -----------------------------------------------
            # Read length
            # -----------------------------------------------

            len = length($10)
            lengths[len]++
            total++
            if (total == 1 || len < min_len)
                min_len = len
            if (total == 1 || len > max_len)
                max_len = len


            # -----------------------------------------------
            # Sequencing layout from SAM flags
            # -----------------------------------------------

            flag = $2 + 0
            # FLAG 0x1: paired template
            if (flag % 2 == 1) {
                paired++
                # FLAG 0x40: first segment in template
                if (int(flag / 64) % 2 == 1)
                    read1++
                # FLAG 0x80: last segment in template
                if (int(flag / 128) % 2 == 1)
                    read2++
            } else {
                unpaired++
            }
            # Stop after N alignments when requested.
            if (limit > 0 && total >= limit)
                exit
        }
        END {
            if (total == 0) {
                print "ERROR: No usable alignments found." \
                      > "/dev/stderr"
                exit 1
            }


            # -----------------------------------------------
            # Report read length distribution
            # -----------------------------------------------

            print ""
            print "Read length distribution:"

            for (len = min_len; len <= max_len; len++) {
                if (len in lengths) {
                    printf "  %d nt: %d reads (%.1f%%)\n",
                        len,
                        lengths[len],
                        100 * lengths[len] / total
                }
            }


            # -----------------------------------------------
            # Report paired-end flags
            # -----------------------------------------------

            print ""
            print "Sequencing layout (SAM flags):"
            printf "  Total alignments: %d\n", total
            printf "  Flagged as paired: %d (%.1f%%)\n",
                paired,
                100 * paired / total
            printf "  Flagged as unpaired: %d (%.1f%%)\n",
                unpaired,
                100 * unpaired / total

            if (paired > 0) {
                printf "  First-in-pair flags: %d\n", read1
                printf "  Last-in-pair flags: %d\n", read2
            }


            # -----------------------------------------------
            # Interpret sequencing layout
            # -----------------------------------------------

            print ""
            print "Interpretation:"

            if (paired == 0) {
                print "  Consistent with SINGLE-END alignment."
            } else if (unpaired == 0) {
                print "  Consistent with PAIRED-END alignment."
            } else {
                print "  WARNING: Mixed paired/unpaired flags."
                print "  Investigate before running rMATS."
            }
            print ""
        }
    '

    # Capture exit codes from both commands in the pipeline.
    statuses=("${PIPESTATUS[@]}")

    set -e

    # Verify AWK execution.
    if (( statuses[1] != 0 )); then
        echo "ERROR: Failed to analyze $id." >&2
        return 1
    fi

    # Allow SIGPIPE when processing a limited subset.
    if (( statuses[0] != 0 )); then
        if (( limit > 0 && statuses[0] == 141 )); then
            : # Expected when AWK stops early.
        else
            echo "ERROR: samtools failed for $id." >&2
            return 1
        fi
    fi
}


# ------------------------------------------------------------
# 6. Run selected analysis mode
# ------------------------------------------------------------

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

echo "Read validation completed."
