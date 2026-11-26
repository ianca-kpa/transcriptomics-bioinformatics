
#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Module II — Summarize rMATS skipped exon (SE) results
#
# Input:
#   SE.MATS.JC.txt
#
# Analysis:
#   - Count detected skipped exon events
#   - Identify events with valid statistical results
#   - Count significant events using FDR and delta PSI
#   - Summarize the direction of exon inclusion changes
#
# Comparison:
#   Group 1: 1-month samples
#   Group 2: 4-month samples
#
# IncLevelDifference = IncLevel1 - IncLevel2
# ============================================================


# ------------------------------------------------------------
# 1. Project directories
# ------------------------------------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

source "${PROJECT_DIR}/configs/paths.sh"

# rMATS junction-count output for skipped exons.
FILE="${RMATS_DIR}/Ctrl-vs-Test/SE.MATS.JC.txt"


# ------------------------------------------------------------
# 2. Statistical thresholds
# ------------------------------------------------------------

# False discovery rate threshold.
FDR_THRESHOLD=0.05

# Minimum absolute difference in exon inclusion.
# 0.05 corresponds to 5 percentage points.
DELTA_PSI_THRESHOLD=0.05


# ------------------------------------------------------------
# 3. Validate input file
# ------------------------------------------------------------

if [[ ! -s "$FILE" ]]; then
    echo "ERROR: rMATS output file not found or empty:" >&2
    echo "$FILE" >&2
    exit 1
fi

echo "========================================"
echo "rMATS — Skipped exon analysis"
echo "========================================"
echo "Input file: $FILE"
echo "FDR threshold: $FDR_THRESHOLD"
echo "Minimum |delta PSI|: $DELTA_PSI_THRESHOLD"
echo


# ------------------------------------------------------------
# 4. Summarize skipped exon events
# ------------------------------------------------------------

awk -F '\t' \
    -v fdr_threshold="$FDR_THRESHOLD" \
    -v delta_threshold="$DELTA_PSI_THRESHOLD" '

# Check whether a field contains a valid numeric value.
# This prevents NA, empty strings and other invalid values
# from being treated as numbers.
function valid_number(x) {
    return x ~ /^[-+]?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$/
}


# ------------------------------------------------------------
# 4.1. Identify columns from the header
# ------------------------------------------------------------

NR == 1 {
    for (i = 1; i <= NF; i++) {
        # Remove Windows carriage returns if present.
        sub(/\r$/, "", $i)
        cols[$i] = i
    }

    # Required statistical columns.
    if (!("FDR" in cols) ||
        !("IncLevelDifference" in cols)) {
        print "ERROR: Required rMATS columns not found." \
            > "/dev/stderr"
        invalid_header = 1
        exit 1
    }

    fdr_col = cols["FDR"]
    delta_col = cols["IncLevelDifference"]
    required_col = (fdr_col > delta_col) ?
                   fdr_col : delta_col
    printf "Detected columns: FDR=%d, IncLevelDifference=%d\n",
           fdr_col, delta_col
    next
}

# Ignore empty lines.
NF == 0 {
    next
}


# ------------------------------------------------------------
# 4.2. Process individual SE events
# ------------------------------------------------------------

{
    # Reject incomplete rows.
    if (NF < required_col) {
        printf "ERROR: Incomplete row at line %d.\n", NR \
            > "/dev/stderr"
        invalid_row = 1
        exit 1
    }

    total++

    # Display one event as an example.
    if (total == 1) {
        gene = ("geneSymbol" in cols) ?
               $(cols["geneSymbol"]) : "NA"
        printf "First event: gene=%s, FDR=%s, deltaPSI=%s\n",
               gene, $fdr_col, $delta_col
    }

    fdr = $fdr_col
    delta = $delta_col
    sub(/\r$/, "", delta)

    # Skip events without valid statistical values.
    if (!valid_number(fdr) || !valid_number(delta)) {
        missing_stats++
        next
    }
    tested++

    # Convert the difference to a numeric value.
    effect = delta + 0

    # Calculate absolute delta PSI.
    abs_effect = (effect < 0) ? -effect : effect

    # --------------------------------------------------------
    # 4.3. Apply significance thresholds
    # --------------------------------------------------------

    if (fdr + 0 < fdr_threshold &&
        abs_effect >= delta_threshold) {
        significant++
        # Positive delta PSI:
        # higher exon inclusion in group 1 (1 month).
        if (effect > 0) {
            higher_group1++
        }
        # Negative delta PSI:
        # higher exon inclusion in group 2 (4 months).
        else if (effect < 0) {
            higher_group2++
        }
    }
}

# ------------------------------------------------------------
# 4.4. Print summary
# ------------------------------------------------------------

END {
    if (invalid_header || invalid_row) {
        exit 1
    }
    print ""
    print "========================================"
    print "Summary"
    print "========================================"
    printf "Total SE events: %d\n", total
    printf "Events with numeric FDR and deltaPSI: %d\n",
           tested
    printf "Events with missing/invalid statistics: %d\n",
           missing_stats
    print ""
    printf "Significant SE events (FDR < %.2f, |deltaPSI| >= %.2f): %d\n",
           fdr_threshold, delta_threshold, significant
    printf "Higher exon inclusion in group 1 (1mo): %d\n",
           higher_group1
    printf "Higher exon inclusion in group 2 (4mo): %d\n",
           higher_group2
}
' "$FILE"

echo
echo "rMATS summary completed."