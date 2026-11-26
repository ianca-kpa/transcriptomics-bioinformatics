
#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Module II — Summarize significant rMATS skipped exon events
#
# Input:  SE.MATS.JC.txt (junction counts, skipped exon events)
# Output: SE.MATS.JC.significant_FDR0.05_dPSI0.05.tsv
#
# The output keeps all original rMATS columns and is sorted by:
#   1. FDR (ascending)
#   2. Absolute IncLevelDifference (descending)
#
# Group 1: 1-month samples
# Group 2: 4-month samples
# IncLevelDifference = mean inclusion level 1 - level 2
# ============================================================

# 1. Load machine-specific paths from the repository root.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
source "${PROJECT_DIR}/configs/paths.sh"

OUTDIR="${RMATS_DIR}/Ctrl-vs-Test"
INPUT="${OUTDIR}/SE.MATS.JC.txt"
OUTPUT="${OUTDIR}/SE.MATS.JC.significant_FDR0.05_dPSI0.05.tsv"

# 2. Significance thresholds.
FDR_THRESHOLD=0.05
DELTA_PSI_THRESHOLD=0.05

# 3. Validate input and prevent accidental replacement
# of the raw file.
if [[ ! -s "$INPUT" ]]; then
    echo "ERROR: rMATS output missing or empty: $INPUT" >&2
    exit 1
fi

if [[ "$INPUT" == "$OUTPUT" ]]; then
    echo "ERROR: Input and output paths must differ." >&2
    exit 1
fi

# Temporary files allow us to publish an output only after successful processing.
# The records file has two temporary sorting keys before each original rMATS row.
TMP_RECORDS="$(mktemp "${OUTDIR}/.se-records.XXXXXX")"
TMP_OUTPUT="$(mktemp "${OUTDIR}/.se-output.XXXXXX")"

# Remove temporary files when the script exits.
trap 'rm -f "$TMP_RECORDS" "$TMP_OUTPUT"' EXIT

printf '========================================\n'
printf 'rMATS — Skipped exon analysis\n'
printf '========================================\n'
printf 'Input: %s\n' "$INPUT"
printf 'FDR threshold: %s\n' "$FDR_THRESHOLD"
printf 'Minimum |delta PSI|: %s\n\n' "$DELTA_PSI_THRESHOLD"


# ------------------------------------------------------------
# 4. Read the rMATS output and filter significant events
# ------------------------------------------------------------

# Locate columns by their names rather than fixed positions.
# Preserve all original columns in the exported TSV.

awk -F '\t' -v OFS='\t' \
    -v fdr_threshold="$FDR_THRESHOLD" \
    -v delta_threshold="$DELTA_PSI_THRESHOLD" \
    -v records_file="$TMP_RECORDS" \
    -v output_file="$TMP_OUTPUT" '

# Check whether a field contains a valid numeric value.
function numeric(x) {
    return x ~ /^[-+]?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$/
}

{
    # Handle Windows line endings.
    sub(/\r$/, "", $0)
}


# ------------------------------------------------------------
# 4.1. Identify required columns
# ------------------------------------------------------------

NR == 1 {
    for (i = 1; i <= NF; i++) {
        cols[$i] = i
    }

    if (!("FDR" in cols) ||
        !("IncLevelDifference" in cols)) {
        print "ERROR: Required columns FDR / IncLevelDifference not found." > "/dev/stderr"
        failed = 1
        exit 1
    }

    fdr_col = cols["FDR"]
    delta_col = cols["IncLevelDifference"]
    required_col = (fdr_col > delta_col) ?
                   fdr_col : delta_col
    # Preserve the original header.
    print $0 > output_file
    close(output_file)
    printf "Detected columns: FDR=%d, IncLevelDifference=%d\n",
           fdr_col, delta_col
    next
}

# Skip empty lines.
NF == 0 || (NF == 1 && $1 == "") {
    next
}


# ------------------------------------------------------------
# 4.2. Process skipped exon events
# ------------------------------------------------------------

{
    # Check for incomplete rows.
    if (NF < required_col) {
        printf "ERROR: Missing columns at line %d.\n", NR \
            > "/dev/stderr"
        failed = 1
        exit 1
    }

    total++

    fdr = $fdr_col
    delta = $delta_col

    # Display the first event as a basic reference.
    if (total == 1) {
        gene = ("geneSymbol" in cols) ?
               $(cols["geneSymbol"]) : "NA"
        printf "First event: gene=%s, FDR=%s, deltaPSI=%s\n",
               gene, fdr, delta
    }

    # Exclude events with missing or invalid statistics.
    if (!numeric(fdr) || !numeric(delta)) {
        invalid_stats++
        next
    }

    tested++

    # Convert delta PSI to a number.
    dpsi = delta + 0
    # Calculate absolute delta PSI.
    abs_dpsi = (dpsi < 0) ? -dpsi : dpsi

    # --------------------------------------------------------
    # 4.3. Apply significance thresholds
    # --------------------------------------------------------

    if ((fdr + 0) < fdr_threshold &&
        abs_dpsi >= delta_threshold) {
        significant++

        # Positive delta PSI:
        # Higher exon inclusion in group 1 (1 month).
        if (dpsi > 0) {
            higher_group1++
        }

        # Negative delta PSI:
        # Higher exon inclusion in group 2 (4 months).
        else if (dpsi < 0) {
            higher_group2++
        }

        # Store:
        #   Column 1: FDR (ascending sort key)
        #   Column 2: |delta PSI| (descending sort key)
        #   Remaining columns: original rMATS record
        # The temporary columns will be removed after sorting.

        printf "%s\t%.12g\t%s\n",
               fdr, abs_dpsi, $0 > records_file
    }
}


# ------------------------------------------------------------
# 4.4. Print statistical summary
# ------------------------------------------------------------

END {
    close(records_file)
    if (failed) {
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
           invalid_stats
    printf "Significant SE events (FDR < %.2f, |deltaPSI| >= %.2f): %d\n",
           fdr_threshold, delta_threshold, significant
    printf "Higher exon inclusion in group 1 (1mo): %d\n",
           higher_group1
    printf "Higher exon inclusion in group 2 (4mo): %d\n",
           higher_group2
}

' "$INPUT"


# ------------------------------------------------------------
# 5. Sort and export significant events
# ------------------------------------------------------------

# Sort by:
#   1. FDR: smallest to largest
#   2. Absolute delta PSI: largest to smallest
#
# Remove the two temporary sorting columns.
# Keep the original rMATS columns unchanged.
# If no events pass the filters, the exported TSV
# will contain only the header.

LC_ALL=C sort -t $'\t' -k1,1g -k2,2gr "$TMP_RECORDS" \
    | cut -f3- >> "$TMP_OUTPUT"

# Publish the completed output file.
mv -f "$TMP_OUTPUT" "$OUTPUT"


# ------------------------------------------------------------
# 6. Final report
# ------------------------------------------------------------

printf '\nSaved significant events (original rMATS columns):\n'
printf '%s\n' "$OUTPUT"
printf 'Exported event rows: %s\n' \
    "$(($(wc -l < "$OUTPUT") - 1))"
printf '\nrMATS summary completed.\n'
