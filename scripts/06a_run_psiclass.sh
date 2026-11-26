
#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Transcript assembly using PsiCLASS
#
# Based on the Johns Hopkins University transcriptomics course.
# Adapted for local execution on Ubuntu/WSL2.
#
# Supports:
#   Module I  - 6 RNA-seq samples
#   Module II - 9 RNA-seq samples
#
# Input:  STAR coordinate-sorted BAM files
# Output: PsiCLASS transcript assemblies (GTF)
# ============================================================


# ------------------------------------------------------------
# 1. Project directories and configuration
# ------------------------------------------------------------
# Locate this script and the repository root.
# This allows execution from any working directory.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

# Load machine-specific paths.
# The local configuration is excluded from Git.
CONFIG_FILE="${PROJECT_DIR}/configs/paths.sh"

if [[ ! -f "$CONFIG_FILE" ]]; then
    echo "ERROR: Configuration file not found: $CONFIG_FILE" >&2
    exit 1
fi

source "$CONFIG_FILE"

# STAR alignment directory.
# BAM files are expected directly inside STAR_BAM_DIR.
ALNDIR="$STAR_BAM_DIR"

# Sample list.
# Default: Module I (6 samples).
# Can be overridden using PSICLASS_RUNS_FILE.
RUNS_FILE="${PSICLASS_RUNS_FILE:-${PROJECT_DIR}/course_materials/Module_I/metadata/selected_runs.txt}"

# PsiCLASS output directory.
# Default: Module I output directory.
# Module II should use a different directory.
WORKDIR="${PSICLASS_OUTPUT_DIR:-$PSICLASS_DIR}"

# Expected number of biological samples.
# Default: 6 for Module I.
EXPECTED_SAMPLES="${PSICLASS_EXPECTED_SAMPLES:-6}"

# Validation-only mode.
# Set PSICLASS_CHECK_ONLY=1 to validate inputs
# without running transcript assembly.
CHECK_ONLY="${PSICLASS_CHECK_ONLY:-0}"


# ------------------------------------------------------------
# 2. Tools and analysis parameters
# ------------------------------------------------------------

# Number of threads used by PsiCLASS.
THREADS=2

# Use the system SAMtools (version 1.22.1).
# The older version in the psiclass Conda environment
# does not support the quickcheck command.
SAMTOOLS="/usr/bin/samtools"

# Check PsiCLASS availability.
if ! command -v psiclass >/dev/null 2>&1; then
    echo "ERROR: PsiCLASS not found." >&2
    exit 1
fi

# Check SAMtools availability.
if [[ ! -x "$SAMTOOLS" ]]; then
    echo "ERROR: SAMtools not found: $SAMTOOLS" >&2
    exit 1
fi

# Validate the expected sample count.
if ! [[ "$EXPECTED_SAMPLES" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: EXPECTED_SAMPLES must be a positive integer." >&2
    exit 1
fi


# ------------------------------------------------------------
# 3. Output files
# ------------------------------------------------------------
# List of BAM files provided to PsiCLASS.
BAMLIST="${WORKDIR}/bamlist.psiclass"

# Prefix used for PsiCLASS output files.
OUTPUT_PREFIX="${WORKDIR}/psiclass"

# File containing PsiCLASS messages and errors.
LOGFILE="${WORKDIR}/psiclass.log"

# Main expected consensus assembly.
CONSENSUS_GTF="${OUTPUT_PREFIX}_vote.gtf"


# ------------------------------------------------------------
# 4. Validate input directories and sample metadata
# ------------------------------------------------------------
# Verify that the BAM directory exists.
if [[ ! -d "$ALNDIR" ]]; then
    echo "ERROR: BAM directory not found: $ALNDIR" >&2
    exit 1
fi

# Verify that the selected run list exists.
if [[ ! -s "$RUNS_FILE" ]]; then
    echo "ERROR: Sample list missing or empty: $RUNS_FILE" >&2
    exit 1
fi

# Create the output directory if needed.
mkdir -p "$WORKDIR"

echo "========================================"
echo "PsiCLASS transcript assembly"
echo "========================================"
echo "BAM directory: $ALNDIR"
echo "Sample list: $RUNS_FILE"
echo "Output directory: $WORKDIR"
echo "Expected samples: $EXPECTED_SAMPLES"
echo


# ------------------------------------------------------------
# 5. Generate BAM list from selected SRA accessions
# ------------------------------------------------------------
# The sample list determines which BAM files are analyzed.
# This is safer than searching for every BAM in a directory,
# which could accidentally include additional samples.

echo "Preparing BAM list..."

# Write first to a temporary file.
# The final BAM list is created only after validation.
TEMP_BAMLIST="$(mktemp "${WORKDIR}/.bamlist.XXXXXX")"

# Remove the temporary file when the script exits.
trap 'rm -f "$TEMP_BAMLIST"' EXIT

# Track sample identifiers to detect duplicates.
declare -A SEEN_SAMPLES=()

N_SAMPLES=0

while IFS= read -r id || [[ -n "$id" ]]; do
    # Remove Windows carriage returns (CRLF).
    id="${id%$'\r'}"
    # Skip empty lines and comments.
    [[ -z "$id" || "$id" == \#* ]] && continue
    # Validate SRA accession format.
    if [[ ! "$id" =~ ^SRR[0-9]+$ ]]; then
        echo "ERROR: Invalid SRA accession: $id" >&2
        exit 1
    fi
    # Reject duplicated samples.
    if [[ -n "${SEEN_SAMPLES[$id]+found}" ]]; then
        echo "ERROR: Duplicate sample: $id" >&2
        exit 1
    fi

    SEEN_SAMPLES["$id"]=1

    # Construct the expected STAR BAM path.
    bam="${ALNDIR}/${id}_Aligned.sortedByCoord.out.bam"
    # Check whether the BAM exists and is not empty.
    if [[ ! -s "$bam" ]]; then
        echo "ERROR: Missing or empty BAM: $bam" >&2
        exit 1
    fi

    # Check basic BAM integrity.
    # quickcheck validates the BAM header and EOF marker,
    # but does not inspect every alignment record.
    if ! "$SAMTOOLS" quickcheck -v "$bam"; then
        echo "ERROR: BAM failed integrity check: $bam" >&2
        exit 1
    fi

    # Confirm that the BAM header declares coordinate sorting.
    if ! "$SAMTOOLS" view -H "$bam" |
         awk -F '\t' '
             $1 == "@HD" {
                 for (i = 2; i <= NF; i++) {
                     if ($i == "SO:coordinate") found = 1
                 }
             }
             END { exit !found }
         '; then
        echo "ERROR: BAM is not marked coordinate-sorted: $bam" >&2
        exit 1
    fi

    # Add the validated BAM to the input list.
    printf '%s\n' "$bam" >> "$TEMP_BAMLIST"
    N_SAMPLES=$((N_SAMPLES + 1))
    echo "[OK] $id"
done < "$RUNS_FILE"


# ------------------------------------------------------------
# 6. Verify the number of selected samples
# ------------------------------------------------------------
echo
echo "Samples found: $N_SAMPLES"

if [[ "$N_SAMPLES" -ne "$EXPECTED_SAMPLES" ]]; then
    echo "ERROR: Expected $EXPECTED_SAMPLES samples," \
         "but found $N_SAMPLES." >&2
    exit 1
fi

echo "[OK] All selected BAMs passed basic validation."


# ------------------------------------------------------------
# 7. Optional validation-only mode
# ------------------------------------------------------------
# This mode checks paths, sample identifiers and BAM files,
# without launching PsiCLASS or replacing existing outputs.

if [[ "$CHECK_ONLY" == "1" ]]; then
    echo
    echo "Validation completed successfully."
    echo "PsiCLASS was not executed."
    exit 0
fi


# ------------------------------------------------------------
# 8. Protect existing analysis results
# ------------------------------------------------------------
# Prevent accidental overwriting of existing assemblies,
# logs, or BAM lists.

# Module I and Module II should have separate output folders.
for existing in \
    "$BAMLIST" \
    "$LOGFILE" \
    "$WORKDIR"/psiclass*.gtf
do
    if [[ -e "$existing" ]]; then
        echo "ERROR: Existing PsiCLASS output detected:" >&2
        echo "$existing" >&2
        echo "Review the output directory before rerunning." >&2
        exit 1
    fi
done

# Publish the validated BAM list.
mv "$TEMP_BAMLIST" "$BAMLIST"


# ------------------------------------------------------------
# 9. Run PsiCLASS transcript assembly
# ------------------------------------------------------------
echo
echo "Starting PsiCLASS transcript assembly..."

# PsiCLASS parameters:
# --lb : File containing BAM paths, one per line.
# -o   : Prefix for generated transcript assemblies.
# -p   : Number of threads.
#
# Both standard output and error messages are saved
# to the analysis log.

if ! psiclass \
    --lb "$BAMLIST" \
    -o "$OUTPUT_PREFIX" \
    -p "$THREADS" \
    > "$LOGFILE" 2>&1
then
    echo "ERROR: PsiCLASS execution failed." >&2
    echo "Last log messages:" >&2
    tail -n 30 "$LOGFILE" >&2
    exit 1
fi


# ------------------------------------------------------------
# 10. Verify transcript assembly output
# ------------------------------------------------------------

# Check that PsiCLASS generated the consensus GTF.
if [[ ! -s "$CONSENSUS_GTF" ]]; then
    echo "ERROR: Consensus GTF missing or empty:" >&2
    echo "$CONSENSUS_GTF" >&2
    exit 1
fi

echo
echo "========================================"
echo "PsiCLASS completed successfully!"
echo "========================================"
echo "Samples analyzed: $N_SAMPLES"
echo "Consensus GTF: $CONSENSUS_GTF"
echo "BAM list: $BAMLIST"
echo "Log file: $LOGFILE"
echo "========================================"
