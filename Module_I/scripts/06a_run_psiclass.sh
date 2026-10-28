
#!/usr/bin/env bash
set -euo pipefail

# ------------------------------------------------------------
# Module I - Transcript assembly with PsiCLASS
# Adapted from the Johns Hopkins University course
# ------------------------------------------------------------

# 1. Project directories
BASEDIR="/mnt/d/Bioinformatics/transcriptomics-data/Module_I"
ALNDIR="${BASEDIR}/results/STAR/BAM"
WORKDIR="${BASEDIR}/results/PsiCLASS"

# 2. Tools
SAMTOOLS="/usr/bin/samtools"

# 3. Output files
BAMLIST="${WORKDIR}/bamlist.psiclass"
OUTPUT_PREFIX="${WORKDIR}/psiclass"
LOGFILE="${WORKDIR}/psiclass.log"

# 4. Number of processes
THREADS=2

# ------------------------------------------------------------
# 5. Check required tools and directories
# ------------------------------------------------------------
# Check if PsiCLASS is available in the current environment
command -v psiclass >/dev/null 2>&1 || {
    echo "ERROR: PsiCLASS is not available in the current environment." >&2
    exit 1
}
# Check if SAMtools is available and executable
if [[ ! -x "$SAMTOOLS" ]]; then
    echo "ERROR: SAMtools not found: $SAMTOOLS" >&2
    exit 1
fi
# Check if the alignment directory exists
if [[ ! -d "$ALNDIR" ]]; then
    echo "ERROR: BAM directory not found: $ALNDIR" >&2
    exit 1
fi

mkdir -p "$WORKDIR"

# ------------------------------------------------------------
# 6. Generate BAM list
# ------------------------------------------------------------

echo "Generating BAM list..."

# Find all BAM files in the alignment directory, sort them, and save to a list file
# -maxdepth 1: Only search in the specified directory, not subdirectories
# -type f: Only include files (not directories)
# -name "*_Aligned.sortedByCoord.out.bam": Match files with the specified pattern

find "$ALNDIR" \
    -maxdepth 1 \
    -type f \
    -name "*_Aligned.sortedByCoord.out.bam" \
    | sort > "$BAMLIST"

echo "BAM files found:"
cat "$BAMLIST"

# Verify the expected number of samples
N_SAMPLES=$(wc -l < "$BAMLIST")

# Check if the number of samples is as expected (6 in this case)
if [[ "$N_SAMPLES" -ne 6 ]]; then
    echo "ERROR: Expected 6 BAM files, found $N_SAMPLES." >&2
    exit 1
fi

# ------------------------------------------------------------
# 7. Validate BAM files
# ------------------------------------------------------------

echo "Validating BAM files..."

# Check that each BAM file exists and is not empty
while IFS= read -r bam; do
    # Check if the BAM file exists and is not empty
    if [[ ! -s "$bam" ]]; then
        echo "ERROR: Missing or empty BAM: $bam" >&2
        exit 1
    fi
    # Validate the BAM file using samtools quickcheck
    if ! "$SAMTOOLS" quickcheck -v "$bam"; then
        echo "ERROR: BAM failed validation: $bam" >&2
        exit 1
    fi

    echo "OK: $(basename "$bam")"

done < "$BAMLIST"

echo "All 6 BAM files passed basic validation."

# ------------------------------------------------------------
# 8. Transcript assembly with PsiCLASS
# ------------------------------------------------------------

echo "Starting PsiCLASS transcript assembly..."
# Run PsiCLASS with the following parameters:
# --lb: Input BAM list file
# -o: Output prefix for PsiCLASS results
# -p: Number of threads to use
# &> : Redirect standard output and errors to a log file

psiclass \
    --lb "$BAMLIST" \
    -o "$OUTPUT_PREFIX" \
    -p "$THREADS" \
    &> "$LOGFILE"

echo "PsiCLASS transcript assembly completed!"
echo "Results directory: $WORKDIR"
echo "Log file: $LOGFILE"
