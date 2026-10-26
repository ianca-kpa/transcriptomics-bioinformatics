
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

command -v psiclass >/dev/null 2>&1 || {
    echo "ERROR: PsiCLASS is not available in the current environment." >&2
    exit 1
}

if [[ ! -x "$SAMTOOLS" ]]; then
    echo "ERROR: SAMtools not found: $SAMTOOLS" >&2
    exit 1
fi

if [[ ! -d "$ALNDIR" ]]; then
    echo "ERROR: BAM directory not found: $ALNDIR" >&2
    exit 1
fi

mkdir -p "$WORKDIR"

# ------------------------------------------------------------
# 6. Generate BAM list
# ------------------------------------------------------------

echo "Generating BAM list..."

find "$ALNDIR" \
    -maxdepth 1 \
    -type f \
    -name "*_Aligned.sortedByCoord.out.bam" \
    | sort > "$BAMLIST"

echo "BAM files found:"
cat "$BAMLIST"

# Verify the expected number of samples
N_SAMPLES=$(wc -l < "$BAMLIST")

if [[ "$N_SAMPLES" -ne 6 ]]; then
    echo "ERROR: Expected 6 BAM files, found $N_SAMPLES." >&2
    exit 1
fi

# ------------------------------------------------------------
# 7. Validate BAM files
# ------------------------------------------------------------

echo "Validating BAM files..."

while IFS= read -r bam; do

    if [[ ! -s "$bam" ]]; then
        echo "ERROR: Missing or empty BAM: $bam" >&2
        exit 1
    fi

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

psiclass \
    --lb "$BAMLIST" \
    -o "$OUTPUT_PREFIX" \
    -p "$THREADS" \
    &> "$LOGFILE"

echo "PsiCLASS transcript assembly completed!"
echo "Results directory: $WORKDIR"
echo "Log file: $LOGFILE"
