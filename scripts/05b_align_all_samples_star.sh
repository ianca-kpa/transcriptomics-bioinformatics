
#!/usr/bin/env bash
set -euo pipefail

# Module II — STAR alignment

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

source "${PROJECT_DIR}/configs/paths.sh"

RUNS_FILE="${PROJECT_DIR}/course_materials/Module_II/metadata/selected_runs_m2.txt"
FASTQ_DIR="${BASEDIR}/FASTQ"

# Check tools and inputs
command -v STAR >/dev/null || {
    echo "ERROR: STAR not found" >&2
    exit 1
}

command -v samtools >/dev/null || {
    echo "ERROR: samtools not found" >&2
    exit 1
}

[[ -s "$RUNS_FILE" ]] || {
    echo "ERROR: Missing run list: $RUNS_FILE" >&2
    exit 1
}

[[ -s "${STAR_INDEX_DIR}/Genome" ]] || {
    echo "ERROR: STAR index not found" >&2
    exit 1
}

mkdir -p "$STAR_BAM_DIR"

# Align samples
while IFS= read -r id || [[ -n "$id" ]]; do

    id="${id%$'\r'}"
    [[ -z "$id" || "$id" == \#* ]] && continue

    fastq="${FASTQ_DIR}/${id}.fastq"
    bam="${STAR_BAM_DIR}/${id}_Aligned.sortedByCoord.out.bam"
    log="${STAR_BAM_DIR}/${id}_Log.final.out"

    # Do not overwrite existing results
    if [[ -e "$bam" || -e "$log" ]]; then

        if [[ -s "$bam" && -s "$log" ]] &&
           samtools quickcheck -v "$bam"; then

            echo "[SKIP] Existing BAM and log: $id"
            continue
        fi

        echo "ERROR: Incomplete existing STAR output: $id" >&2
        exit 1
    fi

    [[ -s "$fastq" ]] || {
        echo "ERROR: Missing FASTQ: $fastq" >&2
        exit 1
    }

    echo "Aligning $id"

  # Run STAR alignment
  # --runThreadN: Number of threads to use
  # --genomeDir: Path to the STAR genome index
  # --readFilesIn: Input FASTQ file
  # --outSAMtype: Output format (BAM, sorted by coordinate)
  # --outSAMstrandField: Include strand information in the output BAM
  # --outFileNamePrefix: Prefix for output files

    STAR \
        --runThreadN 2 \
        --genomeDir "$STAR_INDEX_DIR" \
        --readFilesIn "$fastq" \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMstrandField intronMotif \
        --outFileNamePrefix "${STAR_BAM_DIR}/${id}_"

    samtools quickcheck -v "$bam"

    [[ -s "$log" ]] || {
        echo "ERROR: Missing STAR log: $log" >&2
        exit 1
    }

    echo "[OK] Alignment completed: $id"

done < "$RUNS_FILE"

echo "STAR alignment step completed."