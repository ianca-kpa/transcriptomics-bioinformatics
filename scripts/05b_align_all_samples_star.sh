
#!/usr/bin/env bash
set -euo pipefail

## STAR alignment script

STAR=STAR

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

source "${PROJECT_DIR}/configs/paths.sh"

TXT="${PROJECT_DIR}/course_materials/Module_II/metadata/selected_runs_m2.txt"
DATADIR="${BASEDIR}/FASTQ"
WORKDIR="$STAR_DIR"
STAR2IDX="$STAR_INDEX_DIR"

# Create main output directory
# mkdir -p "$WORKDIR"

while IFS= read -r i || [[ -n "$i" ]]; do

  # Remove Windows carriage return if present
  i="${i%$'\r'}"

  [[ -z "$i" ]] && continue  # Skip empty lines
  [[ "$i" == \#* ]] && continue  # Skip comment lines

  echo "Aligning ${i}"

  # Create an output directory for each sample
  mkdir -p "${WORKDIR}/${i}"

  # Check if FASTQ file exists
  if [[ ! -f "${DATADIR}/${i}.fastq" ]]; then
    echo "ERROR: FASTQ file not found: ${i}.fastq" >&2
    exit 1
  fi

  # Run STAR alignment
  # --runThreadN: Number of threads to use
  # --genomeDir: Path to the STAR genome index
  # --readFilesIn: Input FASTQ file
  # --outSAMtype: Output format (BAM, sorted by coordinate)
  # --outSAMstrandField: Include strand information in the output BAM
  # --outFileNamePrefix: Prefix for output files

  time "$STAR" \
    --runThreadN 2 \
    --genomeDir "${STAR2IDX}/" \
    --readFilesIn "${DATADIR}/${i}.fastq" \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMstrandField intronMotif \
    --outFileNamePrefix "${WORKDIR}/${i}/${i}_"

  echo "Completed: ${i}"

done < "$TXT"

echo "All alignments completed!"