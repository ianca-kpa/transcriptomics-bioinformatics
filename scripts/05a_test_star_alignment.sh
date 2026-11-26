
#!/usr/bin/env bash
set -euo pipefail

## STAR alignment script

STAR=STAR

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

source "${PROJECT_DIR}/configs/paths.sh"

DATADIR="${BASEDIR}/FASTQ"
WORKDIR="${STAR_DIR}/test"
STAR2IDX="$STAR_INDEX_DIR"

mkdir -p "$WORKDIR"

for i in SRR3734796
do
  echo "Aligning ${i}"

  # Run STAR alignment
  # --runThreadN: Number of threads to use
  # --genomeDir: Path to the STAR genome index
  # --readFilesIn: Input FASTQ file
  # --outSAMtype: Output format (BAM, sorted by coordinate)
  # --outSAMstrandField: Include strand information in the output BAM
  # --outFileNamePrefix: Prefix for output files
  
  time \
  $STAR --runThreadN 2 --genomeDir $STAR2IDX/ \
        --readFilesIn $DATADIR/${i}.fastq \
        --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif \
        --outFileNamePrefix "${WORKDIR}/${i}_"

done

echo "All alignments completed!"