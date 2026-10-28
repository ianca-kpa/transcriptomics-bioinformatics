#!/usr/bin/env bash
set -euo pipefail

## STAR alignment script

STAR=STAR

BASEDIR=/mnt/d/Bioinformatics/transcriptomics-data/Module_I
DATADIR=${BASEDIR}/FASTQ
WORKDIR=${BASEDIR}/results/STAR
STAR2IDX=${BASEDIR}/reference/star_index

mkdir -p $WORKDIR

for i in SRR3734796
do
  echo "Aligning ${i}"

  time \

  # Run STAR alignment
  # --runThreadN: Number of threads to use
  # --genomeDir: Path to the STAR genome index
  # --readFilesIn: Input FASTQ file
  # --outSAMtype: Output format (BAM, sorted by coordinate)
  # --outSAMstrandField: Include strand information in the output BAM
  # --outFileNamePrefix: Prefix for output files
  
  $STAR --runThreadN 2 --genomeDir $STAR2IDX/ \
        --readFilesIn $DATADIR/${i}.fastq \
        --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif \
        --outFileNamePrefix ${WORKDIR}/${i}_

done

echo "All alignments completed!"
