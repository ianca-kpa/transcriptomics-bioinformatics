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
  $STAR --runThreadN 2 --genomeDir $STAR2IDX/ \
        --readFilesIn $DATADIR/${i}.fastq \
        --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif \
        --outFileNamePrefix ${WORKDIR}/${i}_

done

echo "All alignments completed!"
