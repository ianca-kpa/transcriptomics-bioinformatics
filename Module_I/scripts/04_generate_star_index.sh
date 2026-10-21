
#!/usr/bin/env bash
set -e

## STAR genome index generation

DATADIR=/mnt/d/Bioinformatics/transcriptomics-data/Module_I/reference
mkdir -p "${DATADIR}/star_index"

STAR --runThreadN 2 \
     --runMode genomeGenerate \
     --genomeDir ${DATADIR}/star_index \
     --genomeFastaFiles ${DATADIR}/genome/GRCm39.primary_assembly.genome.fa \
     --genomeSAsparseD 16 \
     --genomeSAindexNbases 12

echo "STAR genome index completed!"
