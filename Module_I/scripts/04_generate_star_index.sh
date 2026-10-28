
#!/usr/bin/env bash
set -e

## STAR genome index generation

DATADIR=/mnt/d/Bioinformatics/transcriptomics-data/Module_I/reference
mkdir -p "${DATADIR}/star_index"

# Generate STAR genome index
# --runThreadN: Number of threads to use
# --runMode genomeGenerate: Generate genome index
# --genomeDir: Directory to store the genome index
# --genomeFastaFiles: Path to the genome FASTA file
# --genomeSAsparseD: Sparse index parameter (default: 16)
# --genomeSAindexNbases: Number of bases for the suffix array index (default

STAR --runThreadN 2 \
     --runMode genomeGenerate \
     --genomeDir ${DATADIR}/star_index \
     --genomeFastaFiles ${DATADIR}/genome/GRCm39.primary_assembly.genome.fa \
     --genomeSAsparseD 16 \
     --genomeSAindexNbases 12

echo "STAR genome index completed!"
