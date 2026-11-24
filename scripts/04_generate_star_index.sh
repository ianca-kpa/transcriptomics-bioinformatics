
#!/usr/bin/env bash
set -e

## STAR genome index generation

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

source "${PROJECT_DIR}/configs/paths.sh"

DATADIR="$(dirname "$(dirname "$STAR_INDEX_DIR")")"
mkdir -p "$STAR_INDEX_DIR"

# Generate STAR genome index
# --runThreadN: Number of threads to use
# --runMode genomeGenerate: Generate genome index
# --genomeDir: Directory to store the genome index
# --genomeFastaFiles: Path to the genome FASTA file
# --genomeSAsparseD: Sparse index parameter
# --genomeSAindexNbases: Number of bases for the suffix array index

STAR --runThreadN 2 \
     --runMode genomeGenerate \
     --genomeDir "$STAR_INDEX_DIR" \
     --genomeFastaFiles "${BASEDIR}/reference/genome/GRCm39.primary_assembly.genome.fa" \
     --genomeSAsparseD 16 \
     --genomeSAindexNbases 12

echo "STAR genome index completed!"
