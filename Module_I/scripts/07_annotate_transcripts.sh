
#!/usr/bin/env bash
set -euo pipefail

# Module I - Annotation of PsiCLASS transcript assemblies
# Adapted from the Johns Hopkins University course

BASEDIR="/mnt/d/Bioinformatics/transcriptomics-data/Module_I"
WORKDIR="${BASEDIR}/results/PsiCLASS"

# Reference gene annotation (to be downloaded/verified)
ANNOT="${BASEDIR}/reference/annotation/gencode.vM30.annotation.gtf"

# Course helper scripts (locations to verify)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

ADDSTRAND="${PROJECT_DIR}/docs/add-strand"

# The add-genename executable must be located or installed.
ADDGENENAME="add-genename"

# Add annotated gene names
mkdir -p "${WORKDIR}/WithGeneNames"

# Add gene names to the PsiCLASS GTF file using the reference annotation
"${ADDGENENAME}" \
    "${ANNOT}" \
    "${WORKDIR}/psiclass_gtf.list" \
    -o "${WORKDIR}/WithGeneNames"

# Add strand information using the reference annotation
# -r: Reference annotation GTF file
# -i: Input GTF file (PsiCLASS output)
# <: Input redirection to read from the specified file
# >: Output redirection to write to the specified file

"${ADDSTRAND}" \
    "${ANNOT}" \
    -r \ 
    -i \ 
    < "${WORKDIR}/WithGeneNames/psiclass_vote.gtf" \ 
    > "${WORKDIR}/WithGeneNames/psiclass_vote.withStrand.gtf"  

echo "PsiCLASS transcript annotation completed!"
