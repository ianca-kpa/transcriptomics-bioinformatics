
#!/usr/bin/env bash

# External data directory
BASEDIR="/path/to/transcriptomics-data"

# Alignment files
STAR_DIR="${BASEDIR}/STAR"
STAR_BAM_DIR="${BASEDIR}/STAR/BAM"
STAR_INDEX_DIR="${BASEDIR}/reference/star_index"

# PsiCLASS outputs
PSICLASS_DIR="${BASEDIR}/PsiCLASS"
PSICLASS_M2_DIR="${BASEDIR}/PsiCLASS_Module_II"

# Cuffdiff outputs
CUFFDIFF_DIR="${BASEDIR}/Cuffdiff"

# Reference annotation
REFERENCE_DIR="${BASEDIR}/reference/annotation"
REFERENCE_GTF="${REFERENCE_DIR}/gencode.vM30.annotation.gtf"

# Sorted alignments and differential expression results
SORTED_BAM_DIR="${BASEDIR}/results-proj1/Sorted"
DESEQ2_DIR="${BASEDIR}/results-proj1/DESeq2"
DESEQ2_RESULTS_DIR="${DESEQ2_DIR}/Results"

# rMATS outputs
RMATS_DIR="${BASEDIR}/rMATS"