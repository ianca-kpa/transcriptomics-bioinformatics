
#!/usr/bin/env bash
set -euo pipefail
module load cufflinks/2.2.1

# Module I - Differential expression analysis with Cuffdiff2
# Adapted from the Johns Hopkins University course

# Project directories
BASEDIR="/mnt/d/Bioinformatics/transcriptomics-data/Module_I"
DATADIR="${BASEDIR}/results/STAR/BAM"
WORKDIR="${BASEDIR}/results/Cuffdiff"

# PsiCLASS transcript annotation
ANNOT="${BASEDIR}/results/PsiCLASS/WithGeneNames/psiclass_vote.withStrand.gtf"

# Control group: 1 month
CTRLBAMS="${DATADIR}/SRR3734796_Aligned.sortedByCoord.out.bam,${DATADIR}/SRR3734797_Aligned.sortedByCoord.out.bam,${DATADIR}/SRR3734798_Aligned.sortedByCoord.out.bam"

# Test group: 4 months
TESTBAMS="${DATADIR}/SRR3734820_Aligned.sortedByCoord.out.bam,${DATADIR}/SRR3734824_Aligned.sortedByCoord.out.bam,${DATADIR}/SRR3734825_Aligned.sortedByCoord.out.bam"

# Create output directory
mkdir -p "$WORKDIR"

# Differential expression analysis
# -o: Output directory
# -p: Number of threads
# ANNOT: PsiCLASS transcript annotation (GTF)
# CTRLBAMS: Control group BAM files (1 month)
# TESTBAMS: Test group BAM files (4 months)
# &>: Redirect standard output and errors to a log file

cuffdiff \
    -o "$WORKDIR" \ 
    -p 2 \ 
    "$ANNOT" \ 
    "$CTRLBAMS" \ 
    "$TESTBAMS" \ 
    &> "${WORKDIR}/cuffdiff.log" 

echo "Cuffdiff2 analysis completed!"
