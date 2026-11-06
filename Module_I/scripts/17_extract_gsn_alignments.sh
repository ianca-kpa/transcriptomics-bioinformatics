
#!/usr/bin/env bash
set -euo pipefail

# Directories
BASEDIR="/mnt/d/Bioinformatics/transcriptomics-data/Module_I"

BAMDIR="${BASEDIR}/results/STAR/BAM"
WORKDIR="${BASEDIR}/results/Visualization"

# Annotation files
P_ANNOT="${BASEDIR}/results/PsiCLASS/WithGeneNames/psiclass_vote.gtf"
G_ANNOT="${BASEDIR}/reference/annotation/gencode.vM30.annotation.gtf"

mkdir -p "$WORKDIR"
cd "$WORKDIR"

### Gsn: 2:35146392-35197904
# The course uses "2", but my STAR BAM files use "chr2".

# The coordinates were obtained from the GENCODE annotation:
# smallest start coordinate of Gsn gene = 35146392
# grep -w Gsn "$G_ANNOT" | cut -f4 | sort -nk1 | head -1

# largest end coordinate of Gsn gene = 35197904
# grep -w Gsn "$G_ANNOT" | cut -f5 | sort -nrk1 | head -1

# Extract GENCODE gene annotations
grep -w Gsn "$G_ANNOT" > Gsn.gencode.vM30.gtf

# Extract PsiCLASS gene and transcript assemblies
grep -w Gsn "$P_ANNOT" > Gsn.psiclass.gtf

# Index the original coordinate-sorted BAM files
# Useful for STAR alignments that may not have been indexed yet
for i in SRR3734796 SRR3734797 SRR3734798 SRR3734820 SRR3734824 SRR3734825
do
    BAM="${BAMDIR}/${i}_Aligned.sortedByCoord.out.bam"

    if [ ! -f "${BAM}.bai" ]; then
        samtools index "$BAM"
    fi
done

# Extract alignments in the Gsn region
# Useful for visualizing the alignments in IGV or other genome browsers
for i in SRR3734796 SRR3734797 SRR3734798 SRR3734820 SRR3734824 SRR3734825
do
    samtools view -b \
        "${BAMDIR}/${i}_Aligned.sortedByCoord.out.bam" \
        "chr2:35146392-35197904" > "Gsn.${i}.bam"

    samtools index "Gsn.${i}.bam"
done

echo "Gsn visualization files are ready!"
