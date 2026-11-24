
#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

source "${PROJECT_DIR}/configs/paths.sh"

# Directories
BAMDIR="$STAR_BAM_DIR"
WORKDIR="${BASEDIR}/results-proj1/Visualization"

# Annotation files
P_ANNOT="${PSICLASS_DIR}/WithGeneNames/psiclass_vote.gtf"
G_ANNOT="$REFERENCE_GTF"

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
