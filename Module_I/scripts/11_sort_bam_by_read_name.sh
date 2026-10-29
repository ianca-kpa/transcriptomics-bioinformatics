
#!/usr/bin/env bash
set -e

# Input and output directories
DATADIR="/mnt/d/Bioinformatics/transcriptomics-data/Module_I/results/STAR/BAM"
WORKDIR="/mnt/d/Bioinformatics/transcriptomics-data/Module_I/results/Sorted"

mkdir -p "$WORKDIR"
cd "$WORKDIR"

# samtools sort: a command-line tool for sorting BAM files by read name.
# The command takes the following parameters:
# -n: sort by read name
# -@ 2: use 2 threads for sorting
# -O BAM: output format is BAM
# -o: output file name
# The variable $i stores the current sample ID in the loop.
# ${i} is used to build the input and output file names.

for i in SRR3734796 SRR3734797 SRR3734798 SRR3734820 SRR3734824 SRR3734825
do
    echo "Sorting $i..."

    /usr/bin/samtools sort \
    -n \
    -@ 2 \
    -O BAM \
    -o "${i}.so.bam" \
    "../STAR/BAM/${i}_Aligned.sortedByCoord.out.bam"
done

# The variable {i} is used to represent the current sample name in the loop. 
# It is used to construct the input and output file names for the sorting operation. 
# For example, if {i} is "SRR3734796", the input file will be "../STAR/BAM/SRR3734796_Aligned.sortedByCoord.out.bam" and the output file will be "SRR3734796.so.bam".

echo "All BAM files have been sorted!"
