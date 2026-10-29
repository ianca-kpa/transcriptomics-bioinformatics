
#!/usr/bin/env bash
set -e

# Directory containing BAM files sorted by read name
WORKDIR="/mnt/d/Bioinformatics/transcriptomics-data/Module_I/results/Sorted"

# Check all six samples
for i in SRR3734796 SRR3734797 SRR3734798 SRR3734820 SRR3734824 SRR3734825
do
    BAM="${WORKDIR}/${i}.so.bam"

    echo "Checking $i..."

    # Check if the BAM file exists
    if [ ! -f "$BAM" ]; then
        echo "ERROR: BAM file not found: $BAM"
        exit 1
    fi

    # Check basic BAM file integrity
    /usr/bin/samtools quickcheck "$BAM"
    # quickcheck is a command-line tool that checks the integrity of BAM files. 
    # It verifies that the file is not truncated and that the header and alignment sections are consistent.

    # Check the BAM sorting order
    /usr/bin/samtools view -H "$BAM" | grep '^@HD'
    # The command extracts the header lines from the BAM file and searches for the '@HD' line, which contains information about the file's header.

    echo "Finished checking $i"
    echo "---------------------------"
done

echo "All six BAM files passed the basic checks!"
