
#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

source "${PROJECT_DIR}/configs/paths.sh"

# Directory containing BAM files sorted by read name
WORKDIR="$SORTED_BAM_DIR"

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

    # Require the BAM header to declare queryname sorting.
    if ! /usr/bin/samtools view -H "$BAM" |
         awk -F '\t' '
             $1 == "@HD" {
                 for (i = 2; i <= NF; i++) {
                     if ($i == "SO:queryname") found = 1
                 }
             }
             END { exit !found }
         '; then
        echo "ERROR: BAM is not marked queryname-sorted: $BAM" >&2
        exit 1
    fi

    echo "Sort order: queryname"

    echo "Finished checking $i"
    echo "---------------------------"
done

echo "All six BAM files passed the basic checks!"
