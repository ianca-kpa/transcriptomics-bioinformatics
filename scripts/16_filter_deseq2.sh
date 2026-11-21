
#!/usr/bin/env bash
set -euo pipefail

# Directories
BASEDIR="/mnt/d/Bioinformatics/transcriptomics-data/Module_I"
WORKDIR="${BASEDIR}/results/DESeq2/Results"

INPUT="${WORKDIR}/WithGeneNames/gn_results_1mo_4mo.psiclass.csv"
OUTPUT="${WORKDIR}/deseq_filter_summary.csv"

# Check input file
if [ ! -f "$INPUT" ]; then
    echo "Error: input file not found: $INPUT" >&2
    exit 1
fi

# Create summary table - columns: baseMean_min, abs_log2FC_min, pvalue_0.05, padj_0.05
echo "baseMean_min,abs_log2FC_min,pvalue_0.05,padj_0.05" > "$OUTPUT"

# Loop through baseMean and abs(log2FC) thresholds
# 0, 10, 15 for baseMean and 0, 0.5, 1 for abs(log2FC)
for base in 0 10 15 
do
    for fc in 0 0.5 1
    do
        # Count genes using unadjusted p-value - p-value < 0.05, baseMean >= base, abs(log2FC) >= fc, and not NA.
        p_count=$(awk -F',' -v b="$base" -v f="$fc" '
            NR > 1 &&
            $3 != "NA" && $4 != "NA" && $7 != "NA" &&
            ($3+0) >= b &&
            (($4+0) <= -f || ($4+0) >= f) &&
            ($7+0) <= 0.05 {count++}
            END {print count+0}
        ' "$INPUT")

        # Count genes using adjusted p-value - padj < 0.05, baseMean >= base, abs(log2FC) >= fc, and not NA.
        q_count=$(awk -F',' -v b="$base" -v f="$fc" '
            NR > 1 &&
            $3 != "NA" && $4 != "NA" && $8 != "NA" &&
            ($3+0) >= b &&
            (($4+0) <= -f || ($4+0) >= f) &&
            ($8+0) <= 0.05 {count++}
            END {print count+0}
        ' "$INPUT")

        echo "$base,$fc,$p_count,$q_count" >> "$OUTPUT"
    done
done

echo "DESeq2 filtering summary:"
column -s, -t "$OUTPUT"

echo
echo "Summary saved to: $OUTPUT"
