#!/usr/bin/env bash

# Module I - Cuffdiff2 differential expression summary

# Comparison: Month 4 versus Month 1
# Positive log2FC: higher expression at Month 4
# Negative log2FC: lower expression at Month 4

FILE="results/Cuffdiff/gene_exp.diff"
ISOFORM="results/Cuffdiff/isoform_exp.diff"

# Count significant genes and isoforms
echo "Significant genes:"
grep -c $'\tyes$' "$FILE"

echo "Significant isoforms:"
grep -c $'\tyes$' "$ISOFORM"

# Count upregulated genes
echo "Upregulated genes:"
awk -F '\t' '$14 == "yes" && $10 > 0' "$FILE" | wc -l

# Count downregulated genes
echo "Downregulated genes:"
awk -F '\t' '$14 == "yes" && $10 < 0' "$FILE" | wc -l

# Count genes with p-value < 0.05 and |log2FC| > 1
echo "Genes with p < 0.05 and |log2FC| > 1:"
awk -F '\t' '$7 == "OK" && $12 < 0.05 && ($10 < -1 || $10 > 1)' "$FILE" | wc -l