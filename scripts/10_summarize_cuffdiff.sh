#!/usr/bin/env bash
set -euo pipefail

# Module I - Cuffdiff2 differential expression summary

# Comparison: Month 4 versus Month 1
# Positive log2FC: higher expression at Month 4
# Negative log2FC: lower expression at Month 4

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

source "${PROJECT_DIR}/configs/paths.sh"

FILE="${CUFFDIFF_DIR}/gene_exp.diff"
ISOFORM="${CUFFDIFF_DIR}/isoform_exp.diff"

# Count genes with p-value < 0.05
echo "Genes with p < 0.05:"
awk -F '\t' '$7 == "OK" && $12 < 0.05' "$FILE" | wc -l

# Count significant genes and isoforms
echo "Significant genes (q < 0.05):"
grep -c $'\tyes$' "$FILE"

echo "Significant isoforms (q < 0.05):"
grep -c $'\tyes$' "$ISOFORM"

# Count upregulated genes
echo "Upregulated genes (q < 0.05 and log2FC > 0):"
awk -F '\t' '$14 == "yes" && $10 > 0' "$FILE" | wc -l

# Count downregulated genes
echo "Downregulated genes (q < 0.05 and log2FC < 0):"
awk -F '\t' '$14 == "yes" && $10 < 0' "$FILE" | wc -l

# Count genes with p-value < 0.05 and |log2FC| > 1
echo "Genes with p < 0.05 and |log2FC| > 1:"
awk -F '\t' '$7 == "OK" && $12 < 0.05 && ($10 < -1 || $10 > 1)' "$FILE" | wc -l