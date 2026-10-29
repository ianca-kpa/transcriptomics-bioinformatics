
#!/usr/bin/env bash
set -euo pipefail

# This script summarizes the DESeq2 differential expression results.
# Comparison: 4-month-old vs 1-month-old samples.

# awk reads the CSV files and counts the number of lines.
# NR represents the total number of lines.
# NR-1 excludes the header from the gene count.

{
    echo "DESeq2 - Differential Expression Summary"
    echo "Comparison: 4mo vs 1mo"
    echo "----------------------------------------"

    awk 'END {print "Total genes:", NR-1}' \
        results/DESeq2/results_1mo_4mo.psiclass.csv

    awk 'END {print "Genes with p-value <= 0.05:", NR-1}' \
        results/DESeq2/results_1mo_4mo.pval05.psiclass.csv

    awk 'END {print "Genes with padj <= 0.05:", NR-1}' \
        results/DESeq2/results_1mo_4mo.qval05.psiclass.csv

    # Count significant upregulated and downregulated genes
    Rscript -e '
    res <- read.csv("results/DESeq2/results_1mo_4mo.psiclass.csv")

    sig <- res[!is.na(res$padj) & res$padj <= 0.05, ]

    cat("Upregulated genes:", sum(sig$log2FoldChange > 0), "\n")
    cat("Downregulated genes:", sum(sig$log2FoldChange < 0), "\n")
    '

} | tee results/DESeq2/deseq2_summary.txt
