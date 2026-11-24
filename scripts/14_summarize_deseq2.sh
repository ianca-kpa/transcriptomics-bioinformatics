
#!/usr/bin/env bash
set -euo pipefail

# This script summarizes the DESeq2 differential expression results.
# Comparison: 4-month-old vs 1-month-old samples.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

source "${PROJECT_DIR}/configs/paths.sh"

WORKDIR="$DESEQ2_RESULTS_DIR"

# awk reads the CSV files and counts the number of lines.
# NR represents the total number of lines.
# NR-1 excludes the header from the gene count.

{
    echo "DESeq2 - Differential Expression Summary"
    echo "Comparison: 4mo vs 1mo"
    echo "----------------------------------------"

    awk 'END {print "Total genes:", NR-1}' \
        "${WORKDIR}/results_1mo_4mo.psiclass.csv"

    awk 'END {print "Genes with p-value <= 0.05:", NR-1}' \
        "${WORKDIR}/results_1mo_4mo.pval05.psiclass.csv"

    echo "Genes with q-value <= 0.05:"
    awk 'END {print "Genes with padj <= 0.05:", NR-1}' \
        "${WORKDIR}/results_1mo_4mo.qval05.psiclass.csv"

    # Count significant upregulated and downregulated genes
    export DESEQ2_RESULTS_CSV="${WORKDIR}/results_1mo_4mo.psiclass.csv"
    Rscript -e '
    res <- read.csv(Sys.getenv("DESEQ2_RESULTS_CSV"))

    sig <- res[!is.na(res$padj) & res$padj <= 0.05, ]

    cat("Upregulated genes:", sum(sig$log2FoldChange > 0), "\n")
    cat("Downregulated genes:", sum(sig$log2FoldChange < 0), "\n")
    '

} | tee "${WORKDIR}/deseq2_summary.txt"
