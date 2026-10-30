
# DESeq2 differential expression analysis
# Comparison: 4 months vs 1 month
# Annotation: PsiCLASS

# Required packages
packages <- c(
    "DESeq2",
    "GenomicFeatures",
    "txdbmaker",
    "Rsamtools",
    "GenomicAlignments",
    "EnhancedVolcano",
    "ggplot2"
)

for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        stop(paste("Missing R package:", pkg))
    }
    library(pkg, character.only = TRUE)
}

cat("All required R packages are available!\n")

# Directories
BASEDIR <- "/mnt/d/Bioinformatics/transcriptomics-data/Module_I"

DATADIR <- file.path(BASEDIR, "results/Sorted")
OUTDIR <- file.path(BASEDIR, "results/DESeq2")
RESULTSDIR <- file.path(OUTDIR, "Results")

ANNOT <- file.path(
    BASEDIR,
    "results/PsiCLASS/WithGeneNames/psiclass_vote.withStrand.gtf"
)

PHENO <- file.path(OUTDIR, "pheno_1mo_4mo.csv")

dir.create(RESULTSDIR, recursive = TRUE, showWarnings = FALSE)

# Check input files
if (!file.exists(ANNOT)) {
    stop("PsiCLASS annotation not found")
}

if (!file.exists(PHENO)) {
    stop(paste("Metadata file not found:", PHENO))
}

# Find name-sorted BAM files
fls <- list.files(
    DATADIR,
    pattern = "\\.so\\.bam$",
    full.names = TRUE
)

fls <- sort(fls)

stopifnot(length(fls) == 6)

print(basename(fls))

# Import PsiCLASS annotation
hse <- txdbmaker::makeTxDbFromGFF(
    ANNOT,
    format = "gtf"
)

exonsByGene <- exonsBy(hse, by = "gene")

# Count reads overlapping genes
bamLst <- BamFileList(
    fls,
    yieldSize = 100000
)

se <- summarizeOverlaps(
    features = exonsByGene,
    reads = bamLst,
    mode = "Union",
    singleEnd = TRUE,
    ignore.strand = TRUE,
    fragments = FALSE
)

colnames(se) <- basename(fls)

# Export raw counts
write.csv(
    assay(se),
    file.path(
        RESULTSDIR,
        "raw_counts.1mo_4mo.psiclass.csv"
    )
)

# Read sample metadata
sampleInfo <- read.csv(
    PHENO,
    stringsAsFactors = FALSE
)

stopifnot(all(c("sampleID", "condition") %in% names(sampleInfo)))
stopifnot(!anyDuplicated(sampleInfo$sampleID))

# Match metadata to BAM files
seIdx <- match(
    colnames(se),
    sampleInfo$sampleID
)

stopifnot(!anyNA(seIdx))
stopifnot(setequal(colnames(se), sampleInfo$sampleID))

sampleInfo <- sampleInfo[seIdx, , drop = FALSE]

sampleInfo$condition <- factor(
    sampleInfo$condition,
    levels = c("1mo", "4mo")
)

stopifnot(!anyNA(sampleInfo$condition))

colData(se) <- cbind(
    colData(se),
    sampleInfo
)

# Create the DESeq2 dataset
dds <- DESeqDataSet(
    se,
    design = ~ condition
)

# Set 1 month as the reference group
dds$condition <- relevel(
    dds$condition,
    ref = "1mo"
)

# Run differential expression analysis
dds <- DESeq(dds)

# Comparison: 4mo vs 1mo
res <- results(
    dds,
    contrast = c("condition", "4mo", "1mo")
)

# Export normalized counts
write.csv(
    counts(dds, normalized = TRUE),
    file.path(
        RESULTSDIR,
        "normalized_counts.1mo_4mo.psiclass.csv"
    )
)

# Export all DESeq2 results
res_df <- as.data.frame(res)

write.csv(
    res_df,
    file.path(
        RESULTSDIR,
        "results_1mo_4mo.psiclass.csv"
    )
)

# Filter by unadjusted p-value
res_pval <- res_df[
    !is.na(res_df$pvalue) &
    res_df$pvalue <= 0.05,
    ,
    drop = FALSE
]

write.csv(
    res_pval,
    file.path(
        RESULTSDIR,
        "results_1mo_4mo.pval05.psiclass.csv"
    )
)

# Filter by adjusted p-value
res_qval <- res_df[
    !is.na(res_df$padj) &
    res_df$padj <= 0.05,
    ,
    drop = FALSE
]

write.csv(
    res_qval,
    file.path(
        RESULTSDIR,
        "results_1mo_4mo.qval05.psiclass.csv"
    )
)

# Variance-stabilizing transformation
vsd <- vst(dds, blind = TRUE)

# PCA plot
pca_plot <- plotPCA(
    vsd,
    intgroup = "condition"
)

ggsave(
    filename = file.path(
        OUTDIR,
        "PCA_1mo_4mo_byCondition.png"
    ),
    plot = pca_plot,
    width = 7,
    height = 5,
    dpi = 300
)

# Export PCA coordinates
pca <- prcomp(
    t(assay(vsd))
)

percentVar <- 100 * pca$sdev^2 / sum(pca$sdev^2)

pca_table <- data.frame(
    sampleID = rownames(pca$x),
    condition = as.character(dds$condition),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    row.names = NULL
)

write.csv(
    pca_table,
    file.path(
        OUTDIR,
        "1mo_4mo_pca_table.psiclass.csv"
    ),
    row.names = FALSE
)

# Volcano plot
png(
    filename = file.path(
        OUTDIR,
        "volcano_1mo_4mo.psiclass.png"
    ),
    width = 1400,
    height = 1100,
    res = 150
)

print(
    EnhancedVolcano(
        res_df,
        lab = rownames(res_df),
        x = "log2FoldChange",
        y = "padj",
        pCutoff = 0.05,
        FCcutoff = 1.5,
        title = "DESeq2: 4mo vs 1mo",
        subtitle = "PsiCLASS annotation"
    )
)

dev.off()

# Summary
cat("\nDESeq2 analysis completed!\n")
cat("Total genes:", nrow(res_df), "\n")
cat("Genes with p-value <= 0.05:", nrow(res_pval), "\n")
cat("Genes with padj <= 0.05:", nrow(res_qval), "\n")
cat("Results saved to:", RESULTSDIR, "\n")
