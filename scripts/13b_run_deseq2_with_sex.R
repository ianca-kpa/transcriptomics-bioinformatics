
# DESeq2 differential expression analysis
# Comparison: 4 months vs 1 month
# Adjusted for biological sex
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

# Metadata including biological sex
PHENO <- file.path(
    OUTDIR,
    "pheno_1mo_4mo.w-sex.csv"
)

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
        "raw_counts.1mo_4mo.A.psiclass.csv"
    )
)

# Read sample metadata with biological sex
sampleInfo <- read.csv(
    PHENO,
    stringsAsFactors = FALSE
)

stopifnot(
    all(c("sampleID", "condition", "sex") %in% names(sampleInfo))
)

stopifnot(!anyDuplicated(sampleInfo$sampleID))

# Match metadata to BAM files
seIdx <- match(
    colnames(se),
    sampleInfo$sampleID
)

stopifnot(!anyNA(seIdx))
stopifnot(setequal(colnames(se), sampleInfo$sampleID))

sampleInfo <- sampleInfo[seIdx, , drop = FALSE]

# Define experimental variables
sampleInfo$condition <- factor(
    sampleInfo$condition,
    levels = c("1mo", "4mo")
)

# Accept common representations of biological sex
sex_values <- trimws(tolower(sampleInfo$sex))

sex_values[sex_values %in% c("f", "female")] <- "F"
sex_values[sex_values %in% c("m", "male")] <- "M"

if (anyNA(sex_values) ||
    !all(sex_values %in% c("F", "M"))) {
    stop("Unexpected sex values in metadata. Check the sex column.")
}

sampleInfo$sex <- factor(
    sex_values,
    levels = c("F", "M")
)

stopifnot(!anyNA(sampleInfo$condition))
stopifnot(!anyNA(sampleInfo$sex))
stopifnot(nlevels(droplevels(sampleInfo$sex)) == 2)

# Check sample distribution
cat("\nSample distribution by condition and sex:\n")
print(table(sampleInfo$condition, sampleInfo$sex))

# Add metadata to the count object
colData(se) <- cbind(
    colData(se),
    sampleInfo
)

# Create the DESeq2 dataset
# Adjust for sex while testing the effect of condition
dds <- DESeqDataSet(
    se,
    design = ~ sex + condition
)

# Set reference levels
dds$sex <- relevel(
    dds$sex,
    ref = "F"
)

dds$condition <- relevel(
    dds$condition,
    ref = "1mo"
)

# Run differential expression analysis
dds <- DESeq(dds)

# Comparison: 4mo vs 1mo, adjusted for sex
res <- results(
    dds,
    contrast = c("condition", "4mo", "1mo")
)

# Export normalized counts
write.csv(
    counts(dds, normalized = TRUE),
    file.path(
        RESULTSDIR,
        "normalized_counts.1mo_4mo.A.psiclass.csv"
    )
)

# Export all DESeq2 results
res_df <- as.data.frame(res)

write.csv(
    res_df,
    file.path(
        RESULTSDIR,
        "results_1mo_4mo.A.psiclass.csv"
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
        "results_1mo_4mo.pval05.A.psiclass.csv"
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
        "results_1mo_4mo.qval05.A.psiclass.csv"
    )
)

# Variance-stabilizing transformation
vsd <- vst(dds, blind = TRUE)

# PCA plot by condition
pca_plot <- plotPCA(
    vsd,
    intgroup = "condition"
)

ggsave(
    filename = file.path(
        OUTDIR,
        "PCA_1mo_4mo_byCondition.A.png"
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

pca_table <- data.frame(
    sampleID = rownames(pca$x),
    condition = as.character(dds$condition),
    sex = as.character(dds$sex),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    row.names = NULL
)

write.csv(
    pca_table,
    file.path(
        OUTDIR,
        "1mo_4mo_pca_table.A.psiclass.csv"
    ),
    row.names = FALSE
)

# Volcano plot
png(
    filename = file.path(
        OUTDIR,
        "volcano_1mo_4mo.A.psiclass.png"
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
        subtitle = "Adjusted for sex - PsiCLASS"
    )
)

dev.off()

# Summary
cat("\nDESeq2 analysis adjusted for sex completed!\n")
cat("Total genes:", nrow(res_df), "\n")
cat("Genes with p-value <= 0.05:", nrow(res_pval), "\n")
cat("Genes with padj <= 0.05:", nrow(res_qval), "\n")

cat(
    "Upregulated genes:",
    sum(res_qval$log2FoldChange > 0, na.rm = TRUE),
    "\n"
)

cat(
    "Downregulated genes:",
    sum(res_qval$log2FoldChange < 0, na.rm = TRUE),
    "\n"
)

cat("Results saved to:", RESULTSDIR, "\n")
