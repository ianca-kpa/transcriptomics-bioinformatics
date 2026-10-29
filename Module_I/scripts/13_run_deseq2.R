
# Module I - Differential gene expression with DESeq2
# Comparison: 1-month versus 4-month mouse hippocampus samples
#
# This script performs a differential expression analysis using DESeq2 to compare gene expression between 1-month-old and 4-month-old mouse hippocampus samples.

# Required packages for DESeq2 analysis
packages <- c(
    "DESeq2",
    "GenomicFeatures",
    "Rsamtools",
    "GenomicAlignments",
    "EnhancedVolcano",
    "txdbmaker"
)

# Check if each package is installed
for (pkg in packages) {

    if (!requireNamespace(pkg, quietly=TRUE)) {
        stop(paste("Package not installed:", pkg))
    }

    library(pkg, character.only=TRUE)
}

cat("All required R packages are available!\n")

# --------------------------------------------------
# 1. Input files
# --------------------------------------------------

BASEDIR <- "/mnt/d/Bioinformatics/transcriptomics-data/Module_I"

# BAM files sorted by read name
DATADIR <- file.path(BASEDIR, "results/Sorted")

# PsiCLASS transcript annotation
ANNOT <- file.path(
    BASEDIR,
    "results/PsiCLASS/WithGeneNames/psiclass_vote.withStrand.gtf"
)

# Sample groups:
# 1mo: SRR3734796, SRR3734797, SRR3734798
# 4mo: SRR3734820, SRR3734824, SRR3734825

# List all BAM files
fls <- list.files(
    DATADIR,
    pattern="\\.bam$",
    full.names=TRUE
)

# Check that exactly six BAM files were found
stopifnot(length(fls) == 6)
print(basename(fls))

# Select the six BAM files
fls_subset <- fls[1:6]

# --------------------------------------------------
# 2. Import gene annotation
# --------------------------------------------------

# Create a transcript database from the PsiCLASS GTF
hse <- txdbmaker::makeTxDbFromGFF(
    ANNOT,
    format="gtf"
)

# Group exonic regions by gene
exonsByGene <- exonsBy(hse, by="gene")

# --------------------------------------------------
# 3. Count reads overlapping gene exons
# --------------------------------------------------

# Create a list of BAM files
bamLst <- BamFileList(
    fls_subset,
    yieldSize=100000
)

# Count reads assigned to each gene
# singleEnd=TRUE: single-end sequencing data
# ignore.strand=TRUE: count without strand specificity
# fragments=FALSE: no paired-end fragment handling

se <- summarizeOverlaps(
    features=exonsByGene,
    reads=bamLst,
    mode="Union",
    singleEnd=TRUE,
    ignore.strand=TRUE,
    fragments=FALSE
)

# Use BAM filenames as sample identifiers
colnames(se) <- basename(fls_subset)

# Save raw gene counts
write.csv(
    assay(se),
    file="raw_counts.1mo_4mo.psiclass.csv"
)

# --------------------------------------------------
# 4. Read sample metadata
# --------------------------------------------------

# Load sample information
sampleInfo <- read.csv("pheno_1mo_4mo.csv")

# Match BAM filenames with metadata
seIdx <- match(colnames(se), sampleInfo$sampleID)

# Check for missing, duplicated, or unexpected samples
stopifnot(!anyNA(seIdx))
stopifnot(!anyDuplicated(sampleInfo$sampleID))
stopifnot(setequal(colnames(se), sampleInfo$sampleID))

# Arrange metadata in the same order as the count matrix
sampleInfo <- sampleInfo[seIdx, , drop=FALSE]

# Set the reference condition
sampleInfo$condition <- factor(
    sampleInfo$condition,
    levels=c("1mo", "4mo")
)

stopifnot(!anyNA(sampleInfo$condition))

# Add metadata to the count object
colData(se) <- cbind(colData(se), sampleInfo)

# Display sample information
print(as.data.frame(colData(se)))

# --------------------------------------------------
# 5. Differential expression with DESeq2
# --------------------------------------------------

# Create the DESeq2 dataset
dds <- DESeqDataSet(
    se,
    design=~ condition
)

# Set 1 month as the reference group
dds$condition <- relevel(dds$condition, ref="1mo")

# Run differential expression analysis
dds <- DESeq(dds)

# Display the available comparison
print(resultsNames(dds))

# Extract results: 4 months versus 1 month
res <- results(
    dds,
    contrast=c("condition", "4mo", "1mo")
)

# Filter results by p-value and adjusted p-value
res_pval0.05 <- res[
    which(res$pvalue <= 0.05),
]

res_qval0.05 <- res[
    which(res$padj <= 0.05),
]

# --------------------------------------------------
# 6. Export results
# --------------------------------------------------

# Save normalized counts
write.csv(
    round(counts(dds, normalized=TRUE), digits=2),
    file="normalized_counts.1mo_4mo.psiclass.csv"
)

# Save results for all genes
write.csv(
    as.data.frame(res),
    file="results_1mo_4mo.psiclass.csv"
)

# Save genes with p-value <= 0.05
write.csv(
    as.data.frame(res_pval0.05),
    file="results_1mo_4mo.pval05.psiclass.csv"
)

# Save genes with adjusted p-value <= 0.05
write.csv(
    as.data.frame(res_qval0.05),
    file="results_1mo_4mo.qval05.psiclass.csv"
)

# Print a summary of the results
print(summary(res))

cat("Genes with p-value <= 0.05:", nrow(res_pval0.05), "\n")
cat("Genes with padj <= 0.05:", nrow(res_qval0.05), "\n")

# --------------------------------------------------
# 7. PCA visualization
# --------------------------------------------------

# Variance-stabilizing transformation
vsd <- vst(dds, blind=TRUE)

# Save PCA plot
png(
    "PCA_1mo_4mo_byCondition.png",
    width=720,
    height=720
)

print(plotPCA(vsd, intgroup="condition"))

dev.off()

# Calculate PCA coordinates
vsd_mat <- assay(vsd)
pca <- prcomp(t(vsd_mat))

# Save PCA coordinates
write.csv(
    pca$x,
    file="1mo_4mo_pca_table.psiclass.csv"
)

print(summary(pca))

# --------------------------------------------------
# 8. Volcano plot
# --------------------------------------------------

# Use empty labels to avoid overcrowding
gene_labels <- rep("", nrow(res))

png(
    "volcano_1mo_4mo.psiclass.png",
    width=720,
    height=720
)

myV <- EnhancedVolcano(
    as.data.frame(res),
    lab=gene_labels,
    title="1mo_4mo.psiclass",
    x="log2FoldChange",
    y="padj",
    pCutoff=0.05,
    FCcutoff=1.5,
    labSize=3.0,
    titleLabSize=24,
    cutoffLineWidth=0.2,
    ylab=bquote(~-Log[10]~italic(Padj)),
    xlim=c(-10, 10)
)

print(myV)

dev.off()

cat("DESeq2 analysis completed!\n")
