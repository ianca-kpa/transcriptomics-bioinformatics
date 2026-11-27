
# ============================================================
# DESeq2 differential expression analysis
# Comparison: 4 months vs 1 month
# Adjusted for biological sex
# Annotation: PsiCLASS (genome-wide)
# ============================================================

# Required packages
packages <- c(
    "DESeq2", "GenomicFeatures", "txdbmaker", "Rsamtools",
    "GenomicAlignments", "EnhancedVolcano", "ggplot2"
)

for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        stop(paste("Missing R package:", pkg))
    }
    library(pkg, character.only = TRUE)
}

cat("All required R packages are available!\n")


# ============================================================
# 1. Directories and configuration
# ============================================================

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", script_args[grepl("^--file=", script_args)])

SCRIPT_PATH <- normalizePath(file_arg, winslash = "/", mustWork = TRUE)
SCRIPT_DIR <- dirname(SCRIPT_PATH)
PROJECT_DIR <- dirname(SCRIPT_DIR)
config_file <- file.path(PROJECT_DIR, "configs", "paths.sh")

if (!file.exists(config_file)) {
    stop(paste("Configuration file not found:", config_file))
}

base_dir_command <- sprintf(
    "unset BASEDIR; source %s && printf '%%s' \"$BASEDIR\"",
    shQuote(config_file)
)
BASELINE <- system2(
    "bash",
    c("-c", shQuote(base_dir_command)),
    stdout = TRUE,
    stderr = TRUE
)

base_dir_status <- attr(BASELINE, "status")
if (!is.null(base_dir_status) && base_dir_status != 0L) {
    stop("Failed to load BASEDIR from configs/paths.sh; check the configuration file.")
}

BASEDIR <- paste(BASELINE, collapse = "")

if (!nzchar(BASEDIR) || !dir.exists(BASEDIR)) {
    stop("configs/paths.sh must define BASEDIR as an existing directory.")
}
BASEDIR <- normalizePath(BASEDIR, winslash = "/", mustWork = TRUE)

DATADIR <- file.path(BASEDIR, "STAR/BAM")
OUTDIR <- file.path(BASEDIR, "DESeq2_genomewide/sex_adjusted")
RESULTSDIR <- file.path(OUTDIR, "Results")

ANNOT <- file.path(BASEDIR, "PsiCLASS/WithGeneNames/psiclass_vote.withStrand.gtf")

# Existing metadata source (read-only input)
PHENO <- file.path(BASEDIR, "results-proj1/DESeq2/pheno_1mo_4mo.w-sex.csv")


# ============================================================
# 2. Validate input files
# ============================================================

if (!file.exists(ANNOT)) {
    stop(paste("PsiCLASS annotation not found:", ANNOT))
}

if (!file.exists(PHENO)) {
    stop(paste("Metadata file not found:", PHENO))
}

if (!dir.exists(DATADIR)) {
    stop(paste("STAR BAM directory not found:", DATADIR))
}


# ============================================================
# 3. Select genome-wide STAR BAM files
# ============================================================

sample_ids <- c(
    "SRR3734796", "SRR3734797", "SRR3734798",
    "SRR3734820", "SRR3734824", "SRR3734825"
)

fls <- file.path(DATADIR, paste0(sample_ids, "_Aligned.sortedByCoord.out.bam"))

stopifnot(length(fls) == 6)
stopifnot(all(file.exists(fls)))

cat("\nSelected STAR BAM files:\n")
print(basename(fls))


# ============================================================
# 4. Read and validate sample metadata
# ============================================================

sampleInfo <- read.csv(PHENO, stringsAsFactors = FALSE)

required_columns <- c("sampleID", "condition", "sex")

if (!all(required_columns %in% names(sampleInfo))) {
    stop("Metadata must contain sampleID, condition and sex columns.")
}

if (anyDuplicated(sampleInfo$sampleID)) {
    stop("Duplicated sample IDs found in metadata.")
}

if (!setequal(basename(fls), sampleInfo$sampleID)) {
    stop("Metadata sample IDs must match the selected STAR BAM basenames.")
}

# Match metadata order to BAM file order
seIdx <- match(basename(fls), sampleInfo$sampleID)

stopifnot(!anyNA(seIdx))
sampleInfo <- sampleInfo[seIdx, , drop = FALSE]

# Define experimental conditions
sampleInfo$condition <- factor(sampleInfo$condition, levels = c("1mo", "4mo"))

# Accept common representations of biological sex
sex_values <- trimws(tolower(sampleInfo$sex))
sex_values[sex_values %in% c("f", "female")] <- "F"
sex_values[sex_values %in% c("m", "male")] <- "M"

if (anyNA(sex_values) || !all(sex_values %in% c("F", "M"))) {
    stop("Unexpected sex values in metadata. Check the sex column.")
}

sampleInfo$sex <- factor(sex_values, levels = c("F", "M"))

stopifnot(!anyNA(sampleInfo$condition))
stopifnot(!anyNA(sampleInfo$sex))
stopifnot(nlevels(droplevels(sampleInfo$sex)) == 2)


# ============================================================
# 5. Validate experimental design
# ============================================================

# Expected distribution:
#            F   M
# 1mo        3   0
# 4mo        1   2

expected_distribution <- matrix(
    c(3, 1, 0, 2),
    nrow = 2,
    dimnames = list(condition = c("1mo", "4mo"), sex = c("F", "M"))
)

observed_distribution <- table(sampleInfo$condition, sampleInfo$sex)

if (!all(observed_distribution == expected_distribution)) {
    stop("Unexpected condition/sex distribution.")
}

cat("\nSample distribution by condition and sex:\n")
print(observed_distribution)


# ============================================================
# 6. Protect existing results
# ============================================================

# Refuse to overwrite any existing output directory
if (file.exists(OUTDIR) || dir.exists(OUTDIR)) {
    stop(paste("Output directory already exists; refusing to overwrite:", OUTDIR))
}


# ============================================================
# 7. Import genome-wide PsiCLASS annotation
# ============================================================

cat("\nImporting PsiCLASS annotation...\n")

hse <- txdbmaker::makeTxDbFromGFF(ANNOT, format = "gtf")
exonsByGene <- exonsBy(hse, by = "gene")


# ============================================================
# 8. Count reads overlapping genes
# ============================================================

cat("\nCounting reads overlapping genes...\n")

bamLst <- BamFileList(fls, yieldSize = 100000)

se <- summarizeOverlaps(
    features = exonsByGene,
    reads = bamLst,
    mode = "Union",
    singleEnd = TRUE,
    ignore.strand = TRUE,
    fragments = FALSE
)

colnames(se) <- basename(fls)

stopifnot(identical(colnames(se), sampleInfo$sampleID))

# Add validated metadata
colData(se) <- cbind(colData(se), sampleInfo)

cat("Genes counted:", nrow(se), "\n")
cat("Samples counted:", ncol(se), "\n")


# ============================================================
# 9. Create output directories
# ============================================================

# Create the parent directory if necessary
parent_dir <- dirname(OUTDIR)

if (!dir.exists(parent_dir)) {
    dir.create(parent_dir, recursive = TRUE, showWarnings = FALSE)
}

if (!dir.exists(parent_dir)) {
    stop(paste("Could not create parent directory:", parent_dir))
}

# Atomic directory creation: fail if it already exists
if (!dir.create(OUTDIR, recursive = FALSE, showWarnings = FALSE)) {
    stop(paste("Could not create output directory; it may already exist:", OUTDIR))
}

if (!dir.create(RESULTSDIR, showWarnings = FALSE)) {
    stop(paste("Could not create results directory:", RESULTSDIR))
}


# ============================================================
# 10. Export raw counts
# ============================================================

raw_counts_file <- file.path(
    RESULTSDIR,
    "raw_counts.1mo_4mo.A.psiclass.csv"
)
write.csv(assay(se), raw_counts_file)


# ============================================================
# 11. Create DESeq2 dataset
# ============================================================

# Adjust for sex while testing the condition effect
dds <- DESeqDataSet(se, design = ~ sex + condition)

# Define reference levels
dds$sex <- relevel(dds$sex, ref = "F")
dds$condition <- relevel(dds$condition, ref = "1mo")

cat("\nDESeq2 design:\n")
print(design(dds))


# ============================================================
# 12. Run differential expression analysis
# ============================================================

cat("\nRunning DESeq2...\n")

dds <- DESeq(dds)

# Comparison: 4mo vs 1mo, adjusted for sex
# Positive log2FC: higher expression at 4mo
# Negative log2FC: higher expression at 1mo

res <- results(dds, contrast = c("condition", "4mo", "1mo"))


# ============================================================
# 13. Export normalized counts
# ============================================================

write.csv(
    counts(dds, normalized = TRUE),
    file.path(RESULTSDIR, "normalized_counts.1mo_4mo.A.psiclass.csv")
)


# ============================================================
# 14. Export differential expression results
# ============================================================

res_df <- as.data.frame(res)

write.csv(
    res_df,
    file.path(RESULTSDIR, "results_1mo_4mo.A.psiclass.csv")
)

# Filter by unadjusted p-value
res_pval <- res_df[
    !is.na(res_df$pvalue) & res_df$pvalue <= 0.05,
    ,
    drop = FALSE
]

write.csv(
    res_pval,
    file.path(RESULTSDIR, "results_1mo_4mo.pval05.A.psiclass.csv")
)

# Filter by adjusted p-value (FDR)
res_qval <- res_df[
    !is.na(res_df$padj) & res_df$padj <= 0.05,
    ,
    drop = FALSE
]

write.csv(
    res_qval,
    file.path(RESULTSDIR, "results_1mo_4mo.qval05.A.psiclass.csv")
)


# ============================================================
# 15. Variance-stabilizing transformation and PCA
# ============================================================

cat("\nGenerating PCA plot...\n")

vsd <- vst(dds, blind = TRUE)

pca_plot <- plotPCA(vsd, intgroup = "condition")

ggsave(
    filename = file.path(OUTDIR, "PCA_1mo_4mo_byCondition.A.png"),
    plot = pca_plot,
    width = 7,
    height = 5,
    dpi = 300
)

# Export PCA coordinates
pca <- prcomp(t(assay(vsd)))

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
    file.path(OUTDIR, "1mo_4mo_pca_table.A.psiclass.csv"),
    row.names = FALSE
)


# ============================================================
# 16. Volcano plot
# ============================================================

cat("\nGenerating volcano plot...\n")

png(
    filename = file.path(OUTDIR, "volcano_1mo_4mo.A.psiclass.png"),
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


# ============================================================
# 17. Analysis summary
# ============================================================

session_info_file <- file.path(OUTDIR, "sessionInfo.txt")
summary_file <- file.path(OUTDIR, "analysis_summary.txt")
provenance_file <- file.path(OUTDIR, "analysis_provenance.txt")

summary_lines <- c(
    "DESeq2 analysis adjusted for sex completed.",
    paste("Total genes:", nrow(res_df)),
    paste("Genes with p-value <= 0.05:", nrow(res_pval)),
    paste("Genes with padj <= 0.05:", nrow(res_qval)),
    paste(
        "Higher expression at 4mo:",
        sum(res_qval$log2FoldChange > 0, na.rm = TRUE)
    ),
    paste(
        "Higher expression at 1mo:",
        sum(res_qval$log2FoldChange < 0, na.rm = TRUE)
    ),
    "Model: ~ sex + condition (sex reference: F; condition reference: 1mo)",
    "Contrast: condition 4mo vs 1mo",
    paste("Results saved to:", RESULTSDIR)
)
writeLines(summary_lines, summary_file)
writeLines(capture.output(sessionInfo()), session_info_file)

output_paths <- c(
    RESULTSDIR,
    file.path(
        RESULTSDIR,
        c(
            "raw_counts.1mo_4mo.A.psiclass.csv",
            "normalized_counts.1mo_4mo.A.psiclass.csv",
            "results_1mo_4mo.A.psiclass.csv",
            "results_1mo_4mo.pval05.A.psiclass.csv",
            "results_1mo_4mo.qval05.A.psiclass.csv"
        )
    ),
    file.path(
        OUTDIR,
        c(
            "PCA_1mo_4mo_byCondition.A.png",
            "1mo_4mo_pca_table.A.psiclass.csv",
            "volcano_1mo_4mo.A.psiclass.png",
            "sessionInfo.txt",
            "analysis_summary.txt",
            "analysis_provenance.txt"
        )
    )
)

provenance_lines <- c(
    "DESeq2 sex-adjusted analysis provenance",
    paste("Execution time (UTC):", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    paste("R version:", R.version.string),
    paste("DESeq2 version:", as.character(packageVersion("DESeq2"))),
    paste("Script:", SCRIPT_PATH),
    paste("Script MD5:", unname(tools::md5sum(SCRIPT_PATH))),
    paste("Metadata:", normalizePath(PHENO, winslash = "/", mustWork = TRUE)),
    paste("Metadata MD5:", unname(tools::md5sum(PHENO))),
    paste("GTF:", normalizePath(ANNOT, winslash = "/", mustWork = TRUE)),
    paste("GTF MD5:", unname(tools::md5sum(ANNOT))),
    paste("Raw-count matrix:", normalizePath(raw_counts_file, winslash = "/", mustWork = TRUE)),
    paste("Raw-count matrix MD5:", unname(tools::md5sum(raw_counts_file))),
    paste("Design:", paste(deparse(design(dds)), collapse = "")),
    "Contrast: condition 4mo vs 1mo",
    "Counting: mode=Union; singleEnd=TRUE; ignore.strand=TRUE; fragments=FALSE",
    "Samples (in BAM order):",
    capture.output(print(sampleInfo[, required_columns, drop = FALSE])),
    "BAM inputs:",
    paste0("  ", normalizePath(fls, winslash = "/", mustWork = TRUE)),
    "Execution summary:",
    summary_lines,
    "Output paths:",
    paste0("  ", output_paths)
)
writeLines(provenance_lines, provenance_file)

cat("\n========================================\n")
cat(paste(summary_lines, collapse = "\n"), "\n")
cat("Provenance saved to:", provenance_file, "\n")
