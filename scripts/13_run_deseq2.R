# DESeq2 differential expression analysis
# Comparison: 4 months vs 1 month
# Annotation: PsiCLASS

# Required packages
packages <- c("DESeq2", "EnhancedVolcano", "ggplot2")

for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        stop(paste("Missing R package:", pkg))
    }
}

# Resolve project and data paths
script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", script_args[grepl("^--file=", script_args)])
if (length(file_arg) != 1L || !nzchar(file_arg)) {
    stop("Could not determine the script path from commandArgs()")
}

SCRIPT_DIR <- dirname(normalizePath(file_arg, mustWork = TRUE))
PROJECT_DIR <- dirname(SCRIPT_DIR)
config_file <- file.path(PROJECT_DIR, "configs", "paths.sh")
if (!file.exists(config_file)) {
    stop(paste("Path configuration not found:", config_file))
}

base_dir_command <- sprintf(
    "source %s && printf '%%s' \"$BASEDIR\"",
    shQuote(config_file)
)
BASELINE <- system2(
    "bash",
    c("-c", shQuote(base_dir_command)),
    stdout = TRUE,
    stderr = ""
)
if (!is.null(attr(BASELINE, "status")) && attr(BASELINE, "status") != 0) {
    stop("Could not load BASEDIR from configs/paths.sh")
}

BASEDIR <- paste(BASELINE, collapse = "")
if (!nzchar(BASEDIR) || !dir.exists(BASEDIR)) {
    stop(paste("BASEDIR is empty or does not exist:", BASEDIR))
}

COUNTS_FILE <- file.path(
    BASEDIR,
    "results-proj1/DESeq2/Results/raw_counts.1mo_4mo.A.psiclass.csv"
)
PHENO <- file.path(
    BASEDIR,
    "results-proj1/DESeq2/pheno_1mo_4mo.w-sex.csv"
)
OUTDIR <- file.path(BASEDIR, "DESeq2_genomewide/unadjusted")
RESULTSDIR <- file.path(OUTDIR, "Results")

# Validate all inputs before creating output directories
if (!file.exists(COUNTS_FILE)) {
    stop(paste("Raw-count matrix not found:", COUNTS_FILE))
}
if (!file.exists(PHENO)) {
    stop(paste("Metadata file not found:", PHENO))
}

if (file.exists(OUTDIR) || dir.exists(OUTDIR)) {
    stop(paste("Output directory already exists; refusing to overwrite:", OUTDIR))
}

count_table <- read.csv(
    COUNTS_FILE,
    check.names = FALSE,
    stringsAsFactors = FALSE
)
if (ncol(count_table) < 2L) {
    stop("Raw-count matrix must contain a gene-ID column and sample columns")
}
gene_ids <- count_table[[1]]
count_data <- count_table[-1]
original_sample_ids <- colnames(count_data)
expected_sample_ids <- c(
    "SRR3734796",
    "SRR3734797",
    "SRR3734798",
    "SRR3734820",
    "SRR3734824",
    "SRR3734825"
)

if (nrow(count_data) != 14214L) {
    stop(paste("Expected 14214 gene records; found", nrow(count_data)))
}
if (!is.character(gene_ids) || anyNA(gene_ids) ||
    any(!nzchar(gene_ids)) || anyDuplicated(gene_ids)) {
    stop("Gene IDs must be present, non-empty, and unique")
}
if (ncol(count_data) != 6L) {
    stop(paste("Expected six count columns; found", ncol(count_data)))
}
if (anyNA(original_sample_ids) || any(!nzchar(original_sample_ids)) ||
    anyDuplicated(original_sample_ids)) {
    stop("Count-matrix sample identifiers must be present, non-empty, and unique")
}
if (!all(vapply(count_data, is.numeric, logical(1)))) {
    stop("All count-matrix sample columns must contain numeric values")
}

# Map the original BAM-based column names explicitly to SRA accessions.
sample_name_pattern <- "^SRR[0-9]+_Aligned\\.sortedByCoord\\.out\\.bam$"
if (any(!grepl(sample_name_pattern, original_sample_ids))) {
    stop("Count-matrix sample names do not match the expected STAR BAM basename pattern")
}
mapped_sample_ids <- sub(
    "_Aligned\\.sortedByCoord\\.out\\.bam$",
    "",
    original_sample_ids
)
if (anyDuplicated(mapped_sample_ids) ||
    !setequal(mapped_sample_ids, expected_sample_ids)) {
    stop("Count-matrix columns do not map exactly to the six expected SRA accessions")
}

count_matrix <- as.matrix(count_data)
if (!is.numeric(count_matrix) || anyNA(count_matrix) ||
    any(!is.finite(count_matrix))) {
    stop("Counts must be numeric, finite, and contain no missing values")
}
if (any(count_matrix < 0) || any(count_matrix != floor(count_matrix))) {
    stop("Counts must be non-negative integers; values will not be rounded or truncated")
}
if (any(count_matrix > .Machine$integer.max)) {
    stop("Counts exceed the maximum safely representable R integer")
}
storage.mode(count_matrix) <- "integer"
rownames(count_matrix) <- gene_ids
colnames(count_matrix) <- original_sample_ids

sampleInfo <- read.csv(
    PHENO,
    check.names = FALSE,
    stringsAsFactors = FALSE
)
if (!all(c("sampleID", "condition") %in% names(sampleInfo))) {
    stop("Metadata must contain sampleID and condition columns")
}
sampleInfo$sampleID <- as.character(sampleInfo$sampleID)
sampleInfo$condition <- as.character(sampleInfo$condition)
if (nrow(sampleInfo) != length(expected_sample_ids) ||
    anyNA(sampleInfo$sampleID) || any(!nzchar(sampleInfo$sampleID)) ||
    anyDuplicated(sampleInfo$sampleID)) {
    stop("Metadata must contain exactly six unique, non-empty sample IDs")
}
if (any(!grepl(sample_name_pattern, sampleInfo$sampleID))) {
    stop("Metadata sample names do not match the expected STAR BAM basename pattern")
}
metadata_sra_ids <- sub(
    "_Aligned\\.sortedByCoord\\.out\\.bam$",
    "",
    sampleInfo$sampleID
)
if (anyDuplicated(metadata_sra_ids) ||
    !setequal(metadata_sra_ids, expected_sample_ids)) {
    stop("Metadata sample names do not map exactly to the six expected SRA accessions")
}
if (anyNA(sampleInfo$condition) ||
    !all(sampleInfo$condition %in% c("1mo", "4mo"))) {
    stop("Metadata contains missing or unexpected conditions")
}
condition_counts <- table(factor(sampleInfo$condition, levels = c("1mo", "4mo")))
if (!all(condition_counts == 3L)) {
    stop("Expected exactly three samples at 1mo and three samples at 4mo")
}

metadata_index <- match(mapped_sample_ids, metadata_sra_ids)
if (anyNA(metadata_index)) {
    stop("Could not match every count-matrix sample to metadata")
}
sampleInfo <- sampleInfo[metadata_index, , drop = FALSE]
sampleInfo$sra_accession <- metadata_sra_ids[metadata_index]
sampleInfo$matrix_sample_id <- original_sample_ids

# Verify sample matching after metadata reordering
stopifnot(identical(sampleInfo$sra_accession, mapped_sample_ids))
stopifnot(identical(sampleInfo$sampleID, original_sample_ids))

sampleInfo$condition <- factor(sampleInfo$condition, levels = c("1mo", "4mo"))
rownames(sampleInfo) <- original_sample_ids
if (!identical(colnames(count_matrix), rownames(sampleInfo))) {
    stop("Metadata could not be aligned exactly to the count-matrix columns")
}

count_md5 <- unname(tools::md5sum(COUNTS_FILE))
metadata_md5 <- unname(tools::md5sum(PHENO))
sample_records <- paste(
    paste(
        original_sample_ids,
        sampleInfo$sampleID,
        sampleInfo$sra_accession,
        as.character(sampleInfo$condition),
        sep = " | "
    ),
    collapse = "\n"
)
provenance_lines <- c(
    paste("Input count matrix:", COUNTS_FILE),
    paste("Input count matrix MD5:", count_md5),
    paste("Metadata:", PHENO),
    paste("Metadata MD5:", metadata_md5),
    paste("Number of genes:", nrow(count_matrix)),
    paste("Number of samples:", ncol(count_matrix)),
    "Samples (matrix identifier | metadata identifier | SRA accession | condition):",
    sample_records,
    "Statistical model: ~ condition",
    "Reference group: 1mo",
    'Statistical contrast: c("condition", "4mo", "1mo")',
    paste("Execution date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
)

# Create new output directories only after all inputs have passed validation.
if (!dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)) {
    stop(paste("Could not create output directory; it may already exist:", OUTDIR))
}
if (!dir.exists(OUTDIR)) {
    stop(paste("Output directory was not created:", OUTDIR))
}
if (!dir.create(RESULTSDIR, showWarnings = FALSE) || !dir.exists(RESULTSDIR)) {
    stop(paste("Could not create statistical results directory:", RESULTSDIR))
}

output_files <- c(
    file.path(RESULTSDIR, "raw_counts.1mo_4mo.psiclass.csv"),
    file.path(RESULTSDIR, "normalized_counts.1mo_4mo.psiclass.csv"),
    file.path(RESULTSDIR, "results_1mo_4mo.psiclass.csv"),
    file.path(RESULTSDIR, "results_1mo_4mo.pval05.psiclass.csv"),
    file.path(RESULTSDIR, "results_1mo_4mo.qval05.psiclass.csv"),
    file.path(OUTDIR, "PCA_1mo_4mo_byCondition.png"),
    file.path(OUTDIR, "1mo_4mo_pca_table.psiclass.csv"),
    file.path(OUTDIR, "volcano_1mo_4mo.psiclass.png"),
    file.path(OUTDIR, "analysis_provenance.txt"),
    file.path(OUTDIR, "sessionInfo.txt"),
    file.path(OUTDIR, "analysis_summary.txt")
)
if (any(file.exists(output_files) | dir.exists(output_files))) {
    stop("An expected output path already exists; refusing to overwrite")
}

writeLines(provenance_lines, file.path(OUTDIR, "analysis_provenance.txt"))

# These are raw counts from the sex-adjusted workflow, not sex-adjusted values.
write.csv(
    count_matrix,
    file.path(RESULTSDIR, "raw_counts.1mo_4mo.psiclass.csv")
)

# Build the unadjusted DESeq2 dataset from the validated shared count matrix.
dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = sampleInfo,
    design = ~condition
)
dds$condition <- relevel(dds$condition, ref = "1mo")
dds <- DESeq2::DESeq(dds)

# Comparison: 4mo vs 1mo
res <- DESeq2::results(dds, contrast = c("condition", "4mo", "1mo"))
res_df <- as.data.frame(res)

# Export normalized counts and complete DESeq2 results.
write.csv(
    DESeq2::counts(dds, normalized = TRUE),
    file.path(RESULTSDIR, "normalized_counts.1mo_4mo.psiclass.csv")
)
write.csv(
    res_df,
    file.path(RESULTSDIR, "results_1mo_4mo.psiclass.csv")
)

# Filter by unadjusted p-value and adjusted p-value.
res_pval <- res_df[
    !is.na(res_df$pvalue) & res_df$pvalue <= 0.05, ,
    drop = FALSE
]
res_qval <- res_df[
    !is.na(res_df$padj) & res_df$padj <= 0.05, ,
    drop = FALSE
]
write.csv(
    res_pval,
    file.path(RESULTSDIR, "results_1mo_4mo.pval05.psiclass.csv")
)
write.csv(
    res_qval,
    file.path(RESULTSDIR, "results_1mo_4mo.qval05.psiclass.csv")
)

# Variance-stabilizing transformation and PCA visualization.
vsd <- DESeq2::vst(dds, blind = TRUE)
pca_plot <- DESeq2::plotPCA(vsd, intgroup = "condition")
ggplot2::ggsave(
    filename = file.path(OUTDIR, "PCA_1mo_4mo_byCondition.png"),
    plot = pca_plot,
    width = 7,
    height = 5,
    dpi = 300
)

# Export PCA coordinates with both original and SRA sample identifiers.
pca <- prcomp(t(SummarizedExperiment::assay(vsd)))
pca_table <- data.frame(
    matrix_sample_id = original_sample_ids,
    sampleID = sampleInfo$sampleID,
    sra_accession = sampleInfo$sra_accession,
    condition = as.character(dds$condition),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    row.names = NULL
)
write.csv(
    pca_table,
    file.path(OUTDIR, "1mo_4mo_pca_table.psiclass.csv"),
    row.names = FALSE
)

# Volcano plot
png(
    filename = file.path(OUTDIR, "volcano_1mo_4mo.psiclass.png"),
    width = 1400,
    height = 1100,
    res = 150
)
print(
    EnhancedVolcano::EnhancedVolcano(
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

# Save session details and summarize the completed analysis.
writeLines(capture.output(sessionInfo()), file.path(OUTDIR, "sessionInfo.txt"))
summary_lines <- c(
    "DESeq2 analysis completed.",
    paste("Total genes:", nrow(res_df)),
    paste("Genes with p-value <= 0.05:", nrow(res_pval)),
    paste("Genes with padj <= 0.05:", nrow(res_qval)),
    paste("Model: ~ condition (reference: 1mo; contrast: 4mo vs 1mo)"),
    paste("Results saved to:", RESULTSDIR)
)
writeLines(summary_lines, file.path(OUTDIR, "analysis_summary.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")
