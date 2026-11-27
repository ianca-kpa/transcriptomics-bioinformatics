# Module I — Sex-adjusted DESeq2 and method comparison

*Dataset: GSE83931 — mouse hippocampus RNA-seq*

*Comparison: 4 months (4mo) vs 1 month (1mo)*

This note compares the genome-wide unadjusted and sex-adjusted DESeq2 results,
and the sex-adjusted results with Cuffdiff2. I keep these genome-wide
comparisons separate from the regional Proj1 analysis described in the
RNA-seq notebook. That regional workflow covers approximately 3.4 Mb of
chromosome 7 and contains 142 loci; these loci are distinct from the 14,214
gene IDs in the genome-wide analyses.

## Step-by-step workflow

1. Identify which files are current unadjusted results, previously verified sex-adjusted results, and the not-yet-run destination for corrected script 13b.
2. Confirm the sample distribution, models, reference groups, contrast, input data, and significance criteria.
3. Reproduce the unadjusted-versus-adjusted DESeq2 comparison by gene ID using the code in Section 5.
4. Reproduce the Cuffdiff2-versus-sex-adjusted DESeq2 comparison by gene ID using the code in Section 6.
5. Check shared and exclusive gene IDs and compare log2FC directions before interpreting the results.
6. Discuss significance and effect-size differences cautiously, given the unbalanced design and small sample size.
7. Treat rerunning corrected script 13b as optional; do not create its protected output directory in advance.

## 1. Purpose and analysis status

Script 13 completed successfully. Its unadjusted results are available under
`$BASEDIR/DESeq2_genomewide/unadjusted/`. The corrected
`scripts/13b_run_deseq2_with_sex.R` is prepared to write to
`$BASEDIR/DESeq2_genomewide/sex_adjusted/`, but it has not been rerun. The
previously verified sex-adjusted results remain under the legacy
`$BASEDIR/results-proj1/DESeq2/Results/` path.

The corrected script and the previous adjusted execution are distinct. I use
the previous results for the exploratory comparisons below, and I do not
present them as outputs of the corrected script.

## 2. Experimental design and sex imbalance

The six samples have this distribution:

| Condition | Female | Male |
| --- | ---: | ---: |
| 1mo | 3 | 0 |
| 4mo | 1 | 2 |

The 1mo group has no male samples, and the total sample size is small. This
imbalance limits how broadly I can interpret age and sex effects.

## 3. Statistical models and inputs

The unadjusted DESeq2 analysis uses the existing raw-count matrix and the
model `~ condition`. It does not recount BAMs.

The corrected script 13b instead reads the six genome-wide coordinate-sorted
STAR BAMs under `$BASEDIR/STAR/BAM/`, imports the genome-wide PsiCLASS
annotation at
`$BASEDIR/PsiCLASS/WithGeneNames/psiclass_vote.withStrand.gtf`, and uses the
validated metadata at
`$BASEDIR/results-proj1/DESeq2/pheno_1mo_4mo.w-sex.csv`. It builds a
transcript database, groups exons by gene, and runs `summarizeOverlaps()` to
generate raw gene counts before fitting the model:

```r
se <- summarizeOverlaps(
    features = exonsByGene,
    reads = bamLst,
    mode = "Union",
    singleEnd = TRUE,
    ignore.strand = TRUE,
    fragments = FALSE
)
dds <- DESeqDataSet(se, design = ~ sex + condition)
```

The count matrix itself is not statistically adjusted for sex; sex enters
the DESeq2 model, not the count values.

The model specifications are:

```r
# Unadjusted
design = ~ condition

# Sex-adjusted
design = ~ sex + condition
```

Both use `1mo` as the condition reference and test `4mo` versus `1mo`. In
both result tables, a positive `log2FoldChange` indicates higher estimated
expression at 4mo. The adjusted model uses female as the sex reference.

The versioned course metadata example is
`course_materials/Module_I/metadata/pheno_1mo_4mo.A.csv`. It is a course
reference whose sample IDs end in `.so.bam`; it is not the exact metadata
file consumed by corrected script 13b, which matches full STAR BAM basenames
such as `SRR3734796_Aligned.sortedByCoord.out.bam`.

## 4. Provenance and output locations

| Analysis or input | Location and status |
| --- | --- |
| Raw-count matrix used by script 13 | `$BASEDIR/results-proj1/DESeq2/Results/raw_counts.1mo_4mo.A.psiclass.csv` |
| Current unadjusted outputs | `$BASEDIR/DESeq2_genomewide/unadjusted/` — completed |
| Previously verified sex-adjusted results | `$BASEDIR/results-proj1/DESeq2/Results/` |
| Corrected script 13b destination | `$BASEDIR/DESeq2_genomewide/sex_adjusted/` — not rerun |
| Cuffdiff2 gene results | `$BASEDIR/Cuffdiff/gene_exp.diff` |

If I run the corrected script 13b in the future, its adjusted result table
will be written to
`$BASEDIR/DESeq2_genomewide/sex_adjusted/Results/results_1mo_4mo.A.psiclass.csv`.
That file is not present yet; the previously verified adjusted table remains
at the legacy path above.

The future run will also export raw and normalized counts, the p-value and
adjusted-p-value filtered tables, a PCA plot and table, and a volcano plot
under `$BASEDIR/DESeq2_genomewide/sex_adjusted/`. These are expected outputs
of a future execution, not files produced by the corrected script so far.
After a successful run, script 13b will also write `analysis_provenance.txt`,
`analysis_summary.txt`, and `sessionInfo.txt` in the output directory. The
provenance will record the script and metadata paths and checksums, the GTF
path and checksum, the generated raw-count matrix path and checksum, the
model, contrast, sample metadata, execution summary, and output paths. It
will not calculate checksums for the BAM files.

The `.A` raw-count matrix and the raw-count export from the current
unadjusted run have identical gene IDs, sample identifiers, sample order and
count values. The previous adjusted result table and execution log report
the same 14,214-gene universe and six samples. This supports comparing the
two DESeq2 result tables using the shared raw counts.

The previous adjusted execution has incomplete environment provenance. Its
`.Rout` records the results and non-fatal warnings, but does not provide a
complete `sessionInfo()` or establish the exact R and DESeq2 versions used.
The new provenance files will describe a future successful execution; they
do not establish the environment of the previous run. The current unadjusted
run records R 4.5.2 and DESeq2 1.50.2; I cannot assume the previous adjusted
run used the same versions. The available Cuffdiff2 `run.info` records
version 2.2.1.

## 5. Unadjusted versus sex-adjusted DESeq2

Both DESeq2 result tables contain the same 14,214 unique gene IDs and use
the contrast `4mo` versus `1mo`. I define significance as `padj <= 0.05`.

| Result | Unadjusted (`~ condition`) | Sex-adjusted (`~ sex + condition`) |
| --- | ---: | ---: |
| Significant genes | 249 | 124 |
| Positive `log2FoldChange` among significant genes | 67 | 22 |
| Negative `log2FoldChange` among significant genes | 182 | 102 |

| Gene-set comparison | Count |
| --- | ---: |
| Significant in both models | 85 |
| Significant only in the unadjusted model | 164 |
| Significant only in the sex-adjusted model | 39 |
| Shared significant genes with concordant effect direction | 85 of 85 |

I reproduce the comparison from the result tables by matching gene IDs, not
by subtracting aggregate counts:

From the repository root, I load the local path configuration and export
`BASEDIR` so it is available to child R processes. Both R examples below
check that the variable is set and points to an existing directory.

```bash
source configs/paths.sh
export BASEDIR
```

```r
BASEDIR <- Sys.getenv("BASEDIR")
if (!nzchar(BASEDIR) || !dir.exists(BASEDIR)) {
    stop("Set BASEDIR to an existing data directory before running this example.")
}
BASEDIR <- normalizePath(BASEDIR, winslash = "/", mustWork = TRUE)

unadjusted_file <- file.path(
    BASEDIR, "DESeq2_genomewide", "unadjusted", "Results",
    "results_1mo_4mo.psiclass.csv"
)
adjusted_file <- file.path(
    BASEDIR, "results-proj1", "DESeq2", "Results",
    "results_1mo_4mo.A.psiclass.csv"
)
unadjusted <- read.csv(unadjusted_file, row.names = 1, check.names = FALSE)
adjusted <- read.csv(adjusted_file, row.names = 1, check.names = FALSE)

unadjusted_ids <- rownames(unadjusted)
adjusted_ids <- rownames(adjusted)
stopifnot(
    !anyDuplicated(unadjusted_ids), !anyNA(unadjusted_ids),
    all(nzchar(unadjusted_ids)), length(unadjusted_ids) == 14214L
)
stopifnot(
    !anyDuplicated(adjusted_ids), !anyNA(adjusted_ids),
    all(nzchar(adjusted_ids)), length(adjusted_ids) == 14214L
)
stopifnot(setequal(unadjusted_ids, adjusted_ids))

unadjusted_sig_ids <- unadjusted_ids[
    !is.na(unadjusted$padj) & unadjusted$padj <= 0.05
]
adjusted_sig_ids <- adjusted_ids[
    !is.na(adjusted$padj) & adjusted$padj <= 0.05
]
shared_ids <- intersect(unadjusted_sig_ids, adjusted_sig_ids)
unadjusted_only_ids <- setdiff(unadjusted_sig_ids, adjusted_sig_ids)
adjusted_only_ids <- setdiff(adjusted_sig_ids, unadjusted_sig_ids)

unadjusted_lfc <- setNames(unadjusted$log2FoldChange, unadjusted_ids)
adjusted_lfc <- setNames(adjusted$log2FoldChange, adjusted_ids)
direction_agrees <- sign(unadjusted_lfc[shared_ids]) ==
    sign(adjusted_lfc[shared_ids])
stopifnot(!anyNA(direction_agrees))

cat("Unadjusted significant:", length(unadjusted_sig_ids), "\n")
cat("Sex-adjusted significant:", length(adjusted_sig_ids), "\n")
cat("Shared:", length(shared_ids), "\n")
cat("Unadjusted only:", length(unadjusted_only_ids), "\n")
cat("Sex-adjusted only:", length(adjusted_only_ids), "\n")
cat("Concordant directions:", sum(direction_agrees), "/", length(shared_ids), "\n")

stopifnot(
    length(unadjusted_sig_ids) == 249L,
    length(adjusted_sig_ids) == 124L,
    length(shared_ids) == 85L,
    length(unadjusted_only_ids) == 164L,
    length(adjusted_only_ids) == 39L,
    sum(direction_agrees) == 85L
)
```

## 6. Cuffdiff2 versus sex-adjusted DESeq2

The Cuffdiff2 result file and adjusted DESeq2 table each contain 14,214
unique gene IDs. I define Cuffdiff2 significance using
`status == "OK"` and `significant == "yes"`, and DESeq2 significance using
`padj <= 0.05`. Cuffdiff2 used 1mo as the control group and 4mo as the test
group, so I compare its `log2(fold_change)` in the same 4mo-versus-1mo
direction as the DESeq2 `log2FoldChange`.

| Result | Cuffdiff2 | Sex-adjusted DESeq2 |
| --- | ---: | ---: |
| Significant genes | 261 | 124 |

| Gene-set comparison | Count |
| --- | ---: |
| Significant in both methods | 76 |
| Significant only in Cuffdiff2 | 185 |
| Significant only in sex-adjusted DESeq2 | 48 |
| Shared significant genes with concordant effect direction | 76 of 76 |

This code reproduces the comparison from the gene-level result rows:

```r
BASEDIR <- Sys.getenv("BASEDIR")
if (!nzchar(BASEDIR) || !dir.exists(BASEDIR)) {
    stop("Set BASEDIR to an existing data directory before running this example.")
}
BASEDIR <- normalizePath(BASEDIR, winslash = "/", mustWork = TRUE)

adjusted_file <- file.path(
    BASEDIR, "results-proj1", "DESeq2", "Results",
    "results_1mo_4mo.A.psiclass.csv"
)
adjusted <- read.csv(adjusted_file, row.names = 1, check.names = FALSE)
cuff <- read.delim(
    file.path(BASEDIR, "Cuffdiff", "gene_exp.diff"),
    check.names = FALSE,
    stringsAsFactors = FALSE
)

adjusted_ids <- rownames(adjusted)
stopifnot(
    !anyDuplicated(adjusted_ids), !anyNA(adjusted_ids),
    all(nzchar(adjusted_ids)), length(adjusted_ids) == 14214L
)
stopifnot(
    !anyDuplicated(cuff$gene_id), !anyNA(cuff$gene_id),
    all(nzchar(cuff$gene_id)), length(cuff$gene_id) == 14214L
)
stopifnot(setequal(adjusted_ids, cuff$gene_id))

adjusted_sig_ids <- adjusted_ids[
    !is.na(adjusted$padj) & adjusted$padj <= 0.05
]
cuff_sig_rows <- !is.na(cuff$status) & cuff$status == "OK" &
    !is.na(cuff$significant) & tolower(cuff$significant) == "yes"
cuff_sig_ids <- cuff$gene_id[cuff_sig_rows]
shared_ids <- intersect(cuff_sig_ids, adjusted_sig_ids)
cuff_only_ids <- setdiff(cuff_sig_ids, adjusted_sig_ids)
adjusted_only_ids <- setdiff(adjusted_sig_ids, cuff_sig_ids)

cuff_lfc <- setNames(cuff$`log2(fold_change)`, cuff$gene_id)
adjusted_lfc <- setNames(adjusted$log2FoldChange, adjusted_ids)
direction_agrees <- sign(cuff_lfc[shared_ids]) ==
    sign(adjusted_lfc[shared_ids])
stopifnot(!anyNA(direction_agrees))

cat("Cuffdiff2 significant:", length(cuff_sig_ids), "\n")
cat("Sex-adjusted DESeq2 significant:", length(adjusted_sig_ids), "\n")
cat("Shared:", length(shared_ids), "\n")
cat("Cuffdiff2 only:", length(cuff_only_ids), "\n")
cat("Sex-adjusted DESeq2 only:", length(adjusted_only_ids), "\n")
cat("Concordant directions:", sum(direction_agrees), "/", length(shared_ids), "\n")

stopifnot(
    length(cuff_sig_ids) == 261L,
    length(adjusted_sig_ids) == 124L,
    length(shared_ids) == 76L,
    length(cuff_only_ids) == 185L,
    length(adjusted_only_ids) == 48L,
    sum(direction_agrees) == 76L
)
```

Cuffdiff2 compares the 1mo and 4mo groups without the sex covariate used by
the adjusted DESeq2 model. The two methods also use different statistical
procedures, so their p-values and significance totals are not interchangeable.
The overlap and direction agreement are descriptive comparisons, not proof
of biological validity.

## 7. Interpretation and limitations

I distinguish changes in estimated effect size from changes in statistical
significance. Adding a covariate can change model estimates and their
uncertainty; a gene crossing a significance threshold does not, by itself,
demonstrate a sex effect. Concordant directions among shared genes do not
mean that the effect magnitudes are equal.

The adjusted model assumes an additive sex effect and does not estimate a
sex-by-condition interaction. Since there are no male 1mo samples, the
design cannot establish how the age contrast behaves in males at 1mo.
Together with the small sample size, this limits interpretation of both
models. I do not treat these comparisons as causal evidence.

## 8. Reproducibility and next steps

The previous sex-adjusted result table and its `.Rout` support the reported
comparison, but the historical execution environment is not fully recorded.
A future run of corrected script 13b could improve reproducibility; it is
optional and is not required for the exploratory comparisons documented here.

For a future rerun, after confirming the `transcriptomics` Conda environment
and required R packages, I would run the following from the repository root:

```bash
conda activate transcriptomics
Rscript scripts/13b_run_deseq2_with_sex.R
```

The corrected script resolves `BASEDIR`, validates its inputs, and creates
its output directory itself. I must not create
`$BASEDIR/DESeq2_genomewide/sex_adjusted/` in advance, because the script
refuses an existing destination. A successful run will save
`analysis_provenance.txt`, `analysis_summary.txt`, and `sessionInfo.txt`
alongside the results. I can capture the console log separately, keeping it
out of the legacy Proj1 directories.

*Course exercise: Johns Hopkins University / Coursera. The comparisons and paths in this note describe my local genome-wide workflows.*
