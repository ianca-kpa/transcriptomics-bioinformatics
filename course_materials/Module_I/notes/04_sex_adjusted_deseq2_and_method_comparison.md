# Module I — Sex-adjusted DESeq2 and method comparison

*Dataset: GSE83931 — mouse hippocampus RNA-seq*

*Comparison: 4 months (4mo) vs 1 month (1mo)*

This note follows the original DESeq2 analysis in [Note 3](03_deseq2_differential_expression.md). It records how to validate the sex-adjusted model and compare its results with the unadjusted DESeq2 and Cuffdiff2 analyses.

**Status:** the sex-adjusted output files were not present under the currently configured `$BASEDIR/results-proj1/DESeq2/` when this note was prepared. No sex-adjusted result counts or conclusions are reported here.

## Why fit a sex-adjusted model?

The six samples have this distribution:

| Condition | Female | Male |
|---|---:|---:|
| 1mo | 3 | 0 |
| 4mo | 1 | 2 |

The original model is `~ condition`. The follow-up script, `scripts/13b_run_deseq2_with_sex.R`, uses the additive model `~ sex + condition`, with female and 1mo as reference levels. Its reported contrast is 4mo vs 1mo, so a positive `log2FoldChange` indicates higher estimated expression at 4mo.

This model does not include a sex-by-condition interaction. Since there are no male 1mo samples, the design is unbalanced; adding sex as a covariate does not make the groups balanced or eliminate all uncertainty about age and sex effects. Interpret the adjusted estimates cautiously.

## Inputs and output locations

The script expects:

- Six queryname-sorted BAMs under `$BASEDIR/results-proj1/Sorted/`, with names ending in `.so.bam`.
- The annotated PsiCLASS GTF at `$BASEDIR/PsiCLASS/WithGeneNames/psiclass_vote.withStrand.gtf`.
- Metadata at `$BASEDIR/results-proj1/DESeq2/pheno_1mo_4mo.w-sex.csv`, with columns `sampleID`, `condition`, and `sex`.

The course metadata example is `course_materials/Module_I/docs/pheno_1mo_4mo.A.csv`. The script requires the local copy at the configured output location; verify the sample IDs match the BAM filenames exactly.

The script uses single-end counting (`singleEnd = TRUE`, `fragments = FALSE`) and writes its count matrix and DESeq2 tables under:

```text
$BASEDIR/results-proj1/DESeq2/Results/
```

Expected result files include:

- `raw_counts.1mo_4mo.A.psiclass.csv`
- `normalized_counts.1mo_4mo.A.psiclass.csv`
- `results_1mo_4mo.A.psiclass.csv`
- `results_1mo_4mo.pval05.A.psiclass.csv`
- `results_1mo_4mo.qval05.A.psiclass.csv`

The script also writes PCA and volcano plots under `$BASEDIR/results-proj1/DESeq2/`.

## Running and checking the analysis

From the repository root, activate the `transcriptomics` environment and confirm `configs/paths.sh` points to the intended data directory. To retain the run log and preserve the R process exit status:

```bash
source configs/paths.sh
mkdir -p "$DESEQ2_DIR"
set -o pipefail
Rscript scripts/13b_run_deseq2_with_sex.R \
    2>&1 | tee "${DESEQ2_DIR}/deseq2_with_sex.log"
```

Before interpreting results, confirm that:

1. The metadata has six unique sample IDs matching the sorted BAM filenames.
2. The condition/sex table is `1mo: F=3, M=0` and `4mo: F=1, M=2`.
3. The script completes successfully and all expected output files are nonempty.
4. The number of rows and sample columns in the raw count matrix are valid, and its sample order agrees with the metadata.
5. The result table has the expected columns (`baseMean`, `log2FoldChange`, `pvalue`, and `padj`) and no unexpected missing values in the sample or condition metadata.

Do not treat the saved results from the unadjusted model as the output of this script. The `.A.psiclass.csv` suffix identifies the sex-adjusted analysis.

## Comparing significant gene IDs

Compare gene identifiers and effect directions, not only the number of significant rows. For the adjusted DESeq2 model, use `padj <= 0.05`; for Cuffdiff2, restrict to `status == "OK"` and its `significant == "yes"` flag. The result files are produced by different statistical methods, so their p-values and effect estimates are not interchangeable.

Example comparison of significant IDs:

```r
# Set this to the BASEDIR value in configs/paths.sh.
BASEDIR <- "/path/to/transcriptomics-data"
CUFFDIFF_DIR <- file.path(BASEDIR, "Cuffdiff")
DESEQ2_RESULTS_DIR <- file.path(
    BASEDIR,
    "results-proj1/DESeq2/Results"
)

cuff <- read.delim(
    file.path(CUFFDIFF_DIR, "gene_exp.diff"),
    stringsAsFactors = FALSE
)

deseq <- read.csv(
    file.path(
        DESEQ2_RESULTS_DIR,
        "results_1mo_4mo.A.psiclass.csv"
    ),
    row.names = 1,
    check.names = FALSE
)

cuff_sig <- unique(cuff$gene_id[
    cuff$status == "OK" &
    tolower(cuff$significant) == "yes"
])

deseq_sig <- rownames(deseq)[
    !is.na(deseq$padj) & deseq$padj <= 0.05
]

length(intersect(cuff_sig, deseq_sig))
setdiff(cuff_sig, deseq_sig)
setdiff(deseq_sig, cuff_sig)

shared_sig <- intersect(cuff_sig, deseq_sig)
if (length(shared_sig) > 0) {
    cuff_effect <- cuff$`log2(fold_change)`[
        match(shared_sig, cuff$gene_id)
    ]
    deseq_effect <- deseq$log2FoldChange[shared_sig]
    mean(sign(cuff_effect) == sign(deseq_effect), na.rm = TRUE)
} else {
    NA_real_
}
```

The same ID-based comparison can be made between the original and sex-adjusted DESeq2 result tables. Check that the input samples, annotation, and significance criteria are documented for each comparison. Record the overlap, model-specific genes, and direction agreement only after verifying the result files.

## Interpretation limits

- The sex-adjusted model assumes an additive sex effect and does not estimate an interaction.
- There are no male 1mo samples, and the total sample size is small.
- Statistical adjustment does not establish that a gene change is caused by age or sex.
- Cuffdiff2 and DESeq2 differ in their statistical models; compare gene IDs and directions rather than expecting identical significance totals.

*Add numerical results and conclusions only after the sex-adjusted script has been run and its outputs verified.*
