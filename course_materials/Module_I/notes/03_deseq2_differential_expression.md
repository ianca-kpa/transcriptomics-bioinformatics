# Module I — DESeq2 differential expression analysis

**Course:** Bioinformatics Methods for Transcriptomics — Johns Hopkins University (Coursera)  
**Dataset:** GSE83931 — mouse hippocampus RNA-seq  
**Comparison:** 4 months (4mo) vs 1 month (1mo)  
**Recorded working environment:** Ubuntu / WSL2, Conda (`transcriptomics`), R 4.5.2

I continued with the same six RNA-seq samples used for STAR, PsiCLASS and Cuffdiff2. For this part of the course, I followed the DESeq2 exercise and changed the paths and a few settings to run it locally.

**Workflow and result provenance:** I keep the regional Proj1 analysis separate from the current genome-wide analysis. Proj1 contains 142 loci and was performed through the RNA-seq notebook; the current unadjusted genome-wide run contains 14,214 genes and 249 genes with `padj <= 0.05`. These workflows have different gene universes, so I do not combine their results.

Script 13 has completed successfully; its results are available under
`$BASEDIR/DESeq2_genomewide/unadjusted/`. It reads a pre-existing genome-wide raw-count
matrix and does not recount BAMs. Script 13b has been corrected to target
`$BASEDIR/DESeq2_genomewide/sex_adjusted/`, but has not been rerun, so that corrected
script has not produced results there. The previously verified sex-adjusted results
remain under the legacy `$BASEDIR/results-proj1/DESeq2/Results/` path. Rerunning script
13b is optional for improving reproducibility, not a prerequisite for the existing
exploratory comparison.

## Step-by-step workflow for the current unadjusted analysis

1. Validate the existing raw-count matrix and metadata: 14,214 unique gene IDs, six expected samples, integer raw counts, and matching sample order.
2. Set `~ condition`, use `1mo` as the reference, and define the contrast as `4mo` versus `1mo`.
3. The completed run used `Rscript scripts/13_run_deseq2.R` from the repository root. Script 13 imported the matrix; it did not recount BAMs or run `summarizeOverlaps()`. Its output directory now exists, so the overwrite safeguard prevents simply rerunning it there.
4. Review the complete results, filtered tables, normalized counts, PCA, volcano plot, and analysis summary in `$BASEDIR/DESeq2_genomewide/unadjusted/`.
5. Confirm the exported raw counts match the input, then interpret the genome-wide results separately from the regional Proj1 analysis.

The six samples are:

| Condition | SRA runs | Biological sex |
| --- | --- | --- |
| 1mo | SRR3734796, SRR3734797, SRR3734798 | Female, female, female |
| 4mo | SRR3734820, SRR3734824, SRR3734825 | Female, male, male |

**Dataset limitation:** the 1mo samples are all female, whereas the 4mo group contains one female and two males. My first DESeq2 model does not account for sex. Adding sex to the model later does not fix the lack of male 1mo samples.

## 1. Preparing BAM files for gene counting

This section records the earlier genome-wide BAM-sorting exercise in scripts 11 and 12.
Those scripts read coordinate-sorted STAR BAMs from `STAR_BAM_DIR` and write
queryname-sorted BAMs to `SORTED_BAM_DIR`, as defined in `configs/paths.sh`. In the
current configuration, these resolve to `$BASEDIR/STAR/BAM/` and
`$BASEDIR/results-proj1/Sorted/`, respectively. The BAMs currently stored in
`results-proj1/Sorted/` are regional Proj1 data; their presence there does not establish
that they are the outputs of the earlier genome-wide sorting exercise. Since the
configured destination is now used for regional BAMs, I would not rerun script 11 with
this configuration without first choosing a separate output directory. Script 13 does
not use BAMs directly; it reads the existing count matrix. Script 13b uses the original
coordinate-sorted STAR BAMs.

STAR produced BAM files sorted by genomic coordinate, for example:

```text
$BASEDIR/STAR/BAM/SRR3734796_Aligned.sortedByCoord.out.bam
```

The historical destination for that genome-wide sorting exercise is not confirmed. This is an illustrative syntax example only; do not run it with the current `SORTED_BAM_DIR`, which points to the regional `results-proj1/Sorted/` directory:

```bash
/usr/bin/samtools sort -n -@ 2 -O BAM \
    -o "${GENOMEWIDE_SORTED_DIR}/SRR3734796.so.bam" \
    "${STAR_BAM_DIR}/SRR3734796_Aligned.sortedByCoord.out.bam"
```

- `sort` rearranges alignment **records**; it does not realign reads or change their genomic coordinates.
- Without `-n`, `samtools sort` sorts by genomic coordinates (reference and position).
- `-n` sorts by read name (QNAME). Alignments belonging to the same read are grouped together.
- `-@ 2` requests two additional processing threads.
- `-O BAM` specifies BAM output; `-o` gives the destination file.

A useful genome-wide multimapping example was **`SRR3734796.sra.8`**, which had two reported alignments:

| Read | Genomic position | Note |
| --- | --- | --- |
| `SRR3734796.sra.8` | chr1:95426560 | Primary alignment, FLAG 0 |
| `SRR3734796.sra.8` | chr8:71348574 | Secondary reverse-strand alignment, FLAG 272 |

Sorting by QNAME keeps these records together. Coordinate sorting places them according
to their genomic locations. The read was reported as mapping to two locations
(`NH:i:2`); this illustrates **multimapping**. Because the alignments are on chr1 and
chr8, this example cannot come from the current chr7-only regional BAMs in
`results-proj1/Sorted/`.

Scripts:

- `scripts/11_sort_bam_by_read_name.sh`
- `scripts/12_check_sorted_bam.sh`

Script 12 checks the six queryname-sorted files in the configured `SORTED_BAM_DIR` using
`samtools quickcheck` and verifies that their headers declare `SO:queryname`. These
checks establish basic file integrity and sort-order metadata, but do not verify every
record, genomic scope or the provenance of files currently in that directory. Neither
the regional BAMs there now nor these sorting scripts are inputs to the current script
13.

## 2. From alignments to a gene-count matrix

DESeq2 uses a matrix of **raw gene counts** rather than the FPKM-style expression estimates produced by Cuffdiff2.

The gene annotation came from the genome-wide PsiCLASS step:

```text
$BASEDIR/PsiCLASS/WithGeneNames/psiclass_vote.withStrand.gtf
```

The raw-count matrix used by the current unadjusted analysis was generated from the genome-wide BAMs during an earlier execution of the sex-adjusted workflow. Counting happened before DESeq2 fitted the sex-adjusted model: the R script constructed a transcript database, grouped exons by gene and counted alignments overlapping those features:

```r
hse <- txdbmaker::makeTxDbFromGFF(ANNOT, format = "gtf")
exonsByGene <- exonsBy(hse, by = "gene")

bamLst <- BamFileList(fls, yieldSize = 100000)
se <- summarizeOverlaps(
    features = exonsByGene,
    reads = bamLst,
    mode = "Union",
    singleEnd = TRUE,
    ignore.strand = TRUE,
    fragments = FALSE
)
```

What each argument does:

- **GTF:** describes genomic features and their coordinates.
- **`exonsBy(..., by = "gene")`:** groups exon regions by gene ID.
- **`summarizeOverlaps()`:** produces a count matrix, with genes in rows and samples in columns.
- **`mode = "Union"`:** uses the union-style feature-overlap counting rule.
- **`singleEnd = TRUE`:** counts single-end alignments rather than paired-end fragments; this matches the setup used here.
- **`ignore.strand = TRUE`:** strand is not considered when assigning reads to genes. This is an analysis setting and should match the library preparation for a different dataset.
- **`yieldSize = 100000`:** processes the BAM data in chunks.

This workflow history does not mean that the matrix contains sex-adjusted values. `summarizeOverlaps()` produces raw counts; sex is a covariate in the subsequent statistical model, not a transformation of those counts. The current script 13 does not run this counting code. It reads the existing matrix at:

```text
$BASEDIR/results-proj1/DESeq2/Results/raw_counts.1mo_4mo.A.psiclass.csv
```

Although this source file is stored under the legacy `results-proj1/DESeq2/` path, its 14,214 gene IDs match the genome-wide PsiCLASS annotation and the prior adjusted result table. Its values are raw integer counts exported before DESeq2 model fitting; they were not adjusted for sex. The new unadjusted run exported an identical count matrix to `$BASEDIR/DESeq2_genomewide/unadjusted/Results/raw_counts.1mo_4mo.psiclass.csv`.

I had to call `makeTxDbFromGFF()` through **`txdbmaker`** in my installed Bioconductor version. The GTF file and the rest of the earlier counting approach stayed the same.

## 3. Metadata and the first design

The current scripts use `$BASEDIR/results-proj1/DESeq2/pheno_1mo_4mo.w-sex.csv`, which contains `sampleID`, `condition` and `sex`. Each `sampleID` is a full STAR BAM basename, such as `SRR3734796_Aligned.sortedByCoord.out.bam`. Script 13 extracts the SRA accession from the count-matrix and metadata names, matches on that accession, and reorders the metadata to the count columns. Script 13b matches the selected BAM basenames directly. Both preserve the original sample names for traceability.

In the initial course exercise, I also encountered a mismatch between the instructor script's `sampleInfo$sample` reference and a CSV header named `sampleID`. That was a course-file adaptation issue; it is not the metadata mapping used by the current genome-wide scripts.

Before modeling, the metadata rows must be matched to the count-matrix columns. Otherwise, the wrong condition could be assigned to a sample.

The current unadjusted analysis constructs the dataset from the existing raw-count matrix:

```r
dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = sampleInfo,
    design = ~ condition
)
dds$condition <- relevel(dds$condition, ref = "1mo")
dds <- DESeq2::DESeq(dds)
res <- DESeq2::results(dds, contrast = c("condition", "4mo", "1mo"))
```

- `design = ~ condition` models differences between age groups without other covariates.
- `ref = "1mo"` sets the reference group.
- `4mo` vs `1mo` means **positive** `log2FoldChange` corresponds to higher estimated expression in 4mo; a **negative** value corresponds to lower estimated expression in 4mo.

The count-based model is fitted by DESeq2, which estimates normalization factors, dispersions and differential expression using a negative-binomial framework.

## 4. Understanding the result columns

The current unadjusted result table is:

```text
$BASEDIR/DESeq2_genomewide/unadjusted/Results/results_1mo_4mo.psiclass.csv
```

| Field | Meaning |
| --- | --- |
| `baseMean` | Mean of normalized counts across samples |
| `log2FoldChange` | Estimated effect size for 4mo relative to 1mo |
| `lfcSE` | Standard error of the log2 fold-change estimate |
| `stat` | Statistical test statistic |
| `pvalue` | Unadjusted p-value |
| `padj` | P-value adjusted for multiple testing, using Benjamini–Hochberg by default |

Applying a `baseMean` threshold does **not** mean every sample must have that many reads. `baseMean` is the mean of normalized counts across samples.

A small p-value by itself is not enough to control multiple testing across thousands of genes. I use `padj` when interpreting statistically supported differential expression, while retaining the unadjusted p-value for the course's filtering exercise.

### Results from the current genome-wide unadjusted run

I reused the existing genome-wide raw-count matrix described above. I verified the result counts against the full result table and analysis summary, and confirmed that the exported raw counts are identical to the input matrix.

| Result | Genes |
| --- | ---: |
| Total genes in results table | 14,214 |
| `pvalue <= 0.05` | 1,313 |
| `padj <= 0.05` | 249 |
| Positive `log2FoldChange` among `padj <= 0.05` | 67 |
| Negative `log2FoldChange` among `padj <= 0.05` | 182 |
| Missing `pvalue` | 25 |
| Missing `padj` | 1,126 |

The script uses `results()` with its default settings. I do not assign an unverified cause to every missing p-value or adjusted p-value. In the exported table, 24 of the 25 rows with missing p-values have zero `baseMean`; one has a nonzero `baseMean`, and its cause is not established by the exported results. Missing adjusted p-values can arise under DESeq2's default independent filtering, but the exported table does not identify the reason for every missing value.

The script also exports normalized counts, a PCA plot and a volcano plot. The PCA uses a variance-stabilizing transformation (`vst`) to explore sample relationships. The volcano plot combines log2 fold change and statistical significance. These plots are exploratory and do not remove confounding by sex.

Script: `scripts/13_run_deseq2.R`  
Current output directory: `$BASEDIR/DESeq2_genomewide/unadjusted/`

## 5. Adding annotated gene names

The DESeq2 result tables use PsiCLASS gene IDs, which can be difficult to interpret biologically. The course provides `add_gn_name`, a Perl program that reads names from the annotation and adds them to a table; this step annotates results but does not rerun DESeq2.

The `gn_` prefix identifies a file with an added gene-name column. When I use such a table, I inspect its header and column order before filtering. In the annotated table I inspected, the order was:

```text
1: gene name       2: gene ID
3: baseMean        4: log2FoldChange
5: lfcSE           6: stat
7: pvalue          8: padj
```

When I read an annotated CSV with `check.names = FALSE`, the first two column headers can appear blank, so checking the field order matters. Some names may be labeled `novel_...`; adding a name column does not guarantee a familiar gene symbol for every entry.

The current script 13 outputs are not gene-name annotated. `scripts/15_add_gene_names.sh` reads and writes through `DESEQ2_RESULTS_DIR`, which is configured for the legacy Proj1 results directory. I do not treat that script as part of the current genome-wide workflow unless its paths are deliberately adapted.

Course references: `docs/com.add_gn_name` and `docs/add_gn_name`.

## 6. Testing different expression filters

I can narrow a result table using expression level, effect size and statistical significance, but I keep these criteria conceptually separate:

- `baseMean` is the mean of normalized counts across samples; it is not a minimum count required in every sample.
- An absolute `log2FoldChange` threshold selects changes in either direction. A larger absolute value represents a larger estimated effect.
- `pvalue` is unadjusted. `padj` accounts for multiple testing and is the significance measure I use for the current genome-wide result summary.

For the current analysis, I report the verified `pvalue <= 0.05` and `padj <= 0.05` totals above. I do not attach unverified gene counts to additional combinations of `baseMean` and effect-size cutoffs. More restrictive filters produce shorter lists, but that alone does not make a list more biologically valid.

After loading the current results table as `res`, I can define the two significance filters while excluding missing values explicitly:

```r
pvalue_selected <- !is.na(res$pvalue) & res$pvalue <= 0.05
padj_selected <- !is.na(res$padj) & res$padj <= 0.05
```

The course's `scripts/16_filter_deseq2.sh` demonstrates a grid of filters using comma-separated tables, but it reads and writes through `DESEQ2_RESULTS_DIR`, which is configured for the legacy Proj1 results directory. It is not a source of the current genome-wide results. I use `awk -F','` to specify CSV fields in that script, while remembering that `awk` is not a general CSV parser when quoted fields can contain commas; for complex CSVs, I use R or a dedicated CSV parser.

## 7. Biological sex as a covariate — follow-up analysis

At the end of the lecture, the instructor points out the unbalanced sex distribution:

```text
1mo: F, F, F
4mo: F, M, M
```

The original model was:

```r
design = ~ condition
```

The alternative model introduced in class is:

```r
design = ~ sex + condition
```

This model estimates a condition effect while including an **additive** sex term. It does not estimate a sex-by-condition interaction. Because there are no male samples at 1mo, this analysis cannot fully separate biological sex and age effects in a balanced way; results must be interpreted cautiously.

I prepared a separate script rather than modifying the original results:

```text
scripts/13b_run_deseq2_with_sex.R
```

The instructor's supplied alternative metadata file is `pheno_1mo_4mo.A.csv`, with `sampleID`, `condition` and `sex` (F/F/F/F/M/M for the six samples). My local copy uses a descriptive filename instead:

```text
$BASEDIR/results-proj1/DESeq2/pheno_1mo_4mo.w-sex.csv
```

For comparison with the course, the adjusted result files are intended to use the `.A.psiclass.csv` suffix, for example:

```text
results_1mo_4mo.A.psiclass.csv
```

An earlier sex-adjusted execution completed, and its verified `.A` result files remain
under the legacy `$BASEDIR/results-proj1/DESeq2/Results/` path. The corrected
`scripts/13b_run_deseq2_with_sex.R` is prepared to write to
`$BASEDIR/DESeq2_genomewide/sex_adjusted/`, but has not been rerun and has not produced
results there. Rerunning it is optional for improving reproducibility; it is not a
prerequisite for the existing exploratory comparison. Any new count matrix from that
rerun should be checked against the matrix used by the current unadjusted analysis.

## 8. Course commands compared with the local workflow

The course's `com.deseq` runs `R CMD BATCH deseq_1mo_4mo.R` after loading cluster modules. In my Ubuntu/WSL2 setup, the R packages are installed in the Conda environment rather than loaded through cluster `module load` commands.

The instructor's script exports raw counts, normalized counts, three DESeq2 result tables, PCA and a volcano plot. I retained these core outputs and their names. My local scripts add explicit paths and metadata/sample-order checks; they are adaptations of the course exercise, not identical copies of the original script. The current script 13 reads the pre-existing raw-count matrix instead of recounting BAMs.

## 9. Problems encountered and fixes

### Running R from the wrong directory

I first ran:

```bash
Rscript ./scripts/13_run_deseq2.R
```

from `Module_I/`. The script could read the GTF but then failed with:

```text
cannot open file 'pheno_1mo_4mo.csv': No such file or directory
```

Why: `read.csv("pheno_1mo_4mo.csv")` looked in the **current working directory**, not in the script's folder. An earlier version also wrote a stray `raw_counts.1mo_4mo.psiclass.csv` into `Module_I/` before failing.

Fix: define full input/output paths (`PHENO`, `OUTDIR`, `RESULTSDIR`) so execution does not depend on running the script from the results directory. This was an issue in an earlier script version; the current script validates its inputs before creating output directories and writes to the genome-wide unadjusted destination.

### `txdbmaker` versus the course code

The original course script loads `DESeq2`, `GenomicFeatures`, `Rsamtools`, `GenomicAlignments` and `EnhancedVolcano` and calls `makeTxDbFromGFF()` directly. It also uses `singleEnd = FALSE, fragments = TRUE` in `summarizeOverlaps()`. My local version uses `singleEnd = TRUE, fragments = FALSE` because I am treating the six runs as single-end data. This changes how reads are counted, so it is more than a path change. The course calls `makeTxDbFromGFF()` as part of the GenomicFeatures workflow. With the Bioconductor version in my environment, I used:

```r
txdbmaker::makeTxDbFromGFF(ANNOT, format = "gtf")
```

This explicitly calls the package that provides the function in my installed version.

### VS Code R extension and `jsonlite`

The VS Code R extension initially reported that `jsonlite` was missing. After installing it in the Conda environment and checking that the editor used the environment's R executable, the packages loaded successfully.

```bash
which R
# $TRANSCRIPTOMICS_ENV/bin/R
```

I also accidentally entered `library(jsonlite)` and `library(DESeq2)` directly at the Bash prompt. Those are **R expressions**, so I needed to start `R` first, or run them through `Rscript -e '...'`.

The message that some R packages were built under R 4.5.3 while I was using R 4.5.2 was a **warning**, not a failure in the successful package-loading test.

## 10. Repository and result locations

The checked-in repository keeps the scripts and notes separate from the machine-specific analysis data:

```text
repository/
├── configs/
│   └── paths.example.sh
├── course_materials/Module_I/
│   └── notes/
└── scripts/
    ├── 11_sort_bam_by_read_name.sh
    ├── 12_check_sorted_bam.sh
    ├── 13_run_deseq2.R
    ├── 13b_run_deseq2_with_sex.R
    ├── 14_summarize_deseq2.sh
    ├── 15_add_gene_names.sh
    └── 16_filter_deseq2.sh

$BASEDIR/
├── STAR/BAM/
├── PsiCLASS/WithGeneNames/
├── results-proj1/
│   ├── Sorted/  # current regional Proj1 BAMs
│   └── DESeq2/
│       └── Results/  # legacy inputs and previously verified sex-adjusted results
└── DESeq2_genomewide/
    ├── unadjusted/  # completed script 13 outputs
    │   └── Results/
    └── sex_adjusted/  # planned destination; not created, corrected script not rerun
```

Script 13 completed successfully, and its outputs are available under
`$BASEDIR/DESeq2_genomewide/unadjusted/`. The corrected script 13b is prepared to write
to `$BASEDIR/DESeq2_genomewide/sex_adjusted/`, but has not been rerun; this target
therefore has no results from the corrected script. The previously verified
sex-adjusted results remain under `$BASEDIR/results-proj1/DESeq2/Results/`. The
raw-count matrix and metadata used by script 13 remain read-only inputs at their
existing legacy paths.

Large BAMs, count tables, and figures are stored under `$BASEDIR` rather than committed. The local `configs/paths.sh` is machine-specific and is not tracked; create it from `configs/paths.example.sh` and set `BASEDIR` for the machine used.

## Next steps

- If improved reproducibility is needed, rerun the corrected `13b_run_deseq2_with_sex.R` and verify that its regenerated raw counts match the matrix used here.
- Keep any comparison with the previous sex-adjusted output exploratory and compare gene IDs, not only the total number of significant genes.
- Add gene names to the current outputs only after confirming the annotation and gene-ID mapping for this 14,214-gene universe.
- Interpret any biological findings in light of the sex imbalance and small sample size.

See [Note 4 — Sex-adjusted DESeq2 and Method Comparison](04_sex_adjusted_deseq2_and_method_comparison.md) for the run prerequisites and comparison checks.

*Course exercise: Johns Hopkins University / Coursera. Commands, adjustments and results above are from my local run.*
