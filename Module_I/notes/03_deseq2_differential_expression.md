# Module I — DESeq2 differential expression analysis

**Course:** Bioinformatics Methods for Transcriptomics — Johns Hopkins University (Coursera)  
**Dataset:** GSE83931 — mouse hippocampus RNA-seq  
**Comparison:** 4 months (4mo) vs 1 month (1mo)  
**Working environment:** Ubuntu / WSL2, Conda (`transcriptomics`), R 4.5.2

I continued with the same six RNA-seq samples used for STAR, PsiCLASS and Cuffdiff2. For this part of the course, I followed the DESeq2 exercise and changed the paths and a few settings to run it locally.

The six samples are:

| Condition | SRA runs | Biological sex |
| --- | --- | --- |
| 1mo | SRR3734796, SRR3734797, SRR3734798 | Female, female, female |
| 4mo | SRR3734820, SRR3734824, SRR3734825 | Female, male, male |

**Dataset limitation:** the 1mo samples are all female, whereas the 4mo group contains one female and two males. My first DESeq2 model does not account for sex. Adding sex to the model later does not fix the lack of male 1mo samples.

## 1. Preparing BAM files for gene counting

STAR produced BAM files sorted by genomic coordinate, for example:

```text
results/STAR/BAM/SRR3734796_Aligned.sortedByCoord.out.bam
```

For this exercise, I also created **queryname-sorted** BAM files using SAMtools:

```bash
samtools sort -n -@ 2 -O BAM \
    -o results/Sorted/SRR3734796.so.bam \
    results/STAR/BAM/SRR3734796_Aligned.sortedByCoord.out.bam
```

- `sort` rearranges alignment **records**; it does not realign reads or change their genomic coordinates.
- Without `-n`, `samtools sort` sorts by genomic coordinates (reference and position).
- `-n` sorts by read name (QNAME). Alignments belonging to the same read are grouped together.
- `-@ 2` requests two additional processing threads.
- `-O BAM` specifies BAM output; `-o` gives the destination file.

A useful example from the actual data was **`SRR3734796.sra.8`**, which had two reported alignments:

| Read | Genomic position | Note |
| --- | --- | --- |
| `SRR3734796.sra.8` | chr1:95426560 | Primary alignment, FLAG 0 |
| `SRR3734796.sra.8` | chr8:71348574 | Secondary reverse-strand alignment, FLAG 272 |

Sorting by QNAME keeps these records together. Coordinate sorting places them according to their genomic locations. The read was reported as mapping to two locations (`NH:i:2`); this illustrates **multimapping**.

Scripts:

- `scripts/11_sort_bam_by_read_name.sh`
- `scripts/12_check_sorted_bam.sh`

I checked all six output BAMs with `samtools quickcheck` and inspected their headers for `SO:queryname`. These checks showed that the files were readable and marked as queryname-sorted. They do not check every record individually.

## 2. From alignments to a gene-count matrix

DESeq2 uses a matrix of **raw gene counts** rather than the FPKM-style expression estimates produced by Cuffdiff2.

The gene annotation came from the PsiCLASS step:

```text
results/PsiCLASS/WithGeneNames/psiclass_vote.withStrand.gtf
```

The R script constructs a transcript database, groups exons by gene, and counts alignments overlapping those gene features:

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

The raw counts were exported as:

```text
results/DESeq2/Results/raw_counts.1mo_4mo.psiclass.csv
```

I had to call `makeTxDbFromGFF()` through **`txdbmaker`** in my installed Bioconductor version. The GTF file and the rest of the counting approach stayed the same.

## 3. Metadata and the first design

The initial phenotype file, `pheno_1mo_4mo.csv`, contains `sampleID` and `condition`. In the instructor's R script the matching expression refers to `sampleInfo$sample`, whereas the supplied CSV header is `sampleID`. My local adaptation matches on `sampleInfo$sampleID`, which is consistent with the metadata file. Its sample IDs match the name-sorted BAM filenames, such as `SRR3734796.so.bam`.

Before modeling, the metadata rows must be matched to the count-matrix columns. Otherwise, the wrong condition could be assigned to a sample.

The original analysis uses:

```r
dds <- DESeqDataSet(se, design = ~ condition)
dds$condition <- relevel(dds$condition, ref = "1mo")
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "4mo", "1mo"))
```

- `design = ~ condition` models differences between age groups without other covariates.
- `ref = "1mo"` sets the reference group.
- `4mo` vs `1mo` means **positive** `log2FoldChange` corresponds to higher estimated expression in 4mo; a **negative** value corresponds to lower estimated expression in 4mo.

The count-based model is fitted by DESeq2, which estimates normalization factors, dispersions and differential expression using a negative-binomial framework.

## 4. Understanding the result columns

The original result table is:

```text
results/DESeq2/Results/results_1mo_4mo.psiclass.csv
```

| Field | Meaning |
| --- | --- |
| `baseMean` | Mean of normalized counts across samples |
| `log2FoldChange` | Estimated effect size for 4mo relative to 1mo |
| `lfcSE` | Standard error of the log2 fold-change estimate |
| `stat` | Statistical test statistic |
| `pvalue` | Unadjusted p-value |
| `padj` | P-value adjusted for multiple testing, using Benjamini–Hochberg by default |

`baseMean >= 10` does **not** mean every sample has at least 10 reads. It means the mean normalized count across samples meets that threshold.

A small p-value by itself is not enough to control multiple testing across thousands of genes. I use `padj` when interpreting statistically supported differential expression, while retaining the unadjusted p-value for the course's filtering exercise.

### Results from my original DESeq2 run

| Result | Genes |
| --- | ---: |
| Total genes in results table | 14,293 |
| `pvalue <= 0.05` | 1,328 |
| `padj <= 0.05` | 271 |
| Upregulated among `padj <= 0.05` | 77 |
| Downregulated among `padj <= 0.05` | 194 |

My totals are different from the instructor's. The PsiCLASS annotation and counting settings may help explain the difference, but I have not tested that yet. Cuffdiff2 and DESeq2 also gave similar numbers of significant genes; I still need to check whether they identified the same genes.

The script also exports normalized counts, a PCA plot and a volcano plot. The PCA uses a variance-stabilizing transformation (`vst`) to explore sample relationships. The volcano plot combines log2 fold change and statistical significance. These plots are exploratory and do not remove confounding by sex.

Script: `scripts/13_run_deseq2.R`  
Summary script: `scripts/14_summarize_deseq2.sh`

## 5. Adding annotated gene names

The first DESeq2 tables use PsiCLASS gene identifiers, which are not always easy to interpret biologically. The course provides `add_gn_name`, a Perl program that reads gene names from the annotation and adds them to the output tables.

I kept the five original tables in:

```text
results/DESeq2/Results/
```

And generated annotated copies in:

```text
results/DESeq2/Results/WithGeneNames/
```

Example:

```text
results_1mo_4mo.psiclass.csv
        ↓
gn_results_1mo_4mo.psiclass.csv
```

The prefix **`gn_`** identifies files with an added gene-name column. This only adds gene names to the tables. It does not rerun DESeq2.

Script: `scripts/15_add_gene_names.sh`  
Course reference: `docs/com.add_gn_name` and `docs/add_gn_name`

For the `gn_results_1mo_4mo.psiclass.csv` file, the inspected column order was:

```text
1: gene name       2: gene ID
3: baseMean        4: log2FoldChange
5: lfcSE           6: stat
7: pvalue          8: padj
```

When I read the CSV with `check.names = FALSE`, the first two column headers appeared blank. I therefore checked the column order before filtering with `awk`. I checked the uploaded annotated tables: they contain 14,293 data rows, with a gene name and PsiCLASS gene ID in the first two columns. Some names are still labeled `novel_...`, which is expected for transcripts without a familiar gene symbol; adding a column does not guarantee a recognized gene name for every entry.

## 6. Testing different expression filters

The instructor showed how the number of selected genes changes when applying thresholds to expression level, fold change and significance. I reproduced this exercise with `scripts/16_filter_deseq2.sh`.

Values tested:

- Minimum `baseMean`: 0, 10, 15
- Minimum **absolute** log2 fold change: 0, 0.5, 1
- Significance: `pvalue <= 0.05` or `padj <= 0.05`

| Minimum baseMean | Minimum \|log2FC\| | p-value <= 0.05 | padj <= 0.05 |
| ---: | ---: | ---: | ---: |
| 0 | 0 | 1,328 | 271 |
| 0 | 0.5 | 362 | 156 |
| 0 | 1 | 108 | 50 |
| 10 | 0 | 1,307 | 271 |
| 10 | 0.5 | **341** | **156** |
| 10 | 1 | 87 | 50 |
| 15 | 0 | 1,294 | 271 |
| 15 | 0.5 | 328 | 156 |
| 15 | 1 | 74 | 50 |

In the video, the instructor obtained about 380 genes with `baseMean >= 10`, `|log2FC| >= 0.5` and unadjusted `pvalue <= 0.05`; my dataset gave **341** with those same cutoffs. I checked all nine rows in my summary against the uploaded `gn_results_1mo_4mo.psiclass.csv`, and the counts match.

The reference files provided with the course show additional examples rather than this whole grid. In `com.deseq-filter`, the instructor applies `baseMean >= 15` and `|log2FC| >= 1` to the complete, p-value-filtered and q-value-filtered files; `deseq.filter` shows a separate example using `baseMean >= 10`, `|log2FC| >= 1` and `pvalue <= 0.05`. My script extends these examples to compare multiple thresholds in one table.

The instructor uses `tr ',' ' '` before `awk` because the exported tables are comma-separated. My script uses `awk -F','` for the same numeric fields instead.

Observations from the table:

1. The `baseMean` cutoff changed the unadjusted-p-value count only slightly and did not change the 271 genes selected by `padj <= 0.05` when applied alone.
2. The fold-change cutoff made a larger difference: among genes with `padj <= 0.05`, 156 remained at `|log2FC| >= 0.5`, and 50 at `|log2FC| >= 1`.
3. More restrictive filters produce smaller gene lists, but this alone does not make a list more biologically valid.

For interpretation, `|log2FC| >= 1` represents a change of at least twofold in either direction; `|log2FC| >= 0.5` is about 1.41-fold.

The summary is stored locally at:

```text
results/DESeq2/Results/deseq_filter_summary.csv
```

I used `awk -F','` in the filtering script to specify the CSV separator. This works with the numeric columns in these particular files, but `awk` is not a general CSV parser when quoted fields can contain commas. For more complex CSVs, I would use R or a dedicated CSV parser.

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
results/DESeq2/pheno_1mo_4mo.w-sex.csv
```

For comparison with the course, the adjusted result files are intended to use the `.A.psiclass.csv` suffix, for example:

```text
results_1mo_4mo.A.psiclass.csv
```

**Status when writing this note:** I prepared the script and metadata for the sex-adjusted analysis, but I have not checked the `.Rout` or its result tables yet. I will add the gene counts after verifying the run.

## 8. Course commands compared with the local workflow

The course's `com.deseq` runs `R CMD BATCH deseq_1mo_4mo.R` after loading cluster modules. In my Ubuntu/WSL2 setup, the R packages are installed in the Conda environment rather than loaded through cluster `module load` commands.

The instructor's script exports the raw counts, normalized counts, three DESeq2 result tables, PCA and volcano plot. I retained these core outputs and their names. My local scripts add explicit paths and metadata/sample-order checks; they are adaptations of the course exercise, not identical copies of the original script.

Instead of running each filter separately, I used `16_filter_deseq2.sh` to save the nine combinations in `deseq_filter_summary.csv`. This made it easier to compare the cutoffs.

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

Fix: define full input/output paths (`PHENO`, `OUTDIR`, `RESULTSDIR`) so execution does not depend on running the script from the results directory. The stray count file should be compared against the expected output before removal.

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
# /home/ianca/miniconda3/envs/transcriptomics/bin/R
```

I also accidentally entered `library(jsonlite)` and `library(DESeq2)` directly at the Bash prompt. Those are **R expressions**, so I needed to start `R` first, or run them through `Rscript -e '...'`.

The message that some R packages were built under R 4.5.3 while I was using R 4.5.2 was a **warning**, not a failure in the successful package-loading test.

## 10. File organization and version control

The relevant part of the local project is:

```text
Module_I/
├── scripts/
│   ├── 11_sort_bam_by_read_name.sh
│   ├── 12_check_sorted_bam.sh
│   ├── 13_run_deseq2.R
│   ├── 13b_run_deseq2_with_sex.R
│   ├── 14_summarize_deseq2.sh
│   ├── 15_add_gene_names.sh
│   └── 16_filter_deseq2.sh
├── notes/
│   └── 03_deseq2_differential_expression.md
├── docs/                     # Course reference materials; not versioned
└── results/                  # Symlink to local analysis storage
    └── DESeq2/
        ├── Results/
        │   ├── *.psiclass.csv
        │   ├── deseq_filter_summary.csv
        │   └── WithGeneNames/
        │       └── gn_*.psiclass.csv
        ├── *.png
        └── *.Rout
```

I keep large BAMs, count tables and figures on the D: drive and commit the scripts and notes to GitHub. Some scripts still use absolute WSL paths, so they would need path changes on another machine.

## Next steps

- Verify the run of `13b_run_deseq2_with_sex.R` and compare its results with the original model.
- Add gene names to the `.A.psiclass.csv` result files once the sex-adjusted run is confirmed.
- Examine which gene IDs overlap between the original and adjusted results, not only the total number of significant genes.
- Interpret any biological findings in light of the sex imbalance and small sample size.

*Course exercise: Johns Hopkins University / Coursera. Commands, adjustments and results above are from my local run.*
