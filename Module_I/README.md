# Module I — RNA-seq Analysis

**RNA-seq Data Analysis | STAR | PsiCLASS | Cuffdiff2 | DESeq2 | R | Bash**

Scripts and study notes developed while following the **Bioinformatics Methods for Transcriptomics** course offered by Johns Hopkins University through Coursera.

**Course:** [Bioinformatics Methods for Transcriptomics](https://www.coursera.org/learn/bioinformatics-methods-for-transcriptomics)  
**Instructor:** Liliana Florea, PhD  
**Status:** In progress — Module I

## About this module

I worked through the first course module using six mouse hippocampus RNA-seq runs. The analysis starts with SRA downloads and FASTQ conversion, then moves through STAR alignment, PsiCLASS transcript assembly, Cuffdiff2, and DESeq2.

The course commands were designed for a computing cluster, so I had to change paths and some settings to run them on Ubuntu/WSL2. I kept the main steps close to the tutorial and saved my changes in numbered scripts. The notes document what I ran, what went wrong, and how I checked the outputs.

## Dataset and Experimental Design

Module I uses publicly available mouse hippocampus RNA-seq data from [GEO accession GSE83931](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE83931).

Six sequencing runs were analyzed:

| Condition | SRA runs | Biological sex |
|---|---|---|
| 1 month | SRR3734796, SRR3734797, SRR3734798 | Female, female, female |
| 4 months | SRR3734820, SRR3734824, SRR3734825 | Female, male, male |

**Limitation:** all three 1-month samples are female, while two of the 4-month samples are male. Differences attributed to age may also be influenced by sex. I wrote a separate DESeq2 script using `~ sex + condition`, but its results have not yet been checked.

## Workflow

The following steps have been performed locally:

1. **Data retrieval:** selected SRA accessions, downloaded sequencing data with `prefetch`, and converted the files to FASTQ with `fasterq-dump`.
2. **Reference preparation:** used the GRCm39/mm39 mouse genome (GENCODE mouse release M30) and built a memory-reduced STAR index.
3. **Read alignment:** aligned all six samples with STAR, producing coordinate-sorted BAM files; inspected alignment statistics and splice junctions.
4. **Transcript assembly:** assembled transcripts with PsiCLASS and produced an annotated GTF for downstream analyses.
5. **Cuffdiff2:** performed differential expression analysis using the assembled transcript annotation.
6. **BAM sorting and gene counts:** prepared queryname-sorted BAM files with SAMtools and counted reads overlapping annotated genes using `summarizeOverlaps()`.
7. **DESeq2:** compared 4-month and 1-month samples with the model `~ condition`, using 1 month as the reference.
8. **Gene annotation and filtering:** added gene names to DESeq2 output tables and explored thresholds for `baseMean`, absolute log2 fold change, unadjusted p-values, and adjusted p-values.

An additional script, `13b_run_deseq2_with_sex.R`, is prepared for the alternative model `~ sex + condition`. Its successful execution and results still need to be checked.

### STAR Alignment: First Sample

For **SRR3734796**, the STAR alignment summary reported:

| Alignment metric | Result |
|---|---:|
| Input reads | 11,851,553 |
| Average read length | 100 bp |
| Uniquely mapped reads | 89.79% |
| Reads mapped to multiple loci | 6.96% |
| Reads unmapped: too short | 2.15% |
| Mismatch rate per base | 0.27% |
| Splice-junction occurrences | 1,641,050 |

These statistics describe **one sample**, not the aggregate across all six runs.

The STAR index was built from the reference FASTA without incorporating a GTF file, using reduced-memory parameters including `--genomeSAsparseD 16` and `--genomeSAindexNbases 12`.

### Differential Expression Results

The following values come from the completed analyses of the six selected samples, using the PsiCLASS-based annotation.

| Metric | Cuffdiff2 | DESeq2 (`~ condition`) |
|---|---:|---:|
| Genes significant after multiple-testing adjustment | 265 | 271 |
| Upregulated genes (4mo vs 1mo) | 63 | 77 |
| Downregulated genes (4mo vs 1mo) | 202 | 194 |

In the DESeq2 results table, **14,293 genes** were reported, with **1,328 genes at `pvalue ≤ 0.05`** and **271 at `padj ≤ 0.05`**. These results use Benjamini–Hochberg-adjusted p-values for the significance comparison.

The Cuffdiff2 and DESeq2 counts should not be treated as directly interchangeable: the methods use different expression estimation and statistical procedures, and their significant gene sets have not yet been compared gene by gene.

#### Exploring DESeq2 Filters

After adding gene names to the result files, I tested combinations of minimum mean normalized count (`baseMean`), minimum absolute log2 fold change, and statistical significance:

| Minimum `baseMean` | Minimum `|log2FC|` | `pvalue ≤ 0.05` | `padj ≤ 0.05` |
|---:|---:|---:|---:|
| 0 | 0 | 1,328 | 271 |
| 0 | 0.5 | 362 | 156 |
| 0 | 1 | 108 | 50 |
| 10 | 0 | 1,307 | 271 |
| 10 | 0.5 | 341 | 156 |
| 10 | 1 | 87 | 50 |
| 15 | 0 | 1,294 | 271 |
| 15 | 0.5 | 328 | 156 |
| 15 | 1 | 74 | 50 |

The `baseMean` threshold had relatively little effect on the adjusted-significance counts in this dataset. Increasing the absolute log2 fold change threshold reduced the selected gene set more substantially. These are exploratory **post-analysis filters**, not new DESeq2 model fits.

## Module Structure

```text
Module_I/
├── README.md               # This module overview
├── data/                   # Local links to SRA and FASTQ data
├── docs/                   # Course reference files (gitignored)
├── metadata/               # SRA and sample metadata
├── notes/                  # Study notes
├── scripts/                # Numbered analysis scripts
└── results/                # Local link to outputs on another drive
```

Large datasets, reference files, indexes, BAMs, and generated analysis results are stored locally, rather than committed to GitHub. Paths in some scripts are specific to the local WSL2 setup and must be adapted on other machines.

### Main Scripts

| Script | Purpose |
|---|---|
| `01_extract_sra_accessions.sh` | Extract run identifiers from metadata |
| `02_download_sra.sh` | Download SRA runs |
| `03_convert_sra_to_fastq.sh` | Convert SRA to FASTQ |
| `04_generate_star_index.sh` | Build the STAR genome index |
| `05a_test_star_alignment.sh` | Test alignment on one sample |
| `05b_align_all_samples_star.sh` | Align all six samples |
| `06a_run_psiclass.sh` | Assemble transcripts with PsiCLASS |
| `06b_download_reference_annotation.sh` | Prepare reference annotation |
| `07_annotate_transcripts.sh` | Annotate assembled transcripts |
| `08_summarize_psiclass.sh` / `.py` | Summarize assembly outputs |
| `09_run_cuffdiff.sh` | Run Cuffdiff2 |
| `10_summarize_cuffdiff.sh` | Summarize Cuffdiff2 results |
| `11_sort_bam_by_read_name.sh` | Sort BAM records by read name |
| `12_check_sorted_bam.sh` | Validate name-sorted BAM files |
| `13_run_deseq2.R` | Run the original DESeq2 comparison |
| `13b_run_deseq2_with_sex.R` | Prepare a sex-adjusted DESeq2 comparison (results pending verification) |
| `14_summarize_deseq2.sh` | Summarize DESeq2 results |
| `15_add_gene_names.sh` | Add annotation-derived gene names to result tables |
| `16_filter_deseq2.sh` | Compare differential-expression filter thresholds |

### Local DESeq2 Output Layout

```text
results/DESeq2/
├── pheno_1mo_4mo.csv
├── pheno_1mo_4mo.w-sex.csv
├── deseq_1mo_4mo.Rout
├── PCA_1mo_4mo_byCondition.png
├── volcano_1mo_4mo.psiclass.png
└── Results/
    ├── raw_counts.1mo_4mo.psiclass.csv
    ├── normalized_counts.1mo_4mo.psiclass.csv
    ├── results_1mo_4mo.psiclass.csv
    ├── results_1mo_4mo.pval05.psiclass.csv
    ├── results_1mo_4mo.qval05.psiclass.csv
    ├── deseq_filter_summary.csv
    └── WithGeneNames/
        └── gn_*.psiclass.csv
```

The `gn_` tables are created separately using the course-provided `add_gn_name` utility and the PsiCLASS annotation. The `.A.psiclass.csv` files are reserved for the alternative analysis including sex as a covariate.

## Study Notes

- [RNA-seq Data Preparation and STAR Alignment](notes/01_rnaseq_data_preparation_and_star_alignment.md)
- [Cuffdiff2 Differential Expression](notes/02_cuffdiff2_differential_expression.md)

The notebooks record commands, decisions, errors, and explanations from the exercises. Additional notes on DESeq2 and gene annotation will be added as the module progresses.

## Tools and Technologies

**Programming and environment:** Bash, R, Python, Ubuntu/WSL2, Conda, Git/GitHub.  
**Sequencing and alignment:** NCBI SRA Toolkit, STAR, SAMtools.  
**Transcriptomics:** PsiCLASS, Cuffdiff2, DESeq2, Bioconductor (`GenomicFeatures`, `GenomicAlignments`, `Rsamtools`), `txdbmaker`.  
**Analysis and visualization:** read-count matrices, differential expression statistics, BH/FDR correction, PCA, volcano plots, R/ggplot2.

## What remains to check

This is a course exercise, not a validated biological study. The main Module I analyses have run, but I have not yet checked the results from the sex-adjusted DESeq2 script or compared the lists of significant genes from Cuffdiff2 and DESeq2. The imbalance in sex between age groups remains an important limitation.

## Acknowledgements

This repository is based on [Bioinformatics Methods for Transcriptomics](https://www.coursera.org/learn/bioinformatics-methods-for-transcriptomics), taught by Dr. Liliana Florea at Johns Hopkins University via Coursera. Course materials are used as learning references, with local adaptations to scripts and data organization. This is an independent study repository and is not an official university or Coursera resource.

## Back to the Repository

[Main repository README](../README.md)
