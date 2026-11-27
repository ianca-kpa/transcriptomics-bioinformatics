# Module I — RNA-seq Analysis

**RNA-seq Data Analysis | STAR | PsiCLASS | Cuffdiff2 | DESeq2 | R | Bash**

Scripts and study notes developed while following the **Bioinformatics Methods for Transcriptomics** course offered by Johns Hopkins University through Coursera.

**Course:** [Bioinformatics Methods for Transcriptomics](https://www.coursera.org/learn/bioinformatics-methods-for-transcriptomics)  
**Instructor:** Liliana Florea, PhD  
**Status:** Main genome-wide analyses completed; corrected sex-adjusted rerun optional

## About this module

I worked through the first course module using six mouse hippocampus RNA-seq runs. The analysis starts with SRA downloads and FASTQ conversion, then moves through STAR alignment, PsiCLASS transcript assembly, Cuffdiff2, and DESeq2.

The course commands were designed for a computing cluster, so I adapted paths and some settings for Ubuntu/WSL2. I kept the main steps close to the tutorial and saved my changes in numbered scripts. The notes document what I ran, what went wrong, and how I checked the outputs.

## Dataset and Experimental Design

Module I uses publicly available mouse hippocampus RNA-seq data from [GEO accession GSE83931](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE83931).

Six sequencing runs were analyzed:

| Condition | SRA runs | Biological sex |
|---|---|---|
| 1 month | SRR3734796, SRR3734797, SRR3734798 | Female, female, female |
| 4 months | SRR3734820, SRR3734824, SRR3734825 | Female, male, male |

**Design limitation:** all three 1-month samples are female; the 4-month group has one female and two males. The additive DESeq2 model `~ sex + condition` does not resolve this imbalance or estimate a sex-by-condition interaction. I interpret both age comparisons cautiously and do not treat them as causal evidence.

## Workflow

I completed the main genome-wide workflow locally:

1. **Data retrieval:** selected SRA accessions, downloaded sequencing data with `prefetch`, and converted the files to FASTQ with `fasterq-dump`.
2. **Reference preparation:** used the GRCm39/mm39 mouse genome (GENCODE mouse release M30) and built a memory-reduced STAR index.
3. **Read alignment:** aligned all six samples with STAR, producing coordinate-sorted BAM files; inspected alignment statistics and splice junctions.
4. **Transcript assembly:** assembled transcripts with PsiCLASS and produced an annotated GTF for downstream analyses.
5. **Cuffdiff2:** analyzed the six coordinate-sorted genome-wide BAMs against the PsiCLASS annotation; current results are in `$BASEDIR/Cuffdiff/`.
6. **Unadjusted DESeq2:** script 13 imported the existing raw-count matrix and fit `~ condition`, using 1 month as the reference. It did not recount BAMs.
7. **Sex-adjusted DESeq2:** I verified results from a previous `~ sex + condition` execution. The corrected script 13b has not been rerun; a future run is optional for improving reproducibility.

The raw-count matrix contains counts, not values adjusted for sex. The previous counting workflow used genome-wide BAMs and `summarizeOverlaps()` before fitting a DESeq2 model. Scripts 11 and 12 describe a separate BAM-sorting exercise; their configured destination currently points to the regional Proj1 directory and is not an input to script 13.

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

### Genome-wide Differential Expression Results

I report the current verified Cuffdiff2 and unadjusted DESeq2 results, alongside the previously executed and verified sex-adjusted DESeq2 results. Cuffdiff2 significance requires `status == "OK"` and `significant == "yes"`; DESeq2 significance uses `padj <= 0.05`.

| Analysis | Gene records | Significant genes |
|---|---:|---:|
| Cuffdiff2 | 14,214 | 261 |
| DESeq2 unadjusted (`~ condition`) | 14,214 | 249 |
| Previously executed DESeq2 sex-adjusted (`~ sex + condition`) | 14,214 | 124 |

The current Cuffdiff2 results also contain **207 significant isoforms**. The unadjusted DESeq2 result table has **1,313 genes with `pvalue <= 0.05`** and 249 with `padj <= 0.05`.

By gene ID, the two DESeq2 models share **85 significant genes**. Cuffdiff2 and the previously executed sex-adjusted DESeq2 share **76 significant genes**. These are exploratory method comparisons, not evidence of biological validity; I do not infer a sex effect from changes in statistical significance. The full comparisons and their limitations are documented in [Note 04](notes/04_sex_adjusted_deseq2_and_method_comparison.md).

The previously executed adjusted results remain under `$BASEDIR/results-proj1/DESeq2/Results/`; their environment provenance is incomplete. The corrected script 13b has not been rerun and has not produced results in its planned destination.

### Inputs and Output Locations

| Workflow | Inputs | Results or status |
|---|---|---|
| Cuffdiff2 | `$BASEDIR/STAR/BAM/` and `$BASEDIR/PsiCLASS/WithGeneNames/psiclass_vote.withStrand.gtf` | Current results: `$BASEDIR/Cuffdiff/` |
| DESeq2 unadjusted | Existing raw-count matrix: `$BASEDIR/results-proj1/DESeq2/Results/raw_counts.1mo_4mo.A.psiclass.csv` | Completed outputs: `$BASEDIR/DESeq2_genomewide/unadjusted/` |
| Previously executed DESeq2 sex-adjusted | Previously verified results and execution log | `$BASEDIR/results-proj1/DESeq2/Results/` |
| Corrected script 13b | Genome-wide BAMs, PsiCLASS GTF, and validated sex metadata | Planned destination: `$BASEDIR/DESeq2_genomewide/sex_adjusted/`; not rerun or created |

The corrected script 13b will count reads with `summarizeOverlaps()` before fitting the sex-adjusted model; the count matrix itself is not statistically adjusted for sex. The script creates its output directory and refuses to use an existing destination, so I will not create that directory in advance. For explanations of `baseMean`, log2 fold change, p-values, and FDR filtering, see [Note 02](notes/02_cuffdiff2_differential_expression.md) and [Note 03](notes/03_deseq2_differential_expression.md).

## Module Structure

```text
course_materials/Module_I/
├── README.md               # This module overview
├── ENVIRONMENT.md          # Local environment and dependency notes
├── docs/                   # Course reference files and utilities
├── metadata/               # SRA and sample metadata
├── notes/                  # Study notes
├── Final_Project/          # Regional chromosome 7 notebook project
└── raw_counts.*.csv        # Selected count exports

scripts/                    # Numbered repository-level analysis scripts
configs/                    # Local path configuration and example
```

Large datasets, indexes, BAMs, and generated analysis results are stored locally rather than committed to GitHub. The numbered workflows in `scripts/` load machine-specific paths from `configs/paths.sh` and must be adapted on another machine. Script 01 extracts accessions from the Module I SRA table; script 02 defaults to its six-sample list; script 03 defaults to the Module II list and must be given the Module I list for conversion. A Module I six-sample batch-alignment script is still needed; see [the data preparation and STAR note](notes/01_rnaseq_data_preparation_and_star_alignment.md).

### Main Scripts

| Script | Purpose |
|---|---|
| `01_extract_sra_accessions.sh` | Extract run identifiers from Module I metadata |
| `02_download_sra.sh` | Download SRA runs (defaults to Module I; `SRA_RUNS_FILE` overrides) |
| `03_convert_sra_to_fastq.sh` | Convert SRA to FASTQ (defaults to Module II; `SRA_RUNS_FILE` overrides) |
| `04_generate_star_index.sh` | Build the STAR genome index |
| `05a_test_star_alignment.sh` | Test alignment on one Module I sample |
| `05b_align_all_samples_star.sh` | Module II batch alignment; not the Module I six-run workflow |
| `06a_run_psiclass.sh` | Assemble transcripts with PsiCLASS |
| `06b_download_reference_annotation.sh` | Prepare reference annotation |
| `07_annotate_transcripts.sh` | Annotate assembled transcripts |
| `08_summarize_psiclass.sh` / `.py` | Summarize assembly outputs |
| `09_run_cuffdiff.sh` | Run Cuffdiff2 to `$CUFFDIFF_DIR`; check existing results before rerunning because outputs may be replaced |
| `10_summarize_cuffdiff.sh` | Summarize Cuffdiff2 results |
| `11_sort_bam_by_read_name.sh` | Earlier queryname-sorting exercise; current `SORTED_BAM_DIR` points to regional Proj1 BAMs, not inputs to script 13 |
| `12_check_sorted_bam.sh` | Check BAM integrity and declared sort order in the configured directory; this does not establish provenance |
| `13_run_deseq2.R` | Completed unadjusted analysis; imports the existing raw-count matrix and writes to `DESeq2_genomewide/unadjusted/` |
| `13b_run_deseq2_with_sex.R` | Corrected sex-adjusted workflow prepared for a future run; not rerun |
| `14_summarize_deseq2.sh` | Summarizes files under `DESEQ2_RESULTS_DIR`, currently the legacy Proj1 results directory |
| `15_add_gene_names.sh` | Adds gene names under the configured legacy DESeq2 results directory; not wired to current genome-wide outputs |
| `16_filter_deseq2.sh` | Demonstrates filter thresholds under the configured legacy results directory; no current genome-wide filter counts are reported here |

### Regional Final Project

The Final Project notebook analyzes approximately 3.4 Mb of chromosome 7 and
contains 142 loci. This regional teaching workflow is separate from the
genome-wide analyses above. Its saved differential-expression outputs have not
been regenerated with the corrected single-end counting settings.

- [Final Project README](Final_Project/README.md)
- [Final Project analysis notebook](Final_Project/RNAseq_Module_I_Final_Project.ipynb)

## Study Notes

- [RNA-seq Data Preparation and STAR Alignment](notes/01_rnaseq_data_preparation_and_star_alignment.md)
- [Cuffdiff2 Differential Expression](notes/02_cuffdiff2_differential_expression.md)
- [DESeq2 Differential Expression](notes/03_deseq2_differential_expression.md)
- [Sex-adjusted DESeq2 and Method Comparison](notes/04_sex_adjusted_deseq2_and_method_comparison.md)

The notes record commands, decisions, errors, and explanations from the exercises. See Note 04 for the gene-level comparisons and their scientific limitations.

## Tools and Technologies

**Programming and environment:** Bash, R, Python, Ubuntu/WSL2, Conda, Git/GitHub.  
**Sequencing and alignment:** NCBI SRA Toolkit, STAR, SAMtools.  
**Transcriptomics:** PsiCLASS, Cuffdiff2, DESeq2, Bioconductor (`GenomicFeatures`, `GenomicAlignments`, `Rsamtools`), `txdbmaker`.  
**Analysis and visualization:** read-count matrices, differential expression statistics, BH/FDR correction, PCA, volcano plots, R/ggplot2.

## Reproducibility

This is a course exercise, not a validated biological study. The current unadjusted analysis and the previously executed sex-adjusted analysis are distinct runs. Rerunning the corrected script 13b is optional and would improve reproducibility; it is not required for the exploratory comparisons documented in Note 04. The unbalanced sex distribution and small sample size remain important limitations.

## Acknowledgements

This repository is based on [Bioinformatics Methods for Transcriptomics](https://www.coursera.org/learn/bioinformatics-methods-for-transcriptomics), taught by Dr. Liliana Florea at Johns Hopkins University via Coursera. Course materials are used as learning references, with local adaptations to scripts and data organization. This is an independent study repository and is not an official university or Coursera resource.

## Back to the Repository

[Main repository README](../../README.md)
