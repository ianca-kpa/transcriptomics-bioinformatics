# Module I — RNA-seq Final Project

This project was developed while following Bioinformatics Methods for Transcriptomics, taught by Liliana Florea at Johns Hopkins University
through Coursera.

I adapted the workflow to Ubuntu/WSL2 and added checks for input files, sample metadata, intermediate outputs, and result summaries.

## Project scope

The analysis uses six course-provided samples: three from 1-month-old mice and three from 4-month-old mice. The reference covers approximately
3.4 Mb of mouse chromosome 7, using GRCm39 and the corresponding GENCODE M30 annotation.

This regional exercise is separate from the whole-genome analysis described in the [Module I README](../README.md).

## Workflow

The [analysis notebook](RNAseq_Module_I_Final_Project.ipynb) includes:

1. Input and sample metadata validation.
2. STAR reference indexing and read alignment.
3. Alignment inspection with STAR logs and SAMtools.
4. Transcript assembly with PsiCLASS.
5. Annotation using GENCODE.
6. Differential expression analysis with Cuffdiff2.
7. Gene-level counting and differential expression analysis with DESeq2.
8. Result summaries.

## Input files

The notebook expects the following files in `data_raw/`:

- `region.fa`
- `gencode.vM30.region.gtf`
- `metadata.tsv`
- `1mo_Rep1.fastq.gz`
- `1mo_Rep2.fastq.gz`
- `1mo_Rep3.fastq.gz`
- `4mo_Rep1.fastq.gz`
- `4mo_Rep2.fastq.gz`
- `4mo_Rep3.fastq.gz`

Obtain the prepared input bundle from the course resources.
Raw reads and reference files are not included in this repository.

The original study's SRA accessions are not a direct replacement for
the prepared regional dataset.

A small sample table is retained in [metadata.tsv](metadata.tsv).
It reproduces the sample names and age groups displayed in the notebook. Compare it with the metadata supplied with the course bundle before use.

## Running the notebook

1. Install the dependencies described in
   [Environment and execution](../ENVIRONMENT.md).
2. Set `REPO_DIR` and `DATA_DIR` in the notebook configuration.
3. Check the paths to the PsiCLASS and Rscript executables.
4. Place the input files in `data_raw/`.
5. Review the read-layout issue below before accepting counting results.
6. Execute the notebook cells in order.

The notebook reuses existing output files. If inputs or parameters change, use a fresh output directory or deliberately regenerate the
affected steps. Existing files do not necessarily correspond to the current configuration.

## Verification still required

The published notebook describes paired-end data, supplies one FASTQ per sample to STAR, and uses `singleEnd = FALSE` with `fragments = TRUE`
for counting.

The original input layout and BAM pairing flags need to be checked locally so that the description, alignment, and counting settings agree.
The settings should not be changed simply to reproduce expected counts.

## Recorded results

The following values come from the saved notebook outputs. They have not been regenerated during the repository documentation update.

| Metric | Recorded result |
|---|---:|
| PsiCLASS consensus loci | 168 |
| PsiCLASS consensus transcripts | 211 |
| Loci in the DESeq2 result table | 142 |
| DESeq2 loci with unadjusted p ≤ 0.05 | 6 |
| DESeq2 loci with adjusted p ≤ 0.05 | 3 |
| Cuffdiff2 FDR-significant loci | 4 |

PsiCLASS loci and distinct reference gene names are different counts. These results also should not be combined with those from the
whole-genome exercise.

Interpretation is limited to this regional dataset and the course exercise. The read-layout verification remains outstanding.

## Repository contents

The notebook, documentation, and small metadata table are versioned. Raw data, indexes, BAM files, and complete generated results are
excluded from Git.

Selected figures and small result summaries can be added after verification, with their analysis scope and parameters documented.

## Acknowledgements

Based on Bioinformatics Methods for Transcriptomics, taught by Liliana Florea at Johns Hopkins University through Coursera.

This is an independent learning repository, not an official course or university resource.
