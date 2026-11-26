# Module I — RNA-seq Final Project

This project was developed while following Bioinformatics Methods for Transcriptomics, taught by Liliana Florea at Johns Hopkins University
through Coursera.

I adapted the workflow to Ubuntu/WSL2 and added checks for input files, sample metadata, intermediate outputs, and result summaries.

## Project scope

The analysis uses six course-prepared FASTQ samples: three from 1-month-old mice and three from 4-month-old mice. The reference covers approximately
3.4 Mb of mouse chromosome 7, using GRCm39 and the corresponding GENCODE M30 annotation.

This is a regional teaching dataset supplied for the course exercise. It is separate from the six whole-genome SRA runs used in the main Module I analysis.

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

The notebook expects the following files in `data_raw/`, which are included in this repository:

- `region.fa`
- `gencode.vM30.region.gtf`
- `metadata.tsv`
- `1mo_Rep1.fastq.gz`
- `1mo_Rep2.fastq.gz`
- `1mo_Rep3.fastq.gz`
- `4mo_Rep1.fastq.gz`
- `4mo_Rep2.fastq.gz`
- `4mo_Rep3.fastq.gz`

These course-prepared regional inputs are retained with redistribution permission confirmed.

The original study's SRA accessions are not a direct replacement for
the prepared regional dataset.

A copy of the sample table is also retained in [metadata.tsv](metadata.tsv).
Both tables list the sample names and age groups displayed in the notebook.

## Running the notebook

1. Install the dependencies described in
   [Environment and execution](../ENVIRONMENT.md).
2. The notebook locates this project from its current working directory. Optionally set `TRANSCRIPTOMICS_DATA_DIR` to choose where analysis outputs are stored; by default they go under `results-proj1/`.
3. Check the paths to the PsiCLASS and Rscript executables.
4. Confirm the listed input files are present in `data_raw/`.
5. Execute the notebook cells in order. Corrected single-end count outputs are written under `results/DESeq2_single_end/` so prior count results are not reused.

The notebook reuses existing output files. If inputs or parameters change, use a fresh output directory or deliberately regenerate the
affected steps. Existing files do not necessarily correspond to the current configuration.

## Read layout

The six included FASTQ files contain one read per record, and the corresponding source SRA records are marked `SINGLE`. The notebook therefore supplies one FASTQ per sample to STAR and uses `singleEnd = TRUE` with `fragments = FALSE` for counting.

If replacing the supplied FASTQ files, verify their layout and update the STAR alignment and counting settings together.

## Recorded results

The following values come from saved notebook outputs and have not been regenerated with the corrected single-end counting settings. In particular, regenerate the count matrix and DESeq2 analysis before interpreting those DESeq2 values.

| Metric | Recorded result |
|---|---:|
| PsiCLASS consensus loci | 168 |
| PsiCLASS consensus transcripts | 211 |
| Loci in the DESeq2 result table | 142 |
| DESeq2 loci with unadjusted p ≤ 0.05 | 6 |
| DESeq2 loci with adjusted p ≤ 0.05 | 3 |
| Cuffdiff2 FDR-significant loci | 4 |

PsiCLASS loci and distinct reference gene names are different counts. These results are limited to the regional course dataset and should not be combined with those from the whole-genome exercise.

## Repository contents

The notebook, documentation, sample metadata, and the six prepared regional FASTQs and reference files listed above are versioned.
Indexes, BAM files, and generated analysis results are excluded from Git.

Selected figures and small result summaries can be added after verification, with their analysis scope and parameters documented.

## Acknowledgements

Based on Bioinformatics Methods for Transcriptomics, taught by Liliana Florea at Johns Hopkins University through Coursera.

This is an independent learning repository, not an official course or university resource.
