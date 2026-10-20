# Transcriptomics & Bioinformatics

Practical exercises, scripts, and notes developed while following the **Bioinformatics Methods for Transcriptomics** course offered by Johns Hopkins University through Coursera.

**Course:** [Bioinformatics Methods for Transcriptomics](https://www.coursera.org/learn/bioinformatics-methods-for-transcriptomics)  
**Institution:** Johns Hopkins University  
**Instructor:** Liliana Florea, PhD  
**Status:** In progress

## About This Repository

This repository documents my hands-on learning in transcriptomics and RNA-seq data analysis.

The exercises are based on the course materials, with personal adaptations to improve script organization, automation, and reproducibility.

The repository will be updated progressively as I complete the modules and explore different computational methods for transcriptomic data analysis.

## Course Topics

The course covers four main areas:

1. **Gene expression analysis:** RNA-seq data preparation, read alignment, transcript assembly, and differential expression analysis.
2. **Alternative splicing analysis:** computational approaches to identifying and comparing splicing patterns.
3. **Long-read transcriptomics:** transcriptome reconstruction using long-read sequencing technologies.
4. **Differential expression and splicing with long reads:** methods for studying gene and transcript expression and splicing differences.

## Current Progress

### Module I — RNA-seq Data Preparation

Initial exercises focus on retrieving and organizing publicly available sequencing data from the NCBI Sequence Read Archive (SRA).

The current scripts cover:

- Extraction of SRA accession identifiers from metadata.
- Selection of sequencing runs.
- Download of SRA datasets using `prefetch`.
- Conversion of SRA files into FASTQ using `fasterq-dump`.
- Basic Bash scripting and file organization.

The scripts have been adapted to use paths relative to their location, allowing execution from different working directories.

## Repository Structure

```text
transcriptomics-bioinformatics/
├── Module_I/
│   ├── metadata/
│   │   ├── SraRunTable.csv
│   │   ├── Table.txt
│   │   ├── run_id.txt
│   │   └── selected_runs.txt
│   └── scripts/
│       ├── 01_extract_sra_accessions.sh
│       ├── 02_download_sra.sh
│       ├── 03_convert_sra_to_fastq.sh
│       └── tutorial_sra_download_and_conversion.sh
├── .gitignore
└── README.md
```

Raw sequencing data and generated FASTQ files are stored locally and excluded from version control.

## Tools and Technologies

- Bash scripting
- Linux / WSL
- NCBI SRA Toolkit (`prefetch`, `fasterq-dump`)
- Git and GitHub

Additional bioinformatics tools will be documented as they are introduced and used in the exercises.

## Development Status

This repository is a work in progress.

New scripts, analysis steps, and documentation will be added as I advance through the course. The current contents represent learning exercises rather than a complete, validated RNA-seq analysis pipeline.

## Acknowledgements

The exercises in this repository are based on the course [Bioinformatics Methods for Transcriptomics](https://www.coursera.org/learn/bioinformatics-methods-for-transcriptomics), taught by Dr. Liliana Florea at Johns Hopkins University and available through Coursera.

The scripts include personal adaptations for practice, organization, and reproducibility.

This is an independent learning repository and is not an official Johns Hopkins University or Coursera resource.

## Author

**Ianca Karine Prudencio de Albuquerque**

PhD in Tropical Medicine | Master 2 in Bioinformatics and Computational Biology

[GitHub](https://github.com/ianca-kpa) | [LinkedIn](https://www.linkedin.com/in/ianca-albuquerque/)