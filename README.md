# Transcriptomics & Bioinformatics

**RNA-seq | Transcript assembly | Differential expression | R | Bash | Reproducible workflows**

This repository contains scripts and study notes from my practical work with the [Bioinformatics Methods for Transcriptomics](https://www.coursera.org/learn/bioinformatics-methods-for-transcriptomics) course, taught by Dr. Liliana Florea at Johns Hopkins University through Coursera.

**Course status:** In progress

## About the project

I use this repository to keep the scripts and notes from each course exercise. The original commands were designed for a computing cluster, so I adapted them to Ubuntu/WSL2 and organized them into numbered scripts. Sequencing files and generated results are stored separately from the Git repository. My notes include commands I tested, errors I encountered, and changes I made.

## Modules

| Module | Topic | Progress |
|---|---|---|
| [Module I](Module_I/README.md) | RNA-seq preparation, STAR alignment, PsiCLASS assembly, Cuffdiff2 and DESeq2 | Initial analyses completed; sex-adjusted comparison pending verification |
| Later modules | Alternative splicing and long-read transcriptomics | Not yet documented |

## Work completed in Module I

The current analysis uses six mouse hippocampus RNA-seq runs from [GSE83931](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE83931), comparing 1-month and 4-month samples.

- Prepared sequencing reads and aligned all six samples with **STAR**.
- Assembled and annotated transcripts using **PsiCLASS**.
- Performed differential expression analyses using **Cuffdiff2** and **DESeq2**.
- Generated count matrices, PCA and volcano plots, and annotated DESeq2 output tables with gene names.
- Explored the effects of `baseMean`, absolute log2 fold change, and p-value/FDR cutoffs.

In the initial DESeq2 model (`~ condition`), **271 genes** had `padj ≤ 0.05` out of **14,293 genes** in the result table. The six samples are not balanced by sex: the 1-month group contains only females, whereas the 4-month group includes two males. This limits how confidently the differences can be attributed to age. I have written a second DESeq2 script including sex in the design, but have not yet confirmed its output.

For the dataset, methods, exact script list, result tables and limitations, see **[Module I — Workflow and Results](Module_I/README.md)**.

## Repository layout

```text
transcriptomics-bioinformatics/
├── Module_I/
│   ├── README.md       # Module-specific workflow and results
│   ├── scripts/        # Numbered Bash, R and Python scripts
│   ├── metadata/       # Run selection and sample metadata
│   ├── notes/          # Study notes
│   ├── docs/           # Local course reference materials (gitignored)
│   ├── data/           # Local links to sequencing data
│   └── results/        # Local link to analysis outputs
├── .gitignore
└── README.md
```

Large sequencing datasets, indexes, BAMs and generated analysis results are kept outside Git version control. Several scripts contain paths specific to my local WSL2 setup and require adjustment before use elsewhere.

## Tools

Bash, R, Python, Git/GitHub, Ubuntu/WSL2, Conda, SRA Toolkit, STAR, SAMtools, PsiCLASS, Cuffdiff2, DESeq2, and Bioconductor.

## Acknowledgements

The exercises are based on **Bioinformatics Methods for Transcriptomics**, taught by Dr. Liliana Florea at Johns Hopkins University via Coursera. This is an independent learning repository, not an official university or Coursera resource.
