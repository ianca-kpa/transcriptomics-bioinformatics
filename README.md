# Transcriptomics & Bioinformatics

**RNA-seq | Transcript assembly | Differential expression | R | Bash | Reproducible workflows**

This repository contains scripts, notebooks, and study notes from my practical work with the [Bioinformatics Methods for Transcriptomics](https://www.coursera.org/learn/bioinformatics-methods-for-transcriptomics) course, taught by Dr. Liliana Florea at Johns Hopkins University through Coursera.

**Course status:** In progress

## About the project

I use this repository to document the course exercises and the changes I made while running them locally. The original commands were designed for a computing cluster, so I adapted them to Ubuntu/WSL2 and organized them into numbered scripts.

Sequencing files and large generated results are stored separately from Git. My notes include commands I tested, errors I encountered, and how I checked the outputs.

## Modules

| Module | Topic | Progress |
|---|---|---|
| [Module I](Module_I/README.md) | RNA-seq preparation, STAR alignment, PsiCLASS assembly, Cuffdiff2, and DESeq2 | Initial analyses completed; final-project notebook available; verification items remain |
| Later modules | Alternative splicing and long-read transcriptomics | Not yet documented |

## Work completed in Module I

The whole-genome exercise uses six mouse hippocampus RNA-seq runs from [GSE83931](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE83931), comparing 1-month and 4-month samples.

- Prepared sequencing reads and aligned all six samples with **STAR**.
- Assembled and annotated transcripts using **PsiCLASS**.
- Performed differential expression analyses using **Cuffdiff2** and **DESeq2**.
- Generated count matrices, PCA and volcano plots, and annotated DESeq2 output tables with gene names.
- Explored the effects of `baseMean`, absolute log2 fold change, and p-value/FDR cutoffs.

In the initial DESeq2 model (`~ condition`), **271 genes** had `padj ≤ 0.05` out of **14,293 genes** in the result table.

The six samples are not balanced by sex: the 1-month group contains three females, whereas the 4-month group contains one female and two males. This limits how confidently the differences can be attributed to age. A second DESeq2 script includes sex in the design, but its output remains pending verification.

For the dataset, methods, script list, results, and limitations, see **[Module I — Workflow and Results](Module_I/README.md)**.

## Final Project — Regional Dataset

The [final-project notebook](Module_I/Final_Project/RNAseq_Module_I_Final_Project.ipynb) documents a separate analysis of approximately 3.4 Mb from mouse chromosome 7, using six course-provided samples and the corresponding GENCODE M30 annotation.

The workflow includes input validation, STAR alignment, PsiCLASS assembly, annotation, Cuffdiff2, and DESeq2.

The saved DESeq2 result table contains **142 loci**, including **3 with adjusted p-values ≤ 0.05**. These results describe the regional exercise and should be distinguished from the whole-genome results above. PsiCLASS loci and distinct reference gene names are also different counts.

Read layout and counting consistency still require local verification before the final-project results can be considered verified.

## Repository layout

## Repository layout

```text
transcriptomics-bioinformatics/
├── Module_I/
│   ├── README.md
│   ├── ENVIRONMENT.md
│   ├── Final_Project/
│   │   ├── RNAseq_Module_I_Final_Project.ipynb
│   │   ├── metadata.tsv
│   │   ├── data_raw/    # Local input files (gitignored)
│   │   └── results/     # Local link to outputs (gitignored)
│   ├── scripts/        # Numbered Bash, R, and Python scripts
│   ├── metadata/       # SRA metadata, run selection, and phenotype tables
│   ├── notes/          # Study notes
│   ├── docs/           # Local course materials and utilities (gitignored)
│   └── data/           # Local links to sequencing data (gitignored)
├── .gitignore
└── README.md
```

Large sequencing files, reference indexes, BAM files, and generated analysis outputs are excluded from version control.

Several scripts and notebook settings contain paths specific to my local WSL2 setup and require adjustment before use elsewhere.

See [Environment and execution](Module_I/ENVIRONMENT.md) for dependencies and instructions to record the analysis environment.

## Tools

Bash, R, Python, Jupyter, Git/GitHub, Ubuntu/WSL2, Conda, SRA Toolkit, STAR, SAMtools, PsiCLASS, Cuffdiff2, DESeq2, and Bioconductor.

## Acknowledgements

The exercises are based on **Bioinformatics Methods for Transcriptomics**, taught by Dr. Liliana Florea at Johns Hopkins University through Coursera.

This is an independent learning repository, not an official university or Coursera resource.