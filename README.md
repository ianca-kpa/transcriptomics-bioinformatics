
# Transcriptomics & Bioinformatics

**RNA-seq | Transcript assembly | Differential expression | Alternative splicing | R | Bash | Python | Reproducible workflows**

This repository documents my practical work with the [Bioinformatics Methods for Transcriptomics](https://www.coursera.org/learn/bioinformatics-methods-for-transcriptomics) course, taught by Dr. Liliana Florea at Johns Hopkins University through Coursera.

**Course status:** In progress

## About the project

I use this repository to document course exercises, analysis scripts, and the adaptations needed to run transcriptomics workflows locally on Ubuntu/WSL2.

The original exercises were designed for a computing cluster. I reorganized the commands into numbered scripts, introduced configurable data paths, and added checks for input files, BAM integrity, read characteristics, and analysis outputs.

My notes document the commands I used, problems I encountered, and how I investigated the results. This is a learning project, and I distinguish completed analyses from steps that still require verification.

## Modules

| Module | Topic | Status |
|---|---|---|
| [Module I](course_materials/Module_I/README.md) | RNA-seq alignment, transcript assembly, and differential expression | Initial analyses completed; some verification items remain |
| [Module II](course_materials/Module_II/README.md) | Alternative splicing analysis | In progress; initial rMATS analysis completed |

## Module I — RNA-seq and Differential Expression

The whole-genome exercise uses six mouse hippocampus RNA-seq runs from [GSE83931](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE83931), comparing samples from 1-month-old and 4-month-old mice.

The workflow includes:

- SRA data retrieval and FASTQ preparation.
- Read alignment with STAR.
- Transcript assembly and annotation with PsiCLASS.
- Differential expression analysis with Cuffdiff2 and DESeq2.
- Read counting, normalization, and statistical analysis.
- PCA, volcano plots, and exploration of expression filters.

In the initial DESeq2 model (`~ condition`), 271 of 14,293 genes in the results table had adjusted p-values ≤ 0.05.

The comparison has an important limitation: the age groups are not balanced by biological sex. A second DESeq2 script includes sex as a covariate, but its results remain pending verification.

Detailed methods, results, and limitations are documented in the [Module I README](course_materials/Module_I/README.md).

### Module I Final Project

A separate [Jupyter notebook](course_materials/Module_I/Final_Project/RNAseq_Module_I_Final_Project.ipynb) documents an analysis of approximately 3.4 Mb from mouse chromosome 7.

The workflow includes STAR, PsiCLASS, transcript annotation, Cuffdiff2, and DESeq2.

The saved DESeq2 output contains 142 loci, including 3 with adjusted p-values ≤ 0.05.

These regional results are separate from the whole-genome analysis. Read layout and counting consistency still require verification.

## Module II — Alternative Splicing

Module II explores alternative splicing using RNA-seq data from the same mouse hippocampus dataset.

The current rMATS analysis compares three 1-month samples with three 4-month samples.

I added scripts to:

- Validate BAM files before downstream analysis.
- Examine read length distributions and paired-end alignment flags.
- Prepare BAM lists for the selected conditions.
- Run event-level alternative splicing analysis using rMATS.
- Summarize skipped exon events and export significant results.

### Initial rMATS Results

The junction-count analysis (`SE.MATS.JC.txt`) produced:

| Metric | Result |
|---|---:|
| Skipped exon events reported | 20,588 |
| Events with FDR < 0.05 and \|ΔPSI\| ≥ 0.05 | 111 |
| Higher exon inclusion at 1 month | 56 |
| Higher exon inclusion at 4 months | 55 |

The direction is based on `IncLevelDifference`, calculated as the mean inclusion level in group 1 (1 month) minus that in group 2 (4 months).

These are alternative splicing events, not counts of distinct genes.

The analysis remains exploratory. Read-layout validation, full read-length assessment, and the imbalance in biological sex between groups are important considerations before biological interpretation.

Further details are available in the [Module II README](course_materials/Module_II/README.md).

## Repository Structure

```text
transcriptomics-bioinformatics/
├── configs/
│   └── paths.example.sh
├── scripts/
│   ├── 01_extract_sra_accessions.sh
│   ├── ...
│   ├── 06a_run_psiclass.sh
│   ├── ...
│   ├── 17a_check_read_lengths.sh
│   ├── 18_run_rmats.sh
│   └── 19_summarize_rmats.sh
├── course_materials/
│   ├── Module_I/
│   │   ├── README.md
│   │   ├── ENVIRONMENT.md
│   │   ├── metadata/
│   │   ├── notes/
│   │   └── Final_Project/
│   │       ├── README.md
│   │       ├── RNAseq_Module_I_Final_Project.ipynb
│   │       └── data_raw/
│   └── Module_II/
│       ├── README.md
│       └── metadata/
├── .gitignore
└── README.md
```

## Local Configuration and Data

Analysis scripts use paths defined in `configs/paths.sh`, which is machine-specific and excluded from version control.

The repository includes a [configuration example](configs/paths.example.sh). Users can adapt it to their own data storage and software environment.

Whole-genome sequencing files, STAR indexes, BAM files, and generated analysis outputs are stored outside Git.

The Module I final-project folder currently contains a prepared regional teaching dataset and supporting reference files. Redistribution permissions for course-provided files should be checked separately.

See [Environment and execution](course_materials/Module_I/ENVIRONMENT.md) for software dependencies and environment information.

## Tools and Technologies

**Programming:** R, Python, Bash.

**Environment and version control:** Ubuntu/WSL2, Conda, Jupyter, Git/GitHub.

**RNA-seq processing:** SRA Toolkit, STAR, SAMtools, PsiCLASS.

**Statistical analysis:** Cuffdiff2, DESeq2, rMATS, Bioconductor.

**Visualization and results:** PCA, volcano plots, differential expression tables, alternative splicing summaries.

## Limitations and Next Steps

This repository documents practical training rather than a fully validated biological study.

Current priorities include:

- Verifying sequencing layout and read lengths for downstream analyses.
- Reviewing the sex-adjusted DESeq2 comparison.
- Examining significant alternative splicing events at the gene level.
- Improving environment documentation and workflow portability.
- Adding selected figures and reproducible result summaries.

## Acknowledgements

The exercises are based on **Bioinformatics Methods for Transcriptomics**, taught by Dr. Liliana Florea at Johns Hopkins University through Coursera.

This is an independent learning repository, not an official university or Coursera resource.
