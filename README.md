
# Transcriptomics & Bioinformatics

**RNA-seq | Transcript assembly | Differential expression | Alternative splicing | R | Bash | Python**

This repository is my practical learning portfolio in transcriptomics and bioinformatics. It brings together course exercises, analysis scripts, workflow adaptations, verified results, and notes on methodological limitations.

Module I follows [Bioinformatics Methods for Transcriptomics](https://www.coursera.org/learn/bioinformatics-methods-for-transcriptomics), taught by Dr. Liliana Florea at Johns Hopkins University through Coursera. Module II documents work from *Command Line Tools for Transcriptomics*.

**Portfolio status:** In progress

## About the project

I use this repository to document course exercises, analysis scripts, and the adaptations needed to run transcriptomics workflows locally on Ubuntu/WSL2.

The original exercises were designed for a computing cluster. I organized commands into numbered scripts, introduced configurable data paths, and added checks for input files, BAM integrity, read characteristics, and analysis outputs.

My notes document the commands I used, problems I encountered, and how I checked the results. Some runs have incomplete environment provenance, so I distinguish verified results from work that still needs validation or a reproducibility rerun.

## Modules

| Module | Topic | Status |
|---|---|---|
| [Module I](course_materials/Module_I/README.md) | Genome-wide RNA-seq and a separate regional Final Project | Main genome-wide results verified; corrected sex-adjusted rerun optional |
| [Module II](course_materials/Module_II/README.md) | Alternative-splicing analysis with rMATS | Initial skipped-exon analysis completed; validation and event review remain |

## Module I — RNA-seq and Differential Expression

The genome-wide workflow uses six mouse hippocampus RNA-seq runs from [GSE83931](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE83931), comparing 1-month-old and 4-month-old mice. The 1mo group contains three female samples; the 4mo group contains one female and two male samples.

I used:

- SRA data retrieval and FASTQ preparation.
- Read alignment with STAR.
- Transcript assembly and annotation with PsiCLASS.
- Differential expression analysis with Cuffdiff2 and DESeq2.
- Raw-count matrix reuse, normalization, and statistical analysis with DESeq2.
- PCA, volcano plots, and data validation.

### Verified genome-wide results

| Analysis | Model / significance criterion | Verified results |
|---|---|---|
| Cuffdiff2 | `status == "OK"` and `significant == "yes"` | 14,214 gene records; 261 significant genes; 207 significant isoforms |
| DESeq2 unadjusted | `~ condition`; `padj <= 0.05` | 14,214 genes; 249 significant |
| Previously executed DESeq2 sex-adjusted | `~ sex + condition`; `padj <= 0.05` | 14,214 genes; 124 significant |

By gene ID, the two DESeq2 models share 85 significant genes. Cuffdiff2 and the previously executed sex-adjusted DESeq2 share 76. These are exploratory comparisons, not evidence of biological validity.

Script 13 completed the unadjusted analysis by importing an existing raw-count matrix; it did not recount BAMs. The earlier sex-adjusted results are verified, but their execution environment provenance is incomplete. The corrected script 13b has not been rerun and has not produced results in its planned destination.

Current outputs are under `$BASEDIR/Cuffdiff/` and `$BASEDIR/DESeq2_genomewide/unadjusted/`. The previously verified sex-adjusted results remain under `$BASEDIR/results-proj1/DESeq2/Results/`; the corrected script 13b is prepared for `$BASEDIR/DESeq2_genomewide/sex_adjusted/`. See the [Module I README](course_materials/Module_I/README.md) and [study notes](course_materials/Module_I/notes/) for input paths, workflow details, and provenance.

### Module I Regional Final Project

A separate [Jupyter notebook](course_materials/Module_I/Final_Project/RNAseq_Module_I_Final_Project.ipynb) analyzes approximately 3.4 Mb of mouse chromosome 7 and contains 142 loci. This regional workflow is distinct from the genome-wide analyses above.

The saved differential-expression outputs have not been regenerated with the corrected single-end counting settings; I do not present them as newly validated results. See the [Final Project README](course_materials/Module_I/Final_Project/README.md).

## Module II — Alternative Splicing

Module II uses rMATS to explore alternative splicing in the same mouse hippocampus dataset. The Module II STAR workflow covers nine runs across 1-, 2-, and 4-month groups; the current rMATS comparison uses the six 1mo and 4mo samples.

The verified junction-count skipped-exon results are:

| Metric | Result |
|---|---:|
| Skipped-exon events | 20,588 |
| Events with `FDR < 0.05` and `abs(IncLevelDifference) >= 0.05` | 111 |
| Higher inclusion in 1mo | 56 |
| Higher inclusion in 4mo | 55 |

`IncLevelDifference` is the mean inclusion level in group 1 (1mo) minus group 2 (4mo). These are splicing events, not distinct genes.

This analysis remains exploratory. I still need to validate read layout and read lengths across all six samples and review significant events by gene ID and genomic coordinates. The groups are also unbalanced by sex. Further details are in the [Module II README](course_materials/Module_II/README.md) and [rMATS study note](course_materials/Module_II/notes/01_rmats_alternative_splicing.md).

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
│       ├── metadata/
│       └── notes/
├── .gitignore
└── README.md
```

## Local Configuration and Data

Analysis scripts use paths defined in `configs/paths.sh`, which is machine-specific and excluded from version control.

The repository includes a [configuration example](configs/paths.example.sh), which I can adapt to local data storage and software environments.

Whole-genome sequencing files, STAR indexes, BAM files, and generated analysis outputs are stored outside Git.

The Module I final-project folder contains the prepared regional teaching dataset and supporting reference files. Redistribution permission for these course-provided files has been confirmed.

See [Environment and execution](course_materials/Module_I/ENVIRONMENT.md) for software dependencies and environment information. The scripts organize and document workflows, but the available provenance does not establish fully reproducible end-to-end execution for every analysis.

## Tools and Technologies

**Programming and data handling:** R, Bioconductor, Python, Bash; data manipulation and validation.

**Environment and version control:** Linux/WSL2, Conda, Jupyter, Git/GitHub.

**RNA-seq processing:** SRA Toolkit, STAR, SAMtools, PsiCLASS.

**Statistical and event-level analysis:** Cuffdiff2, DESeq2, rMATS; statistical testing and multiple-testing correction.

**Visualization and workflow documentation:** PCA, volcano plots, differential-expression tables, alternative-splicing summaries, and documented workflow checks.

## Limitations and Next Steps

This repository documents practical learning, not a fully validated biological study. In Module I, there are three female samples and no male samples at 1mo, compared with one female and two males at 4mo. Statistical adjustment for sex does not remove this design limitation.

Current priorities include:

- Optionally rerunning corrected script 13b and recording its software environment.
- Validating read layout and read lengths for all six rMATS samples.
- Reviewing significant splicing events by gene ID and genomic coordinates.
- Improving environment documentation and workflow portability.

## Acknowledgements

Module I is based on **Bioinformatics Methods for Transcriptomics**, taught by Dr. Liliana Florea at Johns Hopkins University through Coursera. Module II documents work from **Command Line Tools for Transcriptomics**, as described in its module documentation.

This is an independent learning repository, not an official university or Coursera resource.
