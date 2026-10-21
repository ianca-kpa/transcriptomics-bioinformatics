# Transcriptomics & Bioinformatics

**RNA-seq Data Analysis | STAR | SAMtools | Bash | Reproducible Scientific Workflows**

Practical exercises, scripts, and notes developed while following the **Bioinformatics Methods for Transcriptomics** course offered by Johns Hopkins University through Coursera.

**Course:** [Bioinformatics Methods for Transcriptomics](https://www.coursera.org/learn/bioinformatics-methods-for-transcriptomics)  
**Institution:** Johns Hopkins University  
**Instructor:** Liliana Florea, PhD  
**Status:** In progress

## About This Repository

This repository documents my hands-on learning in transcriptomics and RNA-seq data analysis.

My background combines experimental life sciences, bioinformatics, and omics data analysis. I hold a PhD in Tropical Medicine and a Master 2 in Bioinformatics and Computational Biology from Université Côte d'Azur, France.

My previous computational work has focused on untargeted LC–MS/MS metabolomics, statistical analysis, data visualization, and reproducible R-based workflows. Through this project, I am expanding my practical experience with transcriptomic data processing, sequence alignment, and the interpretation of RNA-seq results.

The exercises are based on the course materials, with personal adaptations to improve script organization, automation, data management, and reproducibility.

The repository will be updated progressively as I complete the modules and explore additional computational methods for transcriptomic data analysis.

## Course Topics

The course covers four main areas:

1. **Gene expression analysis:** RNA-seq data preparation, read alignment, transcript assembly, and differential expression analysis.
2. **Alternative splicing analysis:** Computational approaches to identifying and comparing splicing patterns.
3. **Long-read transcriptomics:** Transcriptome reconstruction using long-read sequencing technologies.
4. **Differential expression and splicing with long reads:** Methods for studying gene and transcript expression and splicing differences.

## Current Progress

### Module I — RNA-seq Data Preparation and Alignment

The first module focuses on retrieving, processing, and aligning publicly available RNA-seq data.

The exercises use mouse hippocampus sequencing data from [GEO accession GSE83931](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE83931).

Six single-end sequencing runs were selected for the tutorial:

- SRR3734796
- SRR3734797
- SRR3734798
- SRR3734820
- SRR3734824
- SRR3734825

**Completed activities:**

- Extraction of SRA accession identifiers from metadata using Bash.
- Selection and organization of sequencing runs.
- Download of SRA datasets using `prefetch`.
- Conversion of SRA files into FASTQ using `fasterq-dump`.
- Organization of sequencing data using symbolic links and separate storage directories.
- Download and preparation of the GRCm39/mm39 mouse reference genome.
- Generation of a memory-reduced STAR genome index.
- Successful alignment of the first RNA-seq sample using STAR.
- Inspection of alignment statistics and splice junctions.
- BAM inspection and indexing using SAMtools.

The scripts were adapted from the course exercises to work in a local Ubuntu/WSL environment.

### STAR Genome Indexing

The reference genome used for alignment is **GRCm39/mm39**, obtained from GENCODE (mouse release M30).

Since the analysis was performed on a computer with limited available RAM, I adapted the STAR indexing parameters to reduce memory requirements.

The index was generated using:

```bash
STAR --runThreadN 2 \
     --runMode genomeGenerate \
     --genomeDir ${DATADIR}/star_index \
     --genomeFastaFiles ${DATADIR}/genome/GRCm39.primary_assembly.genome.fa \
     --genomeSAsparseD 16 \
     --genomeSAindexNbases 12
```

The genome index was successfully generated using STAR 2.7.11b.

The index was built from the reference FASTA without incorporating a GTF annotation file.

### Initial STAR Alignment Results

The first alignment was performed using sample **SRR3734796**, producing a coordinate-sorted BAM file.

| Alignment metric | Result |
|---|---:|
| Number of input reads | 11,851,553 |
| Average input read length | 100 bp |
| Uniquely mapped reads | 89.79% |
| Reads mapped to multiple loci | 6.96% |
| Reads unmapped: too short | 2.15% |
| Mismatch rate per base | 0.27% |
| Splice-junction occurrences | 1,641,050 |

The alignment completed successfully, with a high proportion of uniquely mapped reads.

The STAR output was also used to explore splice-junction motifs, including canonical GT/AG junctions, and to examine how RNA-seq reads can span exon–exon boundaries.

These statistics correspond to the first processed sample. The remaining five alignments are part of the ongoing work.

**Alignment report:** [SRR3734796_Log.final.out](Module_I/results_/SRR3734796_Log.final.out)

## Repository Structure

```text
transcriptomics-bioinformatics/
├── Module_I/
│   ├── metadata/
│   │   ├── SraRunTable.csv
│   │   ├── Table.txt
│   │   ├── run_id.txt
│   │   └── selected_runs.txt
│   ├── scripts/
│   │   ├── 01_extract_sra_accessions.sh
│   │   ├── 02_download_sra.sh
│   │   ├── 03_convert_sra_to_fastq.sh
│   │   ├── 04_generate_star_index.sh
│   │   ├── 05a_test_star_alignment.sh
│   │   ├── 05b_align_all_samples_star.sh
│   │   └── tutorial_sra_download_and_conversion.sh
│   ├── notes/
│   │   └── 01_rnaseq_data_preparation_and_star_alignment.md
│   └── results_/
│       └── SRR3734796_Log.final.out
├── .gitignore
└── README.md
```

Large sequencing datasets, reference genome files, STAR indexes, and BAM files are stored locally and excluded from version control.

The project uses symbolic links to access sequencing data and analysis results stored on a separate drive, keeping the Git repository focused on scripts, metadata, documentation, and selected small output files.

## Scripts and Workflow

The scripts are organized according to the main processing steps.

| Script | Description |
|---|---|
| `01_extract_sra_accessions.sh` | Extract SRA run identifiers from sequencing metadata |
| `02_download_sra.sh` | Download selected sequencing datasets |
| `03_convert_sra_to_fastq.sh` | Convert SRA datasets to FASTQ |
| `04_generate_star_index.sh` | Generate a memory-reduced STAR reference genome index |
| `05a_test_star_alignment.sh` | Test STAR alignment using a single sample |
| `05b_align_all_samples_star.sh` | Run STAR alignment for the six selected samples |

The first three scripts resolve input and output paths relative to their location. The STAR scripts currently use absolute paths configured for the local WSL environment and require path adjustments when executed on another system.

## Learning Notebook

A more detailed record of the exercises is available in the Module I notebook:

**[RNA-seq Data Preparation and STAR Alignment](Module_I/notes/01_rnaseq_data_preparation_and_star_alignment.md)**

The notebook includes:

- Bash commands and script explanations.
- Dataset preparation and FASTQ conversion.
- Data organization and symbolic links.
- Reference genome preparation.
- STAR genome indexing and memory considerations.
- Read alignment and BAM file inspection.
- Mapping statistics and splice-junction interpretation.

The notebook distinguishes completed analyses from examples and planned steps.

## Tools and Technologies

**Programming and workflow**

- Bash scripting
- Linux / Ubuntu (WSL 2)
- Conda
- Git and GitHub

**Bioinformatics**

- NCBI SRA Toolkit (`prefetch`, `fasterq-dump`)
- STAR 2.7.11b
- SAMtools 1.22.1
- FASTQ, SAM, BAM, and BAM indexing
- Splice-aware read alignment
- RNA-seq alignment quality assessment

**Data management and reproducibility**

- Reference genome preparation
- Metadata processing
- File integrity checks
- Symbolic links
- Structured scripts and documentation

## Development Status

This repository is a work in progress.

The current exercises demonstrate practical experience in RNA-seq data retrieval, preprocessing, reference genome indexing, splice-aware alignment, and preliminary alignment quality assessment.

The next steps include completing the alignment of the remaining samples, comparing mapping statistics, and progressing toward transcript assembly, expression quantification, and differential expression analysis.

The repository represents independent practical learning and does not yet constitute a complete, validated RNA-seq analysis pipeline.

## Acknowledgements

The exercises in this repository are based on the course [Bioinformatics Methods for Transcriptomics](https://www.coursera.org/learn/bioinformatics-methods-for-transcriptomics), taught by Dr. Liliana Florea at Johns Hopkins University and available through Coursera.

The scripts include personal adaptations for practice, organization, and reproducibility.

This is an independent learning repository and is not an official Johns Hopkins University or Coursera resource.

## Author

**Ianca Karine Prudencio de Albuquerque**

PhD in Tropical Medicine | Master 2 in Bioinformatics and Computational Biology

[GitHub](https://github.com/ianca-kpa) | [LinkedIn](https://www.linkedin.com/in/ianca-albuquerque/)