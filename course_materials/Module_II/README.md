
# Module II — Alternative Splicing Analysis

Part of the Command Line Tools for Transcriptomics course.

## Overview

This module focuses on identifying and analyzing alternative
splicing patterns from short-read RNA-seq data.

## Dataset

- GEO accession: GSE83931
- Organism: Mus musculus
- Tissue: Hippocampus
- Conditions: 1, 2, and 4 months of development
- Initial rMATS comparison: 3 samples at 1 month versus 3 samples at 4 months

The 1-month and 4-month groups are not balanced by biological sex. Three 1-month samples are female, while the 4-month group contains one female and two males.

## Analysis Workflow

1. RNA-seq data preparation and STAR alignment
2. BAM integrity and read-characteristic validation
3. Isoform-level analysis with Cuffdiff2
4. Event-level alternative splicing analysis with rMATS
5. Intron-level analysis with MntJULiP
6. Intron retention analysis with IRFinder
7. Visualization using IGV and Jutils

### Initial rMATS Results

The junction-count skipped-exon analysis (`SE.MATS.JC.txt`) reported 20,588 events. Using `FDR < 0.05` and `|IncLevelDifference| >= 0.05`, the analysis identified 111 significant events:

| Result | Events |
|---|---:|
| Higher exon inclusion in the 1-month group | 56 |
| Higher exon inclusion in the 4-month group | 55 |

Here, `IncLevelDifference` is the mean inclusion level in group 1 minus the mean inclusion level in group 2. These are alternative splicing events, not counts of distinct genes.

The significant event table is exported to `SE.MATS.JC.significant_FDR0.05_dPSI0.05.tsv` by `scripts/19_summarize_rmats.sh`.

## Repository Organization

- `metadata/`: Sample information and selected run IDs
- `notes/`: Study notes and analysis explanations
- `docs/`: Local course reference materials
- `../scripts/`: Analysis scripts shared across modules

## Data Storage

Raw sequencing data, reference genomes, STAR indexes,
BAM files and generated results are stored outside
the Git repository.

Local paths are configured in `../configs/paths.sh`.
See `../configs/paths.example.sh` for an example.

## Limitations

The rMATS analysis is exploratory. Sequencing layout, read lengths, sex imbalance, and event-level replication should be reviewed before biological interpretation.

## Status

In progress. Initial rMATS analysis completed; downstream event and gene-level review remains to be performed.
