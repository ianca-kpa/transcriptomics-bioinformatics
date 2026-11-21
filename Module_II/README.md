
# Module II — Alternative Splicing Analysis

Part of the Command Line Tools for Transcriptomics course.

## Overview

This module focuses on identifying and analyzing alternative
splicing patterns from short-read RNA-seq data.

## Dataset

- GEO accession: GSE83931
- Organism: Mus musculus
- Tissue: Hippocampus
- Conditions: 1, 2 and 4 months of development
- Main comparison: 1 month vs. 4 months

## Analysis Workflow

1. RNA-seq data preparation and STAR alignment
2. Isoform-level analysis with Cuffdiff2
3. Event-level analysis with rMATS
4. Intron-level analysis with MntJULiP
5. Intron retention analysis with IRFinder
6. Visualization using IGV and Jutils

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

## Status

In progress.
