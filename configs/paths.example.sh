
#!/usr/bin/env bash

# External data directory
BASEDIR="/path/to/transcriptomics-data"

# Input and output directories
DATADIR="${BASEDIR}/FASTQ"
WORKDIR="${BASEDIR}/STAR"
STAR2IDX="${BASEDIR}/reference/star_index"

# Executables
STAR="STAR"
