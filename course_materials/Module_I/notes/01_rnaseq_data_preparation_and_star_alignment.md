# Module I — RNA-seq data preparation and STAR alignment

*Course: Bioinformatics Methods for Transcriptomics (Johns Hopkins University, Coursera)*  
*Dataset: GSE83931, mouse hippocampus RNA-seq*

This note covers the six whole-genome SRA runs used in Module I. They are separate from the regional chromosome 7 dataset in `Final_Project/`.

## Samples and experimental design

| Condition | SRA runs | Sex |
|---|---|---|
| 1 month | `SRR3734796`, `SRR3734797`, `SRR3734798` | Female, female, female |
| 4 months | `SRR3734820`, `SRR3734824`, `SRR3734825` | Female, male, male |

The selected records in `course_materials/Module_I/metadata/SraRunTable.csv` report `LibraryLayout=SINGLE` and an average spot length of 100 bases. The converted data therefore has one FASTQ file per run and STAR is given one input file per sample. The sex distribution is unbalanced; age-related results should be interpreted with that limitation in mind.

The selected accession list and full SRA table are in `course_materials/Module_I/metadata/`.

## Script paths and sample-list scope

The preparation scripts now use the repository's Module I and Module II metadata paths and share the external `$BASEDIR/SRA/` and `$BASEDIR/FASTQ/` directories:

- `scripts/01_extract_sra_accessions.sh` reads `course_materials/Module_I/metadata/SraRunTable.csv` and writes `run_id.txt` beside it.
- `scripts/02_download_sra.sh` defaults to the six Module I accessions in `course_materials/Module_I/metadata/selected_runs.txt`. Set `SRA_RUNS_FILE` to select a different list.
- `scripts/03_convert_sra_to_fastq.sh` defaults to the nine Module II accessions in `course_materials/Module_II/metadata/selected_runs_m2.txt`. It also accepts `SRA_RUNS_FILE`; for Module I, point it to the six-accession list.
- `scripts/05a_test_star_alignment.sh` tests one Module I sample. `scripts/05b_align_all_samples_star.sh` aligns the nine Module II runs; it is not a six-sample Module I batch aligner.

The Module II list has nine runs because it includes the 1-, 2-, and 4-month groups. The current rMATS comparison uses the six 1mo and 4mo samples as a subset. Do not use the nine-run list in place of the six Module I accessions.

## Downloading and converting SRA runs

The intended workflow uses `prefetch` to download each accession and `fasterq-dump` to convert it to FASTQ:

```bash
prefetch "$run" -O "$SRA_DIR"

fasterq-dump "$SRA_FILE" \
    --outdir "$FASTQ_DIR" \
    --outfile "$run.fastq" \
    --threads 4
```

Both scripts source `configs/paths.sh`; SRA files are stored under `$BASEDIR/SRA/<run>/<run>.sra`, and `fasterq-dump` writes uncompressed FASTQs under `$BASEDIR/FASTQ/`. To download the Module II list, run:

```bash
SRA_RUNS_FILE=course_materials/Module_II/metadata/selected_runs_m2.txt \
    ./scripts/02_download_sra.sh
```

The conversion script defaults to that same Module II list. To convert the Module I samples instead, set `SRA_RUNS_FILE` to `course_materials/Module_I/metadata/selected_runs.txt`.

Check that each expected input exists before alignment:

```bash
source configs/paths.sh
FASTQ_DIR="${BASEDIR}/FASTQ"

for id in SRR3734796 SRR3734797 SRR3734798 SRR3734820 SRR3734824 SRR3734825; do
    test -s "${FASTQ_DIR}/${id}.fastq" || {
        echo "Missing FASTQ: ${FASTQ_DIR}/${id}.fastq" >&2
        exit 1
    }
done
```

## Building the STAR index

`scripts/04_generate_star_index.sh` builds an index from the GRCm39 primary-assembly FASTA at:

```text
${BASEDIR}/reference/genome/GRCm39.primary_assembly.genome.fa
```

The current script uses two threads, `--genomeSAindexNbases 12`, and `--genomeSAsparseD 16`:

```bash
source configs/paths.sh
FASTQ_DIR="${BASEDIR}/FASTQ"
STAR --runThreadN 2 \
    --runMode genomeGenerate \
    --genomeDir "$STAR_INDEX_DIR" \
    --genomeFastaFiles "${BASEDIR}/reference/genome/GRCm39.primary_assembly.genome.fa" \
    --genomeSAsparseD 16 \
    --genomeSAindexNbases 12
```

The GENCODE M30 annotation is used in later transcript-analysis steps; this index-generation command does not pass a GTF file to STAR.

## Testing STAR alignment

`scripts/05a_test_star_alignment.sh` tests the first 1-month sample, `SRR3734796`. It reads an uncompressed, single-end FASTQ and writes the coordinate-sorted BAM and STAR logs under `${STAR_DIR}/test/`.

The main STAR options are:

```bash
STAR --runThreadN 2 \
    --genomeDir "$STAR_INDEX_DIR" \
    --readFilesIn "${FASTQ_DIR}/SRR3734796.fastq" \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMstrandField intronMotif \
    --outFileNamePrefix "${STAR_DIR}/test/SRR3734796_"
```

- `--runThreadN 2` uses two STAR threads.
- `--readFilesIn` takes one FASTQ because these runs are single-end.
- `--outSAMtype BAM SortedByCoordinate` writes a coordinate-sorted BAM.
- `--outSAMstrandField intronMotif` adds strand information for spliced alignments.
- `--outFileNamePrefix` controls the BAM and log filenames.

For each completed alignment, inspect `Log.final.out` for read counts and mapping rates. Use `samtools quickcheck` on the BAM before downstream analysis. The Module II batch script performs this check, but the Module I test script does not.

## Alignment scope and interpretation

The Module I sample list contains six runs, but the checked-in `05a` script tests only one and `05b` is configured for Module II. A Module I batch-alignment script is still needed for a reproducible six-sample alignment run from the checked-in scripts.

The sex imbalance (all 1-month samples are female, while two 4-month samples are male) is a limitation of the experimental design, not something alignment can correct. See [the DESeq2 note](03_deseq2_differential_expression.md) for the downstream model and interpretation.

*Commands and paths are summarized from the repository scripts and local course metadata; verify `configs/paths.sh` and input files on the machine used for analysis.*
