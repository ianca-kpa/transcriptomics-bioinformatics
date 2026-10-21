# Module I — RNA-seq Data Preparation and STAR Alignment

**Course:** Bioinformatics Methods for Transcriptomics — Johns Hopkins University (Coursera)  
**Project:** [transcriptomics-bioinformatics](https://github.com/ianca-kpa/transcriptomics-bioinformatics)  
**Status:** In progress — STAR index completed; first sample aligned and assessed; remaining five alignments pending  
**Environment:** Ubuntu (WSL 2), Bash, NCBI SRA Toolkit, STAR, SAMtools

> **Scope:** This notebook records commands, decisions, and observed results from my practical exercises. It is a learning record rather than a validated, end-to-end RNA-seq pipeline. Some commands are simplified examples; paths and existing files should be checked before rerunning them.

## 1. Dataset and learning objectives

The exercise uses publicly available mouse hippocampus RNA-seq data associated with **GEO accession GSE83931**. I started from `SraRunTable.csv`, containing 30 SRA run accessions, and selected six runs for the tutorial:

| Group in tutorial | SRA run | Library layout |
|---|---|---|
| 1 month | `SRR3734796` | Single-end |
| 1 month | `SRR3734797` | Single-end |
| 1 month | `SRR3734798` | Single-end |
| 4 months | `SRR3734820` | Single-end |
| 4 months | `SRR3734824` | Single-end |
| 4 months | `SRR3734825` | Single-end |

**Learning objectives:** retrieve sequencing data; convert SRA to FASTQ; organize a reproducible project; prepare a mouse reference genome; build a STAR genome index under local memory constraints; align reads; and inspect BAM files and alignment reports.

**Design consideration:** The selected age groups differ in sex composition, so sex could confound any downstream comparison between age groups. No differential expression analysis has been performed in this notebook.

## 2. Software and computing environment

My local setup uses **Ubuntu under Windows Subsystem for Linux (WSL 2)**. The Conda environment is called `transcriptomics`.

```bash
conda activate transcriptomics

STAR --version
samtools --version
free -h
nproc
df -h /mnt/d
```

Recorded software versions:

| Tool | Version |
|---|---|
| STAR | `2.7.11b` |
| SAMtools | `1.22.1` |

The computer has approximately 16 GB of physical RAM. For this exercise, I configured WSL with roughly **9.7 GB of RAM, 8 GB of swap, and 8 available logical processors** by creating `C:\Users\ianca\.wslconfig`:

```ini
[wsl2]
memory=10GB
processors=8
swap=8GB
```

The configuration was applied by running `wsl --shutdown` from Windows Command Prompt and reopening Ubuntu. Increasing swap does not replace physical RAM; it can substantially slow memory-intensive analyses.

## 3. Project layout and data storage

The Git repository is stored on the Windows C: drive. Large sequencing inputs, the genome reference, the STAR index, and alignment outputs are stored on the D: drive.

**Code and metadata (repository):**

```text
transcriptomics-bioinformatics/
├── README.md
├── .gitignore
└── Module_I/
    ├── scripts/
    │   ├── 01_extract_sra_accessions.sh
    │   ├── 02_download_sra.sh
    │   ├── 03_convert_sra_to_fastq.sh
    │   ├── 04_generate_star_index.sh
    │   ├── 05a_test_star_alignment.sh
    │   ├── 05b_align_all_samples_star.sh
    │   └── tutorial_sra_download_and_conversion.sh
    ├── metadata/
    │   ├── SraRunTable.csv
    │   ├── Table.txt
    │   ├── run_id.txt
    │   └── selected_runs.txt
    ├── notes/
    │   └── 01_rnaseq_data_preparation_and_star_alignment.md
    └── data/
        ├── SRA   -> /mnt/d/.../Module_I/SRA
        └── FASTQ -> /mnt/d/.../Module_I/FASTQ
```

**Large files (D: drive):**

```text
/mnt/d/Bioinformatics/transcriptomics-data/Module_I/
├── SRA/
├── FASTQ/
├── reference/
│   ├── genome/
│   └── star_index/
└── results/
    └── STAR/
```

To check the links from `Module_I`:

```bash
ls -l data/
readlink -f data/SRA
readlink -f data/FASTQ
```

The large files are not intended for Git version control.

## 4. Extracting SRA run accessions

The current `scripts/01_extract_sra_accessions.sh` finds the project directory relative to the script's own location, so it does not depend on the shell's working directory. It reads the first column of `SraRunTable.csv`, skipping the header.

**Script — `scripts/01_extract_sra_accessions.sh`:**

```bash
#!/usr/bin/env bash
set -euo pipefail

# Get the absolute directory containing this script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Get the project root directory
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

# Define paths
METADATA_DIR="$PROJECT_DIR/metadata"
INPUT_FILE="$METADATA_DIR/SraRunTable.csv"
OUTPUT_FILE="$METADATA_DIR/run_id.txt"

# Extract run IDs from the first column, skipping the header
tail -n +2 "$INPUT_FILE" | cut -d',' -f1 > "$OUTPUT_FILE"

echo "Processing completed!"
echo "Run IDs saved to: $OUTPUT_FILE"
```

Run it from `Module_I`:

```bash
bash scripts/01_extract_sra_accessions.sh
wc -l metadata/run_id.txt
head metadata/run_id.txt
```

**Observed:** 30 SRA run accessions. Six runs were selected and stored in `metadata/selected_runs.txt`. The selection file is an input to scripts 02 and 03; the extraction script does not itself create that selection file.

## 5. Downloading sequencing runs with SRA Toolkit

The `prefetch` utility downloads runs from the NCBI Sequence Read Archive. The current download script reads the six selected accession numbers, resolves project paths relative to its own location, and skips empty/commented lines.

**Script — `scripts/02_download_sra.sh`:**

```bash
#!/usr/bin/env bash
set -euo pipefail

# Get the absolute directory containing this script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Get the project root directory
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

# Define paths
SRA_DIR="$PROJECT_DIR/data/SRA"
RUN_LIST="$PROJECT_DIR/metadata/selected_runs.txt"

# Create output directory
mkdir -p "$SRA_DIR"

# Download each SRA run
while IFS= read -r run || [[ -n "$run" ]]; do

    # Skip empty lines and comments
    [[ -z "$run" || "$run" == \#* ]] && continue
    echo "Downloading $run..."

    prefetch "$run" -O "$SRA_DIR"
    echo "Downloaded: $run"

done < "$RUN_LIST"
echo "All downloads completed!"
```

Execute from `Module_I`:

```bash
bash scripts/02_download_sra.sh
find -L data/SRA -type f -name '*.sra'
```

**Path and reproducibility note:** This version uses `prefetch "$run" -O "$SRA_DIR"`. The corresponding conversion script expects a file at `data/SRA/<RUN>/<RUN>.sra`. After downloading, always inspect the actual directory structure before conversion, especially if changing `prefetch` options or SRA Toolkit versions. The scripts were separated after path-related problems in earlier attempts.

## 6. Converting SRA to FASTQ

The `fasterq-dump` utility converts SRA files into reads and their Phred quality scores. The current script explicitly sets the FASTQ filename using `--outfile`, avoiding the `.sra.fastq` suffix encountered in an earlier version.

**Script — `scripts/03_convert_sra_to_fastq.sh`:**

```bash
#!/usr/bin/env bash
set -euo pipefail

# Get the absolute directory containing this script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Get the project root directory
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

# Define paths
SRA_DIR="$PROJECT_DIR/data/SRA"
FASTQ_DIR="$PROJECT_DIR/data/FASTQ"
RUN_LIST="$PROJECT_DIR/metadata/selected_runs.txt"

# Create output directory
mkdir -p "$FASTQ_DIR"

# Convert each SRA run to FASTQ
while IFS= read -r run || [[ -n "$run" ]]; do

    # Skip empty lines and comments
    [[ -z "$run" || "$run" == \#* ]] && continue
    echo "Converting $run..."

    # Locate the downloaded SRA file
    SRA_FILE="$SRA_DIR/$run/$run.sra"

    if [[ ! -f "$SRA_FILE" ]]; then
        echo "ERROR: SRA file not found: $SRA_FILE" >&2
        exit 1
    fi

    fasterq-dump "$SRA_FILE" \
        --outdir "$FASTQ_DIR" \
        --outfile "$run.fastq" \
        --threads 4
    echo "Converted: $run"

done < "$RUN_LIST"
echo "All conversions completed!"
```

Execute from `Module_I`:

```bash
bash scripts/03_convert_sra_to_fastq.sh
ls -lh data/FASTQ/
```

**Existing FASTQ outputs:**

```text
SRR3734796.fastq
SRR3734797.fastq
SRR3734798.fastq
SRR3734820.fastq
SRR3734824.fastq
SRR3734825.fastq
```

**Observed:** approximately 27 GB of FASTQ data across six *single-end* runs. The command preserves separate data-conversion logic rather than downloading and converting in a single loop. The script should not be rerun over existing outputs without checking available disk space and overwrite behavior.

**Original course-style reference:** The separately uploaded `tutorial_sra_download_and_conversion.sh` uses six explicit `prefetch` commands followed by six `fasterq-dump` calls. The numbered scripts are my adaptation for maintainability and consistent file locations.

## 7. Moving large files to the D: drive

The original SRA and FASTQ files were copied to the D: drive with `rsync` and checked before the local directories were replaced by symbolic links.

```bash
# Commands used from Module_I, before symbolic links were created.
DEST=/mnt/d/Bioinformatics/transcriptomics-data/Module_I

rsync -avh --progress data/SRA/ "$DEST/SRA/"
rsync -avh --progress data/FASTQ/ "$DEST/FASTQ/"
```

Checks included comparing directory contents and calculating file hashes with `md5sum`. The first transfer required an additional synchronization; the resulting SRA and FASTQ directory sizes matched between source and destination. Symbolic links were then created:

```bash
# Historical illustration only: source directories were moved out
# of these paths first. Do not rerun over existing links.
ln -s "$DEST/SRA" data/SRA
ln -s "$DEST/FASTQ" data/FASTQ
```

**Important:** matching total directory sizes alone does not prove that files are identical. File-by-file hashes or equivalent integrity checks are needed before deleting original data.

## 8. Preparing the mouse reference genome

The tutorial's STAR index was built from an mm39 FASTA referred to as `mm39c.fa`, but that specific file was not provided. I obtained the public **GENCODE M30 GRCm39 primary assembly** FASTA instead. It uses chromosome labels such as `chr1`, `chr2`, etc.

**Reference:** [GENCODE mouse release M30](https://www.gencodegenes.org/mouse/release_M30.html)

```bash
REF_DIR=/mnt/d/Bioinformatics/transcriptomics-data/Module_I/reference/genome
mkdir -p "$REF_DIR"
cd "$REF_DIR"

wget -c \
  https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M30/GRCm39.primary_assembly.genome.fa.gz

# Test compressed-file integrity.
gzip -t GRCm39.primary_assembly.genome.fa.gz

# View FASTA sequence names.
zcat GRCm39.primary_assembly.genome.fa.gz | grep '^>' | head

# Decompress while retaining the original archive.
gunzip -k GRCm39.primary_assembly.genome.fa.gz
```

Observed files:

| File | Approximate disk size |
|---|---:|
| `GRCm39.primary_assembly.genome.fa.gz` | 739 MiB |
| `GRCm39.primary_assembly.genome.fa` | 2.6 GiB |

**Important distinction:** this public GENCODE FASTA may not be byte-for-byte identical to the instructor's `mm39c.fa`.

## 9. Building a memory-reduced STAR genome index

The instructor's script used STAR `genomeGenerate` with the mm39 FASTA and no GTF. To work within a local WSL environment with approximately 9.7 GiB of RAM and 8 GiB of swap, the parameters `--genomeSAsparseD 16` and `--genomeSAindexNbases 12` were used experimentally.

**Script - `scripts/04_generate_star_index.sh`:**

```bash
#!/usr/bin/env bash
set -e

## STAR genome index generation

DATADIR=/mnt/d/Bioinformatics/transcriptomics-data/Module_I/reference

STAR --runThreadN 2 \
     --runMode genomeGenerate \
     --genomeDir ${DATADIR}/star_index \
     --genomeFastaFiles ${DATADIR}/genome/GRCm39.primary_assembly.genome.fa \
     --genomeSAsparseD 16 \
     --genomeSAindexNbases 12
echo "STAR genome index completed!"
```

**Before the original run**, the target directory was created separately:

```bash
mkdir -p /mnt/d/Bioinformatics/transcriptomics-data/Module_I/reference/star_index
```

Then, from `Module_I`:

```bash
bash scripts/04_generate_star_index.sh
```

**Observed STAR output:**

```text
Oct 20 19:26:19 ... finished generating suffix array
Oct 20 19:28:00 ... completed Suffix Array index
Oct 20 19:28:50 ..... finished successfully
DONE: Genome generation, EXITING
```

**Before committing:** The current script assumes the index directory already exists; adding `mkdir -p "${DATADIR}/star_index"` makes it self-contained. Do not regenerate the already completed index merely to test that change. No GTF was supplied, so annotated splice-junction counts from an index database are expected to be zero.

## 10. STAR alignment: single-sample test and planned full run

The instructor ran three `.fastq.gz` samples concurrently at a time. I adapted the tutorial for uncompressed local FASTQ inputs, using two threads and processing samples sequentially.

### 10.1 First-sample validation — completed

**Script content supplied for `scripts/05a_test_star_alignment.sh` (renamed from `05-1_align_reads_star.sh`; shebang preserved as submitted):**

```bash
#!/usr/bin/env bash

## STAR alignment script

STAR=STAR

BASEDIR=/mnt/d/Bioinformatics/transcriptomics-data/Module_I
DATADIR=${BASEDIR}/FASTQ
WORKDIR=${BASEDIR}/results/STAR
STAR2IDX=${BASEDIR}/reference/star_index

mkdir -p $WORKDIR

for i in SRR3734796
do
  echo "Aligning ${i}"

  time \
  $STAR --runThreadN 2 --genomeDir $STAR2IDX/ \
        --readFilesIn $DATADIR/${i}.fastq \
        --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif \
        --outFileNamePrefix ${WORKDIR}/${i}_

done

echo "All alignments completed!"
```

Command executed from `Module_I`:

```bash
bash scripts/05a_test_star_alignment.sh
```

**Observed STAR output:**

```text
Oct 20 19:31:57 ..... started STAR run
Oct 20 19:31:57 ..... loading genome
Oct 20 19:32:35 ..... started mapping
Oct 20 19:59:17 ..... finished mapping
Oct 20 20:00:21 ..... started sorting BAM
Oct 20 20:02:16 ..... finished successfully

real    30m19.112s
```

The sample `SRR3734796` was successfully aligned, producing a coordinate-sorted BAM file. The current `05a` script does not have `set -e`; adding `set -euo pipefail` would prevent it from printing a success message after a failed STAR command, provided the rest of the script is checked for strict-mode compatibility.

### 10.2 All six samples — prepared, not yet recorded as completed

**Script content supplied for `scripts/05b_align_all_samples_star.sh` (renamed from `05_align_reads_star.sh`; shebang preserved as submitted):**

```bash
#!/usr/bin/env bash
set -euo pipefail

## STAR alignment script
STAR=STAR

BASEDIR=/mnt/d/Bioinformatics/transcriptomics-data/Module_I
DATADIR=${BASEDIR}/FASTQ
WORKDIR=${BASEDIR}/results/STAR
STAR2IDX=${BASEDIR}/reference/star_index

mkdir -p $WORKDIR

for i in SRR3734796 SRR3734797 SRR3734798 SRR3734820 SRR3734824 SRR3734825
do
  echo "Aligning ${i}"

  time \
  $STAR --runThreadN 2 --genomeDir $STAR2IDX/ \
        --readFilesIn $DATADIR/${i}.fastq \
        --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif \
        --outFileNamePrefix ${WORKDIR}/${i}_
done
echo "All alignments completed!"
```

**Before running:** This loop includes `SRR3734796`, which has already been aligned. Running it unchanged may redo or overwrite the first sample's outputs. To process **only the five remaining samples**, change its loop to:

```bash
for i in SRR3734797 SRR3734798 SRR3734820 SRR3734824 SRR3734825
```

Add `set -euo pipefail` after the shebang to stop on a STAR error, and quote output-directory paths. These changes are **recommendations**, not changes already made to the uploaded file. Unlike the instructor's parallel implementation, both uploaded scripts run STAR synchronously (no `&`/`wait`), so they process one sample at a time.

Neither uploaded alignment script specifies `--limitBAMsortRAM`. The successful run confirms that sorting completed for the first sample under its actual conditions, but memory use should continue to be monitored for the other runs.

## 11. Evaluating alignment quality

STAR produces a final alignment summary for each run:

```bash
BAM_DIR=/mnt/d/Bioinformatics/transcriptomics-data/Module_I/results/STAR
cat "$BAM_DIR/SRR3734796_Log.final.out"
```

**Observed metrics for `SRR3734796`:**

| STAR metric | Value |
|---|---:|
| Number of input reads | 11,851,553 |
| Average input read length | 100 nt |
| Uniquely mapped reads | 10,641,673 |
| Uniquely mapped reads % | **89.79%** |
| Reads mapped to multiple loci | 825,369 (6.96%) |
| Reads mapped to too many loci | 53,176 (0.45%) |
| Unmapped: too short | 254,730 (2.15%) |
| Unmapped: other | 76,605 (0.65%) |
| Mismatch rate per base | 0.27% |
| Number of splices: total | 1,641,050 |
| Number of splices: annotated (sjdb) | 0 |

The high uniquely mapped percentage is encouraging for this first sample. However, a successful run and a high mapping rate do not by themselves establish that all sequencing and downstream analysis quality criteria are met.

The reported splice count represents **splice-junction occurrences in alignments**, not 1.64 million distinct introns. The `Annotated (sjdb) = 0` value is expected because the STAR genome index was generated without a splice-junction database or GTF annotation.

## 12. Inspecting and indexing the BAM with SAMtools

STAR generated files including:

```text
SRR3734796_Aligned.sortedByCoord.out.bam
SRR3734796_Log.out
SRR3734796_Log.progress.out
SRR3734796_Log.final.out
SRR3734796_SJ.out.tab
```

The BAM was approximately **878 MB**. A BAM index (`.bam.bai`) of approximately **2.7 MB** was also generated. The conversation confirms the index file exists; a separate `samtools quickcheck` result has not been shared.

```bash
BAM_DIR=/mnt/d/Bioinformatics/transcriptomics-data/Module_I/results/STAR
BAM="$BAM_DIR/SRR3734796_Aligned.sortedByCoord.out.bam"

# Basic structural check: normally silent when successful.
samtools quickcheck -v "$BAM"

# Inspect reference names and alignment metadata.
samtools view -H "$BAM" | head

# Create a random-access index (already completed for this sample).
samtools index "$BAM"
```

`samtools quickcheck` performs limited structural checks; it is not a full read-by-read validation. `samtools index` generates a companion `.bai` file for efficient retrieval of alignments by genomic region.

## 13. Example: inspecting the *Gsn* gene region

The tutorial later investigates the mouse *Gsn* locus on chromosome 2. Because my GENCODE FASTA uses a `chr` prefix, the corresponding region should be queried with `chr2` in the BAM, not `2`.

**Commands to try next (not yet recorded as successfully executed):**

```bash
BAM_DIR=/mnt/d/Bioinformatics/transcriptomics-data/Module_I/results/STAR
BAM="$BAM_DIR/SRR3734796_Aligned.sortedByCoord.out.bam"
REGION="chr2:35146392-35197904"

# Count alignment records overlapping the region.
samtools view -c "$BAM" "$REGION"

# Extract alignments overlapping the region.
samtools view -b "$BAM" "$REGION" \
    > "$BAM_DIR/Gsn.SRR3734796.bam"

# Index the region-specific BAM.
samtools index "$BAM_DIR/Gsn.SRR3734796.bam"
```

The `-c` option counts **alignment records**, which should not automatically be interpreted as unique molecules or distinct biological reads in every context.

## 14. Progress and next steps

- [x] Extract 30 run accessions from SRA metadata.
- [x] Select six sequencing runs for Module I.
- [x] Download SRA data and convert the six runs to FASTQ.
- [x] Organize large files on the D: drive and create symbolic links.
- [x] Download and decompress the public GRCm39 reference genome.
- [x] Build a reduced-memory STAR index successfully.
- [x] Align `SRR3734796` successfully.
- [x] Review its `Log.final.out` and generate its BAM index.
- [x] Review and commit scripts 04, 05a, and 05b.
- [ ] Run STAR on the other five samples, sequentially (avoid reprocessing the first sample).
- [ ] Compare mapping quality metrics across all six samples.
- [ ] Inspect splice junctions and selected genomic loci.
- [ ] Continue with transcript assembly, quantification, and downstream analyses.

## 15. Key terminology

| Term | Meaning |
|---|---|
| **SRA** | NCBI archive for sequencing reads and associated data |
| **FASTQ** | Sequence records with base-quality scores |
| **Reference genome** | Genome assembly against which reads are aligned |
| **STAR genome index** | Search structures generated from the reference genome for STAR alignment |
| **Spliced alignment** | Alignment that can span exon–exon junctions |
| **BAM** | Binary representation of SAM alignment records |
| **BAI** | Index enabling efficient random access to a BAM |
| **CIGAR** | Compact description of alignment operations |
| **`N` in CIGAR** | Region skipped on the reference, often corresponding to an intron in RNA-seq |
| **Uniquely mapped reads** | Reads mapped to one genomic locus under the aligner's reporting criteria |
| **Multi-mapping reads** | Reads mapped to multiple genomic loci |
| **`Log.final.out`** | STAR alignment and mapping summary |
| **`SJ.out.tab`** | STAR-reported splice junctions |

---

### Summary

This exercise documents an initial RNA-seq workflow on a local Ubuntu/WSL environment. A memory-reduced STAR index for the mouse GRCm39 primary assembly was built successfully, and the first sample (`SRR3734796`) was aligned in approximately 30 minutes with **89.79% uniquely mapped reads**. The remaining five samples and the downstream transcriptomics analyses are still to be completed.

### Acknowledgement

This notebook is based on practical exercises from *Bioinformatics Methods for Transcriptomics*, taught by Liliana Florea, PhD, at Johns Hopkins University and offered through Coursera. The commands and file organization include adaptations made during independent practice. This is not an official course resource.

### Pre-commit review notes (04, 05a, 05b)

- **04:** The index completed successfully. Add an explicit output-directory creation step for future reruns or fresh machines; do **not** regenerate the existing successful index just to check a documentation edit.
- **05a:** Label it as the one-sample validation script. Consider adding strict error handling; it completed successfully for `SRR3734796`.
- **05b:** Label it as the complete six-sample workflow, **not yet executed** as such. Its current loop includes the already processed sample; choose either to retain it as a reproducible full-workflow example or omit the first run when actually processing the remaining five.
- **Paths:** STAR scripts currently use machine-specific absolute paths under `/mnt/d`; document this as a local setup requirement. Scripts 01–03 instead derive their project paths from their script location.
- **Version control:** Commit the Bash scripts and documentation, not the FASTQ, SRA, genome reference, STAR index, BAMs, or alignment log outputs.
