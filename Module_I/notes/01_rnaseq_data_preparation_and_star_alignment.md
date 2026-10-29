# Module I — RNA-seq data preparation and STAR alignment

*Course: Bioinformatics Methods for Transcriptomics (Johns Hopkins University, Coursera)*  
*Project: [transcriptomics-bioinformatics](https://github.com/ianca-kpa/transcriptomics-bioinformatics)*

These are my notes from the first part of Module I. I followed the course exercises using Ubuntu in WSL2 and adapted the commands to my laptop. The main challenge was getting STAR's genome index to fit within the memory available to WSL. I kept the full sequencing files on my D: drive and the scripts and notes in the Git repository on C:.

The six samples have since been aligned and used in PsiCLASS and Cuffdiff2. This note focuses on **data preparation and STAR**, rather than documenting every downstream step again.

## Data used

I worked with mouse hippocampus RNA-seq data from **GEO GSE83931**. The SRA metadata contains 30 runs; for the course exercise I selected these six single-end samples:

| Group | SRA runs |
|---|---|
| 1 month | `SRR3734796`, `SRR3734797`, `SRR3734798` |
| 4 months | `SRR3734820`, `SRR3734824`, `SRR3734825` |

The age groups are the comparison used in the tutorial. Their sex composition may differ, so I should not attribute downstream expression differences entirely to age without checking the metadata.

## Working environment and file organization

- Ubuntu on WSL2, with about **10 GB RAM**, **8 GB swap**, and **8 processors** allocated to WSL.
- STAR **2.7.11b** for alignment; SAMtools **1.22.1** for examining BAM files.
- Conda environment: `transcriptomics`.
- Repository: `/mnt/c/Users/ianca/Documentos/Bioinformatics/transcriptomics-bioinformatics/`.
- Large data: `/mnt/d/Bioinformatics/transcriptomics-data/Module_I/`.

I use symbolic links from `Module_I/data/` to the FASTQ and SRA directories on D:. `Module_I/results` also links to the results directory on D:. This keeps the repository small while letting the scripts refer to files through the project structure.

To check where a link points:

```bash
ls -l data/
readlink -f results
```

## 1. From SRA metadata to FASTQ

I started with `metadata/SraRunTable.csv`. To inspect or extract the accession column, the useful Bash commands were:

```bash
# Skip the CSV header, then take its first field.
tail -n +2 metadata/SraRunTable.csv | cut -d',' -f1 > metadata/run_id.txt

wc -l metadata/run_id.txt
head metadata/run_id.txt
```

I obtained **30 run accessions** and saved the six selected IDs in `metadata/selected_runs.txt`.

The scripts for this part are:

```text
scripts/01_extract_sra_accessions.sh
scripts/02_download_sra.sh
scripts/03_convert_sra_to_fastq.sh
```

The main tools were `prefetch` for downloads and `fasterq-dump` for conversion:

```bash
# Example for one accession; check where prefetch saved the .sra file.
prefetch SRR3734796 -O data/SRA/SRR3734796

# Example of conversion from an existing SRA file.
fasterq-dump data/SRA/SRR3734796/SRR3734796.sra \
    --outdir data/FASTQ \
    --threads 4
```

I encountered a filename issue when `.sra` was included in some generated FASTQ names. I renamed those files instead of repeating the conversion. In the end, I had one `.fastq` per run (single-end sequencing), about **27 GB** in total.

I transferred large files to D: and checked the copies before using symbolic links. Keeping these files outside Git is important: FASTQ and BAM files are much larger than the scripts needed to reproduce the steps.

## 2. Reference genome

The course refers to a mouse mm39 FASTA. I downloaded the public **GENCODE mouse M30 / GRCm39 primary assembly** reference rather than assuming I had the instructor's exact file.

```bash
REFDIR=/mnt/d/Bioinformatics/transcriptomics-data/Module_I/reference/genome
mkdir -p "$REFDIR"
cd "$REFDIR"

wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M30/GRCm39.primary_assembly.genome.fa.gz

gzip -t GRCm39.primary_assembly.genome.fa.gz
gunzip -k GRCm39.primary_assembly.genome.fa.gz
```

The uncompressed FASTA was approximately **2.6 GB**. I checked the chromosome names and saw the `chr` prefix (for example, `chr1` and `chr2`). This matters when querying genomic regions later.

## 3. Building the STAR index

This was the most difficult part of the preparation. Generating an index for the mouse reference required more memory than I could comfortably allocate. I used STAR's memory-saving index settings and two threads:

```bash
STAR --runThreadN 2 \
    --runMode genomeGenerate \
    --genomeDir /mnt/d/Bioinformatics/transcriptomics-data/Module_I/reference/star_index \
    --genomeFastaFiles /mnt/d/Bioinformatics/transcriptomics-data/Module_I/reference/genome/GRCm39.primary_assembly.genome.fa \
    --genomeSAsparseD 16 \
    --genomeSAindexNbases 12
```

- `--runMode genomeGenerate`: build a genome index rather than align reads.
- `--genomeDir`: where STAR writes the index.
- `--genomeFastaFiles`: reference sequence used to build it.
- `--genomeSAsparseD 16` and `--genomeSAindexNbases 12`: settings I changed to reduce index memory requirements.

The index completed successfully on **20 October 2026**, after roughly **28 minutes**. WSL used swap during the process, so I would not expect the same performance as on a computing cluster. I did **not** provide a GTF when building this STAR index; consequently, the STAR report showed zero *annotated* splice junctions in the splice-junction database.

The command is kept in `scripts/04_generate_star_index.sh`.

## 4. Aligning reads with STAR

I first tested STAR with `SRR3734796` before running the other samples. The essential alignment command was:

```bash
STAR --runThreadN 2 \
    --genomeDir /mnt/d/Bioinformatics/transcriptomics-data/Module_I/reference/star_index \
    --readFilesIn /mnt/d/Bioinformatics/transcriptomics-data/Module_I/FASTQ/SRR3734796.fastq \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMstrandField intronMotif \
    --outFileNamePrefix /mnt/d/Bioinformatics/transcriptomics-data/Module_I/results/STAR/SRR3734796_
```

This reads an uncompressed FASTQ, aligns it to the indexed genome, and writes a coordinate-sorted BAM. Because the input was not gzipped, I did not need the course's `--readFilesCommand zcat` option. I also avoided launching several memory-heavy STAR jobs simultaneously.

The single-sample test succeeded. I then completed the remaining five alignments. The current BAM files are organized under `results/STAR/BAM/`, with splice-junction tables in `results/STAR/SpliceJunctions/` and logs kept separately. The scripts used for testing and running the six samples are `05a_test_star_alignment.sh` and `05b_align_all_samples_star.sh`.

### Alignment result checked in detail: SRR3734796

From `Log.final.out`:

| Metric | Value |
|---|---:|
| Input reads | 11,851,553 |
| Uniquely mapped reads | 10,641,673 |
| Uniquely mapped reads (%) | **89.79%** |
| Reads mapped to multiple loci (%) | 6.96% |
| Unmapped: too short (%) | 2.15% |
| Mismatch rate per base (%) | 0.27% |
| Splice-junction occurrences | 1,641,050 |

The first alignment took about **30 minutes**. The unique mapping rate was encouraging, but it is only one part of alignment quality. I have not added a comparable quality-control table for all six samples to this note, because those six per-sample metrics have not yet been collected here.

A useful distinction: the splice count in the STAR report refers to occurrences supporting splicing in alignments, not that number of distinct introns.

## 5. Examining BAM files with SAMtools

I used SAMtools to inspect BAM headers, check basic file structure, and create `.bai` indexes:

```bash
BAM=results/STAR/BAM/SRR3734796_Aligned.sortedByCoord.out.bam

/usr/bin/samtools quickcheck -v "$BAM"
/usr/bin/samtools view -H "$BAM" | head
/usr/bin/samtools index "$BAM"
```

A `.bam` contains aligned reads. A `.bam.bai` is a separate index for quickly retrieving alignments at a particular genomic location. `quickcheck` checks basic file structure; it does not inspect every alignment.

One practical issue I ran into later: an older `samtools` executable inside the `psiclass` Conda environment did not recognize `quickcheck`. The working SAMtools **1.22.1** executable was `/usr/bin/samtools`, so I used its full path where necessary.

### Example locus: Gsn

The course also uses the mouse *Gsn* locus. A region can be inspected with:

```bash
/usr/bin/samtools view -c \
    results/STAR/BAM/SRR3734796_Aligned.sortedByCoord.out.bam \
    'chr2:35146392-35197904'
```

This counts alignment records overlapping the region, not necessarily unique biological molecules. I am keeping it as an inspection example, not claiming a particular count without its output.

## Where this led next

With the six BAM files available, I continued to **PsiCLASS** for transcript assembly and then **Cuffdiff2** for the Month 1 versus Month 4 comparison. Those steps have their own scripts and notes, including `notes/02_cuffdiff2_differential_expression.md`.

What I want to improve next is a compact table of STAR quality metrics for **all six** samples, plus a closer look at splice-junction support. The metadata also needs attention before interpreting expression differences as age-related effects.

*Based on the Johns Hopkins University Coursera exercises, with commands adapted for my local WSL2 setup. These are notes from my practice, not an official course tutorial.*
