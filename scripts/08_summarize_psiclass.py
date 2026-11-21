
#!/usr/bin/env python3

# ============================================================
# Module I - PsiCLASS assembly statistics
#
# This script reads a GTF file and calculates:
# 1. Number of genes
# 2. Number of transcripts
# 3. Genes with multiple transcripts
# 4. Multi-exon transcripts
# 5. Single-exon transcripts
# ============================================================

from pathlib import Path
from collections import defaultdict
import csv


# ------------------------------------------------------------
# 1. Define input and output paths
# ------------------------------------------------------------

# Directory containing the annotated PsiCLASS GTF file
WORKDIR = Path(
    "/mnt/d/Bioinformatics/transcriptomics-data/Module_I/"
    "results/PsiCLASS/WithGeneNames"
)

# Input GTF file
GTF = WORKDIR / "psiclass_vote.withStrand.gtf"

# Output summary table (tab-separated values)
OUTPUT = WORKDIR / "psiclass_summary.tsv"


# ------------------------------------------------------------
# 2. Create dictionaries to store information
# ------------------------------------------------------------

# Store unique transcript IDs for each gene
# Example: {"gene1": {"tx1", "tx2"}}
genes = defaultdict(set)

# Store the number of exons for each transcript
# Example: {"tx1": 3, "tx2": 1}
exons = defaultdict(int)


# ------------------------------------------------------------
# 3. Read the GTF file line by line
# ------------------------------------------------------------

with GTF.open() as file:

    for line in file:

        # Skip comment lines and empty lines
        if line.startswith("#") or not line.strip():
            continue

        # Split the line into 9 tab-separated columns
        columns = line.rstrip("\n").split("\t")

        # Ignore lines that do not contain 9 columns
        if len(columns) != 9:
            continue

        # Column 3: genomic feature (transcript or exon)
        feature = columns[2]

        # Column 9: gene and transcript attributes
        attributes = columns[8].split(";")

        # Extract gene_id and transcript_id
        # Assumes these are the first two attributes
        gene_id = attributes[0].split('"')[1]
        transcript_id = attributes[1].split('"')[1]

        # Register the gene-transcript relationship
        genes[gene_id].add(transcript_id)

        # Initialize the exon count for this transcript
        # setdefault preserves an existing count
        exons.setdefault(transcript_id, 0)

        # Count each exon associated with the transcript
        if feature == "exon":
            exons[transcript_id] += 1


# ------------------------------------------------------------
# 4. Calculate assembly statistics
# ------------------------------------------------------------

# Number of unique genes
num_genes = len(genes)

# Number of unique transcripts
num_transcripts = len(exons)

# Number of genes associated with more than one transcript
genes_multiple = sum(
    len(transcripts) > 1
    for transcripts in genes.values()
)

# Number of transcripts containing two or more exons
multi_exon = sum(
    count > 1
    for count in exons.values()
)

# Number of transcripts containing exactly one exon
single_exon = sum(
    count == 1
    for count in exons.values()
)


# ------------------------------------------------------------
# 5. Display results in the terminal
# ------------------------------------------------------------

print("\nPsiCLASS assembly statistics")
print("----------------------------")

print("Genes:", num_genes)
print("Transcripts:", num_transcripts)
print("Genes with multiple transcripts:", genes_multiple)
print("Multi-exon transcripts:", multi_exon)
print("Single-exon transcripts:", single_exon)


# ------------------------------------------------------------
# 6. Export results to a TSV file
# ------------------------------------------------------------

# Open the output file in writing mode
with OUTPUT.open("w", newline="") as file:

    # Create a tab-separated table
    writer = csv.writer(file, delimiter="\t")

    # Write column headers
    writer.writerow([
        "Sample Name",
        "NumGenes",
        "NumTranscripts",
        "NumGenes >1 txpt",
        "NumTxpts >1 exon",
        "NumSingleExonTxpts"
    ])

    # Write the calculated statistics
    writer.writerow([
        "psiclass_vote.withStrand",
        num_genes,
        num_transcripts,
        genes_multiple,
        multi_exon,
        single_exon
    ])

print(f"\nSummary saved: {OUTPUT}")
