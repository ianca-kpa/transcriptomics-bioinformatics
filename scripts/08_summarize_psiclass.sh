
#!/usr/bin/env bash
# Module I - PsiCLASS assembly statistics

# Input file
GTF="/mnt/d/Bioinformatics/transcriptomics-data/Module_I/results/PsiCLASS/WithGeneNames/psiclass_vote.withStrand.gtf"

# Count unique genes
echo "Number of genes:"
cut -f9 "$GTF" | cut -d ' ' -f2 | sort -u | wc -l  
# -f9 extracts the 9th column (attributes), cut -d ' ' -f2 extracts the gene_id, sort -u gets unique gene_ids, and wc -l counts them.

# Count unique transcripts
echo "Number of transcripts:"
cut -f9 "$GTF" | cut -d ' ' -f4 | sort -u | wc -l 
# -f9 extracts the 9th column (attributes), cut -d ' ' -f4 extracts the transcript_id, sort -u gets unique transcript_ids, and wc -l counts them.

# Count genes with more than one transcript
echo "Genes with multiple transcripts:"
grep -w "transcript" "$GTF" | cut -f9 | cut -d ' ' -f2 | sort | uniq -c | awk '$1 > 1' | wc -l
# -w "transcript" filters lines containing "transcript", cut -f9 extracts the 9th column (attributes), cut -d ' ' -f2 extracts the gene_id, sort and uniq -c counts occurrences of each gene_id, awk '$1 > 1' filters for genes with more than one transcript, and wc -l counts them.

# Count transcripts with more than one exon
echo "Multi-exon transcripts:"
grep -w "exon" "$GTF" | cut -f9 | cut -d ' ' -f4 | sort | uniq -c | awk '$1 > 1' | wc -l
# -w "exon" filters lines containing "exon", cut -f9 extracts the 9th column (attributes), cut -d ' ' -f4 extracts the transcript_id, sort and uniq -c counts occurrences of each transcript_id, awk '$1 > 1' filters for transcripts with more than one exon, and wc -l counts them.

# Count transcripts with exactly one exon
echo "Single-exon transcripts:"
grep -w "exon" "$GTF" | cut -f9 | cut -d ' ' -f4 | sort | uniq -c | awk '$1 == 1' | wc -l
# -w "exon" filters lines containing "exon", cut -f9 extracts the 9th column (attributes), cut -d ' ' -f4 extracts the transcript_id, sort and uniq -c counts occurrences of each transcript_id, awk '$1 == 1' filters for transcripts with exactly one exon, and wc -l counts them.
