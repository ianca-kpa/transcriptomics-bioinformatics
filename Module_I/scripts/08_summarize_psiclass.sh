
#!/usr/bin/env bash
# Module I - PsiCLASS assembly statistics

# Input file
GTF="/mnt/d/Bioinformatics/transcriptomics-data/Module_I/results/PsiCLASS/WithGeneNames/psiclass_vote.withStrand.gtf"

# Count unique genes
echo "Number of genes:"
cut -f9 "$GTF" | cut -d ' ' -f2 | sort -u | wc -l

# Count unique transcripts
echo "Number of transcripts:"
cut -f9 "$GTF" | cut -d ' ' -f4 | sort -u | wc -l

# Count genes with more than one transcript
echo "Genes with multiple transcripts:"
grep -w "transcript" "$GTF" | cut -f9 | cut -d ' ' -f2 | sort | uniq -c | awk '$1 > 1' | wc -l

# Count transcripts with more than one exon
echo "Multi-exon transcripts:"
grep -w "exon" "$GTF" | cut -f9 | cut -d ' ' -f4 | sort | uniq -c | awk '$1 > 1' | wc -l

# Count transcripts with exactly one exon
echo "Single-exon transcripts:"
grep -w "exon" "$GTF" | cut -f9 | cut -d ' ' -f4 | sort | uniq -c | awk '$1 == 1' | wc -l
