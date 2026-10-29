# Module I — Differential expression with Cuffdiff2

*Course: Bioinformatics Methods for Transcriptomics (Johns Hopkins University, Coursera)*  
*Tool: Cuffdiff2 v2.2.1 (Debian package)*

After aligning the six RNA-seq samples with STAR and assembling transcripts with PsiCLASS, I used Cuffdiff2 to compare the **1-month** and **4-month** mouse hippocampus samples. This note records the command I used, a couple of mistakes I had to fix, and the results I obtained. My next step is to try DESeq2 and compare the approaches.

## Input data

| Group | Runs |
|---|---|
| Month 1 (first condition) | `SRR3734796`, `SRR3734797`, `SRR3734798` |
| Month 4 (second condition) | `SRR3734820`, `SRR3734824`, `SRR3734825` |

I used coordinate-sorted STAR BAMs from `results/STAR/BAM/` and the annotated PsiCLASS file:

```text
results/PsiCLASS/WithGeneNames/psiclass_vote.withStrand.gtf
```

My output directory is `results/Cuffdiff/`, which points to the D: drive through the repository's `results` symbolic link.

## Running Cuffdiff2

I adapted the instructor's `com.cuffdiff` command in `scripts/09_run_cuffdiff.sh`. The central command is:

```bash
cuffdiff \
    -o "$WORKDIR" \
    -p 2 \
    "$ANNOT" \
    "$CTRLBAMS" \
    "$TESTBAMS" \
    &> "${WORKDIR}/cuffdiff.log"
```

- `-o`: where to save the results.
- `-p 2`: use two threads on my laptop.
- `ANNOT`: PsiCLASS GTF annotation.
- `CTRLBAMS` / `TESTBAMS`: comma-separated BAM paths, with three replicates in each group.
- `&>`: save program messages and errors in the log file.

The order matters: **Month 1 is the first condition and Month 4 is the second**, so a **positive log2 fold change** means higher estimated expression at 4 months.

### A Bash error I made

My first Cuffdiff2 runs returned:

```text
Error: cuffdiff requires at least 2 SAM files
```

I had accidentally left spaces after the `\` used to continue the command across lines. Bash was not passing all the arguments to Cuffdiff2. Removing those spaces fixed the command. A useful rule for me: **the backslash has to be the last character on the line**.

I monitored the analysis with:

```bash
tail -f results/Cuffdiff/cuffdiff.log
pgrep -af cuffdiff
```

The program also warned that it could not contact its update server. That warning did not explain the earlier command error.

## Understanding `gene_exp.diff`

This was the main output used in the course exercise. I checked the actual header rather than assuming which column contained each result:

```bash
head -n 1 gene_exp.diff | tr '\t' '\n' | nl -ba
```

The file produced in my run has **14 columns**:

| Column | Field | What I used it for |
|---:|---|---|
| 1–4 | `test_id`, `gene_id`, `gene`, `locus` | Gene identifiers, name, genomic location |
| 5–6 | `sample_1`, `sample_2` | Condition names (`q1`, `q2`) |
| 7 | `status` | Whether the statistical test was performed |
| 8–9 | `value_1`, `value_2` | Estimated abundance in each condition |
| 10 | `log2(fold_change)` | Direction and size of the change |
| 11 | `test_stat` | Test statistic |
| 12 | `p_value` | Unadjusted p-value |
| 13 | `q_value` | Multiple-testing-adjusted value |
| 14 | `significant` | Cuffdiff2 `yes` / `no` flag |

For this comparison, log2FC of **+1** corresponds to approximately twice the estimated expression at Month 4; **−1** corresponds to half the expression.

## Checking the tests

```bash
cut -f7 gene_exp.diff | sort | uniq -c
```

This returned (excluding the header):

| Status | Records |
|---|---:|
| `OK` | 14,233 |
| `NOTEST` | 59 |
| `HIDATA` | 1 |
| **Total** | **14,293** |

Not all records received a usable test. I use the `OK` status when reviewing individual results.

## Counting significant results

The `significant` flag is in column 14, and log2FC is in column 10. I used these commands from the `results/Cuffdiff/` directory:

```bash
# Significant gene records
awk -F '\t' '$14 == "yes"' gene_exp.diff | wc -l

# Higher expression at Month 4
awk -F '\t' '$14 == "yes" && $10 > 0' gene_exp.diff | wc -l

# Lower expression at Month 4
awk -F '\t' '$14 == "yes" && $10 < 0' gene_exp.diff | wc -l

# Significant isoform records
grep -c $'\tyes$' isoform_exp.diff
```

My results:

| Result | Count |
|---|---:|
| Significant gene records | **265** |
| Upregulated at Month 4 | **63** |
| Downregulated at Month 4 | **202** |
| Significant isoform records | **212** |

The directions add up: **63 + 202 = 265**. Gene-level and isoform-level counts are different quantities; 212 significant isoforms does not mean 212 of the 265 genes.

## Trying other log2FC cutoffs

The instructor also showed how to combine an unadjusted p-value threshold with a minimum absolute log2FC:

```bash
# p < 0.05 and |log2FC| > 0.5
awk -F '\t' '($10 < -0.5 || $10 > 0.5) && $12 < 0.05' gene_exp.diff | wc -l

# p < 0.05 and |log2FC| > 1
awk -F '\t' '($10 < -1 || $10 > 1) && $12 < 0.05' gene_exp.diff | wc -l
```

| Filter | My count |
|---|---:|
| p < 0.05 and \|log2FC\| > 0.5 | **507** |
| p < 0.05 and \|log2FC\| > 1 | **131** |

I initially made a mistake with the second filter: I wrote `$10 <= 1 || $10 > 1`, which covers essentially every numeric log2FC and returned **1,362** rows. The intended condition is `$10 < -1 || $10 > 1`, which selected **131**. That was a useful reminder to check logical conditions before interpreting a count.

These **p-value filters are exploratory**. They are not the same as selecting genes significant after multiple-testing correction with the Cuffdiff2 `q_value`/`significant` flag.

## Looking at individual genes

I sorted significant results by log2FC:

```bash
# Largest positive changes
awk -F '\t' '$14 == "yes" && $7 == "OK"' gene_exp.diff | sort -k10,10nr | head -n 10

# Largest negative changes
awk -F '\t' '$14 == "yes" && $7 == "OK"' gene_exp.diff | sort -k10,10n | head -n 10
```

Some genes from the output:

| Gene | log2FC | q-value |
|---|---:|---:|
| `Il33` | +1.84754 | 0.00508321 |
| `Cd74` | +1.66516 | 0.00508321 |
| `Serpinb1a` | +1.24378 | 0.00508321 |
| `Eln` | −2.41140 | 0.00918258 |
| `Col3a1` | −2.03349 | 0.00508321 |
| `Enpp6` | −1.91252 | 0.00508321 |
| `Gsn` | −1.71437 | 0.00508321 |
| `Xist` | −1.70287 | 0.00508321 |

I have not investigated these genes in detail yet. One result I want to check carefully is **Xist**, because its expression can be strongly related to sex. The selected groups may differ in sex composition, which complicates an age-only interpretation. This is a course comparison, not enough evidence on its own for a biological conclusion about development.

## Keeping a short summary

I created `scripts/10_summarize_cuffdiff.sh` so I would not need to type the counting commands every time. To print its output and save it:

```bash
./scripts/10_summarize_cuffdiff.sh | tee results/Cuffdiff/cuffdiff_summary.txt
```

The summary returned **265 significant genes**, **212 significant isoforms**, **63 upregulated genes**, **202 downregulated genes**, and **131 records** with p < 0.05 and |log2FC| > 1. Full results stay on D:, while the scripts and notes are versioned in Git.

## What I want to do next

- Check the sample metadata, particularly the sex of each animal.
- Compare the Cuffdiff2 results with **DESeq2** using an appropriate gene-level count matrix.
- Decide which output tables and QC summaries are useful to keep as small reproducible examples in the GitHub repository.

*Based on the Johns Hopkins University Coursera exercise, with local adaptations and results from my own run.*
