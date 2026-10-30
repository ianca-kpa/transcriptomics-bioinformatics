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

## Command reference — how the Bash analysis works

I want to be able to reproduce these filters later, so I am keeping the commands and their explanations together. Unless indicated otherwise, I run them from `results/Cuffdiff/`.

### The Cuffdiff command, argument by argument

```bash
cuffdiff \
    -o "$WORKDIR" \
    -p 2 \
    "$ANNOT" \
    "$CTRLBAMS" \
    "$TESTBAMS" \
    &> "${WORKDIR}/cuffdiff.log"
```

| Element | Meaning |
|---|---|
| `cuffdiff` | Run the Cuffdiff2 program. |
| `-o "$WORKDIR"` | Output directory. Cuffdiff2 writes its result tables there. |
| `-p 2` | Ask for two threads. It is not a memory limit. |
| `"$ANNOT"` | The PsiCLASS GTF file defining transcript models. |
| `"$CTRLBAMS"` | One comma-separated argument holding the three **Month 1** BAM paths. |
| `"$TESTBAMS"` | One comma-separated argument holding the three **Month 4** BAM paths. |
| `&>` | Redirect standard output **and** error messages to `cuffdiff.log` (Bash syntax). |

The commas separate **replicates within a condition**; the spaces between `"$CTRLBAMS"` and `"$TESTBAMS"` separate the two conditions. I forgot that a backslash followed by a space does **not** continue a line: the backslash must be the last character before the newline.

### Inspecting columns before filtering

```bash
head -n 1 gene_exp.diff | tr '\t' '\n' | nl -ba
```

- `head -n 1` reads only the header line.
- `tr '\t' '\n'` replaces tab separators with line breaks, putting each column name on its own line.
- `nl -ba` numbers all lines, starting at 1. This lets me match `$10` to `log2(fold_change)`, `$12` to `p_value` and `$14` to `significant` **in this particular file**.
- In `awk`, `$10` means the **tenth field**, not ten dollars; `$0` means the **entire record/line**.

### Counting test statuses

```bash
cut -f7 gene_exp.diff | sort | uniq -c
```

1. `cut -f7` extracts only the seventh **tab-separated** field (`status`).
2. `sort` groups identical values together.
3. `uniq -c` counts consecutive repetitions; this is why the `sort` step matters.
4. The header `status` appears once and is not a biological result.

### `awk` filters: significance and direction

```bash
awk -F '\t' '$14 == "yes" && $10 > 0' gene_exp.diff | wc -l
```

- `awk` processes one record at a time.
- `-F '\t'` explicitly defines the field separator as a tab; this is more reliable for a TSV than the default whitespace splitting.
- `$14 == "yes"` checks the significance flag. `==` is a **comparison**; it is not variable assignment.
- `&&` means **AND**: both conditions must hold.
- `$10 > 0` selects positive log2FC, so higher estimated expression in **Month 4 compared with Month 1**.
- `wc -l` counts how many matching rows the preceding command printed.

For downregulated results I change only `> 0` to `< 0`:

```bash
awk -F '\t' '$14 == "yes" && $10 < 0' gene_exp.diff | wc -l
```

A shortcut is `grep -c $'\tyes$' gene_exp.diff`, which counts lines ending with a tab followed by `yes`. It relies on the file format; selecting column 14 with `awk` is clearer when I need to combine conditions.

### Combining p-value with a log2FC cutoff

```bash
awk -F '\t' '($10 < -0.5 || $10 > 0.5) && $12 < 0.05' gene_exp.diff | wc -l
```

- `||` means **OR**: I accept a large negative **or** a large positive change.
- Parentheses group the two allowed log2FC directions.
- `&& $12 < 0.05` additionally requires the **unadjusted** p-value to be below 0.05.
- For a stricter effect-size threshold, replace `-0.5` and `0.5` with `-1` and `1`.

I once typed `($10 <= 1 || $10 > 1)`. For numeric values, that condition covers everything (a value is either at most 1 or greater than 1), so I was effectively filtering only on p-value. The correct second condition is:

```bash
awk -F '\t' '($10 < -1 || $10 > 1) && $12 < 0.05' gene_exp.diff | wc -l
```

The result was **131**, while the ±0.5 filter returned **507**. These counts are *not* FDR-controlled findings because they use `p_value`, not `q_value`.

### Sorting and previewing candidate genes

```bash
awk -F '\t' '$14 == "yes" && $7 == "OK"' gene_exp.diff \
    | sort -k10,10nr \
    | head -n 10
```

- `$7 == "OK"` restricts to records for which Cuffdiff2 reported a usable test.
- `sort -k10,10nr` sorts numerically (`n`), in reverse (`r`, largest first), using field 10 as the key. With these tab-delimited output lines, this works because no field contains embedded spaces; to be explicit about tabs, `sort -t $'\t' -k10,10nr` is safer.
- `head -n 10` shows only the first ten records.
- Replace `-k10,10nr` by `-k10,10n` to see the most negative log2FC values first.

Sorting by log2FC is exploratory. I should also inspect expression levels and adjusted p-values before choosing genes to investigate.

### Saving filtered results and summaries

```bash
awk -F '\t' 'NR == 1 || ($14 == "yes" && $10 > 0)' gene_exp.diff > genes_upregulated.tsv
```

- `NR` is the current input record number; `NR == 1` preserves the header.
- `>` **creates or overwrites** the output file. It will replace a previous `genes_upregulated.tsv`.
- `>>` would append instead; I generally do **not** want it for a full filtered table because I might create duplicate results by rerunning the command.

```bash
./scripts/10_summarize_cuffdiff.sh | tee results/Cuffdiff/cuffdiff_summary.txt
```

- `./` runs the script relative to my current directory (this example assumes I am in `Module_I`, not inside `results/Cuffdiff/`).
- `|` sends its displayed output to another command.
- `tee` displays the output **and** saves it to a text file. By default, it overwrites the file; `tee -a` appends.
- If I run a command from a different folder, I must adjust its relative file paths.

### Checks I want to repeat before interpreting results

```bash
head -n 1 results/Cuffdiff/gene_exp.diff | tr '\t' '\n' | nl -ba
cut -f7 results/Cuffdiff/gene_exp.diff | sort | uniq -c
ls -lh results/Cuffdiff/gene_exp.diff
```

These commands are meant to run from `Module_I`. They check column order, test statuses and whether the result file exists. They do not prove that the experimental design is free of confounding factors: I still need to inspect sex metadata before interpreting the Month 1/Month 4 difference biologically.

## Next Steps

After the Cuffdiff2 analysis, I continued with DESeq2 to analyze differential gene expression using read counts.

The next note covers BAM sorting, gene counting, DESeq2 normalization, gene name annotation, filtering, and the use of biological sex as a covariate.

See [Note 3 — DESeq2 Differential Expression](03_deseq2_differential_expression.md).

*Based on the Johns Hopkins University Coursera exercise, with local adaptations and results from my own run.*
