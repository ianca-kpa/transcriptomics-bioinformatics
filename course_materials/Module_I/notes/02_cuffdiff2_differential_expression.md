# Module I — Differential expression with Cuffdiff2

*Course: Bioinformatics Methods for Transcriptomics (Johns Hopkins University, Coursera)*  
*Tool: Cuffdiff2 v2.2.1 (Debian package)*

After aligning the six RNA-seq samples with STAR and assembling transcripts with PsiCLASS, I used Cuffdiff2 to compare the **1-month** and **4-month** mouse hippocampus samples. This note records the current inputs, the analysis steps, how I inspect the output, and a few command-line lessons. The later DESeq2 analyses and their gene-level comparisons are documented in Notes 03 and 04.

The checked-in workflow reads the six genome-wide STAR BAMs from `$BASEDIR/STAR/BAM/`, uses the genome-wide PsiCLASS GTF at `$BASEDIR/PsiCLASS/WithGeneNames/psiclass_vote.withStrand.gtf`, and writes to `$CUFFDIFF_DIR` (currently `$BASEDIR/Cuffdiff/`). I verified the current `gene_exp.diff` and `isoform_exp.diff` directly; this note reports only those current results.

## Step-by-step workflow

1. Check that the six coordinate-sorted STAR BAMs and genome-wide PsiCLASS GTF exist; confirm that the three 1mo samples are the control and the three 4mo samples are the test.
2. Verify `configs/paths.sh` points to the intended input and output locations.
3. If a new run is intended, run `./scripts/09_run_cuffdiff.sh` from the repository root. The script reuses `$CUFFDIFF_DIR` if it already exists and can replace results; do not run it over results you need to preserve.
4. Inspect `$CUFFDIFF_DIR/cuffdiff.log`, `run.info`, and the `gene_exp.diff` header before interpreting fields.
5. Count significant genes using both `status == "OK"` and `significant == "yes"`; the verified current total is 261 of 14,214 unique gene IDs.
6. Use p-value and log2FC cutoffs only as separately labeled exploratory summaries. Compare methods by gene ID, not aggregate counts.
7. See Note 04 for the verified gene-level comparison with sex-adjusted DESeq2 and its design limitations.

## Input data

| Group | Runs |
|---|---|
| Month 1 (first condition) | `SRR3734796`, `SRR3734797`, `SRR3734798` |
| Month 4 (second condition) | `SRR3734820`, `SRR3734824`, `SRR3734825` |

The checked-in script uses the coordinate-sorted STAR BAMs in `$BASEDIR/STAR/BAM/` and the genome-wide annotated PsiCLASS file:

```text
$BASEDIR/PsiCLASS/WithGeneNames/psiclass_vote.withStrand.gtf
```

The output directory is configured as `$CUFFDIFF_DIR` in `configs/paths.sh` (currently `$BASEDIR/Cuffdiff/`). The available `run.info` records Cuffdiff2 version 2.2.1.

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

The order matters: **Month 1 is the first condition and Month 4 is the second**, so a **positive log2 fold change** means higher estimated expression at 4 months. `scripts/09_run_cuffdiff.sh` specifies these groups in that order.

### A Bash error I made

My first Cuffdiff2 runs returned:

```text
Error: cuffdiff requires at least 2 SAM files
```

I had accidentally left spaces after the `\` used to continue the command across lines. Bash was not passing all the arguments to Cuffdiff2. Removing those spaces fixed the command. A useful rule for me: **the backslash has to be the last character on the line**.

I monitored the analysis with:

```bash
source configs/paths.sh
tail -f "${CUFFDIFF_DIR}/cuffdiff.log"
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

The current file has 14,214 gene records: 14,130 with `status == "OK"`, 82 with `status == "NOTEST"`, and 2 with `status == "FAIL"`. Not all records receive a usable test, so I require `status == "OK"` when applying statistical filters.

## Counting significant results

The `significant` flag is in column 14, and log2FC is in column 10. These examples assume I run them from `$CUFFDIFF_DIR`:

```bash
# Significant gene records
awk -F '\t' 'NR > 1 && $7 == "OK" && $14 == "yes"' gene_exp.diff | wc -l

# Significant isoform records
awk -F '\t' 'NR > 1 && $7 == "OK" && $14 == "yes"' isoform_exp.diff | wc -l
```

I use `status == "OK"` and `significant == "yes"` for both the gene- and isoform-level totals. The current results are:

| Result | Current count |
|---|---:|
| Gene records / unique gene IDs | 14,214 |
| `status == "OK"` | 14,130 |
| `status == "NOTEST"` | 82 |
| `status == "FAIL"` | 2 |
| Significant genes | 261 |
| Significant genes with positive `log2(fold_change)` | 60 |
| Significant genes with negative `log2(fold_change)` | 201 |
| Isoform records | 20,989 |
| Significant isoforms | 207 |

## Trying other log2FC cutoffs

The instructor also showed how to combine a usable-test status, an unadjusted p-value threshold, and a minimum absolute log2FC. These filters are exploratory and are not the same as the q-value-based significance flag:

```bash
# p < 0.05 and |log2FC| > 0.5
awk -F '\t' 'NR > 1 && $7 == "OK" && $12 < 0.05 && ($10 < -0.5 || $10 > 0.5)' gene_exp.diff | wc -l

# p < 0.05 and |log2FC| > 1
awk -F '\t' 'NR > 1 && $7 == "OK" && $12 < 0.05 && ($10 < -1 || $10 > 1)' gene_exp.diff | wc -l
```

On the current result file, these filters return 501 and 119 records, respectively. In my troubleshooting notes, I keep the incorrect predicate `($10 <= 1 || $10 > 1)`: every numeric fold change satisfies one side of this OR, so the expression does not impose an effect-size threshold. With `status == "OK"` and `p_value < 0.05`, this incorrect predicate returns 1,365 records in the current file. The corrected conditions above require values below the negative cutoff or above the positive cutoff.

These **p-value filters are exploratory**. They are not the same as selecting genes significant after multiple-testing correction with the Cuffdiff2 `q_value`/`significant` flag.

## Looking at individual genes

I sorted current significant results by `log2(fold_change)`:

```bash
# Largest positive changes
awk -F '\t' 'NR > 1 && $14 == "yes" && $7 == "OK"' gene_exp.diff | sort -k10,10nr | head -n 10

# Largest negative changes
awk -F '\t' 'NR > 1 && $14 == "yes" && $7 == "OK"' gene_exp.diff | sort -k10,10n | head -n 10
```

The current result file has these gene-symbol matches. I retain `gene_id` because one symbol can correspond to multiple Cuffdiff2 records:

| Gene | `gene_id` | `log2(fold_change)` | `q_value` | `significant` |
|---|---|---:|---:|---|
| `Il33` | `chr19.103137` | 1.84387 | 0.00494056 | yes |
| `Cd74` | `chr18.100510` | 1.60748 | 0.00494056 | yes |
| `Serpinb1a` | `chr13.77419` | 1.23373 | 0.00494056 | yes |
| `Eln` | `chr5.34342` | −1.61125 | 0.00494056 | yes |
| `Eln` | `chr5.34345` | −2.40070 | 0.00883125 | yes |
| `Col3a1` | `chr1.1658` | −1.27127 | 0.305538 | no |
| `Col3a1` | `chr1.1659` | −1.73993 | 0.750991 | no |
| `Col3a1` | `chr1.1662` | −1.09552 | 0.429232 | no |
| `Col3a1` | `chr1.1663` | −1.00874 | 0.822282 | no |
| `Col3a1` | `chr1.1664` | −1.98925 | 0.00494056 | yes |
| `Enpp6` | `chr8.50741` | −1.88641 | 0.00494056 | yes |
| `Gsn` | `chr2.8576` | −1.71598 | 0.00494056 | yes |
| `Xist` | `chrX.106929` | −1.70467 | 0.00494056 | yes |

Both current `Eln` records are significant. Of the five `Col3a1` records, only `chr1.1664` is significant. These gene-level observations are descriptive; the unequal sex composition limits an age-only interpretation.

## Keeping a short summary

I use `scripts/10_summarize_cuffdiff.sh` to summarize the current result file:

```bash
source configs/paths.sh
./scripts/10_summarize_cuffdiff.sh | tee "${CUFFDIFF_DIR}/cuffdiff_summary.txt"
```

The current verified totals are **261 significant genes** (60 positive and 201 negative) and **207 significant isoforms**. The script also prints p-value and direction summaries; I cross-check them against the status-aware filters and source tables above before reporting them. `tee` overwrites `cuffdiff_summary.txt`, so I check the destination before saving.

## Interpretation and next steps

- A positive log2 fold change means higher estimated expression at 4mo than at 1mo; it is not evidence of causality.
- Cuffdiff2 and DESeq2 use different statistical methods and model specifications. I compare their gene IDs and effect directions in Note 04 rather than treating their p-values as interchangeable.
- The sex distribution is unbalanced (three females and no males at 1mo; one female and two males at 4mo), which limits an age-only interpretation.

## Command reference — how the Bash analysis works

I want to be able to reproduce these filters later, so I am keeping the commands and their explanations together. The commands that read `gene_exp.diff` assume the current directory is the Cuffdiff output directory; in the current scripts, that is `$CUFFDIFF_DIR`.

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
awk -F '\t' 'NR > 1 && $7 == "OK" && $14 == "yes" && $10 > 0' gene_exp.diff | wc -l
```

- `awk` processes one record at a time.
- `-F '\t'` explicitly defines the field separator as a tab; this is more reliable for a TSV than the default whitespace splitting.
- `$14 == "yes"` checks the significance flag. `==` is a **comparison**; it is not variable assignment.
- `&&` means **AND**: both conditions must hold.
- `$10 > 0` selects positive log2FC, so higher estimated expression in **Month 4 compared with Month 1**.
- `wc -l` counts how many matching rows the preceding command printed.

For negative log2 fold changes I change only `> 0` to `< 0`:

```bash
awk -F '\t' 'NR > 1 && $7 == "OK" && $14 == "yes" && $10 < 0' gene_exp.diff | wc -l
```

A shortcut is `grep -c $'\tyes$' gene_exp.diff`, which counts lines ending with a tab followed by `yes`. It relies on the file format and does not check the status column; selecting fields with `awk` is clearer for the reported counts.

### Combining p-value with a log2FC cutoff

```bash
awk -F '\t' '$7 == "OK" && ($10 < -0.5 || $10 > 0.5) && $12 < 0.05' gene_exp.diff | wc -l
```

- `||` means **OR**: I accept a large negative **or** a large positive change.
- Parentheses group the two allowed log2FC directions.
- `&& $12 < 0.05` additionally requires the **unadjusted** p-value to be below 0.05.
- For a stricter effect-size threshold, replace `-0.5` and `0.5` with `-1` and `1`.

I once typed `($10 <= 1 || $10 > 1)`. For numeric values, that condition covers everything (a value is either at most 1 or greater than 1), so I was effectively filtering only on p-value. The correct second condition is:

```bash
awk -F '\t' '$7 == "OK" && ($10 < -1 || $10 > 1) && $12 < 0.05' gene_exp.diff | wc -l
```

These filters are *not* FDR-controlled because they use `p_value`, not `q_value`.

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
awk -F '\t' 'NR == 1 || ($7 == "OK" && $14 == "yes" && $10 > 0)' gene_exp.diff > genes_upregulated.tsv
```

- `NR` is the current input record number; `NR == 1` preserves the header.
- `>` **creates or overwrites** the output file. It will replace an existing `genes_upregulated.tsv`.
- `>>` would append instead; I generally do **not** want it for a full filtered table because I might create duplicate results by rerunning the command.

```bash
source configs/paths.sh
./scripts/10_summarize_cuffdiff.sh | tee "${CUFFDIFF_DIR}/cuffdiff_summary.txt"
```

- `./` runs the script relative to the current directory; run it from the repository root, where `scripts/10_summarize_cuffdiff.sh` is located.
- `|` sends its displayed output to another command.
- `tee` displays the output **and** saves it to a text file. By default, it overwrites the file; `tee -a` appends.
- If I run a command from a different folder, I must adjust its relative file paths.

### Checks I want to repeat before interpreting results

```bash
source configs/paths.sh
head -n 1 "${CUFFDIFF_DIR}/gene_exp.diff" | tr '\t' '\n' | nl -ba
cut -f7 "${CUFFDIFF_DIR}/gene_exp.diff" | sort | uniq -c
ls -lh "${CUFFDIFF_DIR}/gene_exp.diff"
```

These checks use the configured output path. They verify column order, test statuses and file presence; they do not prove that the experimental design is free of confounding factors. Inspect the sex metadata before interpreting the Month 1/Month 4 difference biologically.

## Related notes

- [Note 03 — DESeq2 differential expression](03_deseq2_differential_expression.md) describes the genome-wide unadjusted DESeq2 workflow and its input count matrix.
- [Note 04 — Sex-adjusted DESeq2 and method comparison](04_sex_adjusted_deseq2_and_method_comparison.md) compares verified gene-level results and explains the limits of that comparison.

*Based on the Johns Hopkins University Coursera exercise, with local adaptations and results from my own run.*
