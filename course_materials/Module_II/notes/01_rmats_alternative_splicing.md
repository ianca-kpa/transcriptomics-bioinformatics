# Module II — Alternative splicing with rMATS

*Dataset: GSE83931 — mouse hippocampus RNA-seq*

*Comparison: 1 month (1mo) vs 4 months (4mo)*

This note documents the initial skipped-exon analysis in Module II, based on the checked-in rMATS scripts and the result files currently available under the configured data directory.

## Study design

The rMATS script compares six samples:

| Group | SRA runs | Sex |
|---|---|---|
| 1mo (group 1) | SRR3734796, SRR3734797, SRR3734798 | Female, female, female |
| 4mo (group 2) | SRR3734820, SRR3734824, SRR3734825 | Female, male, male |

The groups are unbalanced by sex. The rMATS comparison is between age groups and does not adjust for sex, so interpret any age-related differences cautiously.

## Inputs and run settings

[`scripts/18_run_rmats.sh`](../../../scripts/18_run_rmats.sh) reads coordinate-sorted STAR BAMs from `$STAR_BAM_DIR` and the reference GTF configured as `$REFERENCE_GTF`. It creates comma-separated BAM lists for the six samples above and runs:

- `-t single`
- `--readLength 100`
- `--novelSS`
- `--allow-clipping`
- `--nthread 2`

The output directory is `$RMATS_DIR/Ctrl-vs-Test/`; temporary files go to `$RMATS_DIR/tmp_Ctrl-vs-Test/`. These settings describe the script configuration; confirm read layout and read lengths for every sample before rerunning it.

[`scripts/17a_check_read_lengths.sh`](../../../scripts/17a_check_read_lengths.sh) can check read lengths, SAM paired-end flags, and BAM integrity. Its checked-in defaults inspect only the first 10,000 primary mapped alignments from `SRR3734796`. It does not validate all six samples unless its mode is changed to `six_samples`; that mode also inspects only the first 10,000 alignments per sample.

## Event filtering and direction

The primary event table is `SE.MATS.JC.txt`. [`scripts/19_summarize_rmats.sh`](../../../scripts/19_summarize_rmats.sh) filters skipped-exon events using:

```text
FDR < 0.05
|IncLevelDifference| >= 0.05
```

It preserves the original rMATS columns and sorts the exported rows by FDR (ascending), then absolute inclusion-level difference (descending). The output is:

```text
$RMATS_DIR/Ctrl-vs-Test/SE.MATS.JC.significant_FDR0.05_dPSI0.05.tsv
```

`IncLevelDifference` is mean inclusion in group 1 minus mean inclusion in group 2:

- Positive values indicate higher exon inclusion in 1mo.
- Negative values indicate higher exon inclusion in 4mo.

## Results checked

The current configured data directory contains the six BAM inputs, the reference GTF, the rMATS event table, and the filtered event table. I independently checked the event table against the thresholds above:

| Check | Events |
|---|---:|
| Skipped-exon events in `SE.MATS.JC.txt` | 20,588 |
| Events with invalid or missing FDR/delta values | 0 |
| Events passing both thresholds | 111 |
| Higher inclusion in 1mo (positive difference) | 56 |
| Higher inclusion in 4mo (negative difference) | 55 |
| Data rows in the exported significant-event table | 111 |

These are **events**, not counts of distinct genes. I compared the exported rows with the source records passing both thresholds; all 111 rows matched. This verifies the filter output, not the biological interpretation or the suitability of all run settings.

## Reproducibility and repository caveats

The Module II run list contains nine runs across the 1-, 2-, and 4-month groups. `scripts/05b_align_all_samples_star.sh` aligns all nine; `scripts/18_run_rmats.sh` deliberately uses only the six 1mo and 4mo samples for its current comparison. This is a subset analysis, not a run-list error.

The SRA conversion and STAR alignment scripts now use the checked-in Module II run list and the same configured `$BASEDIR/SRA/` and `$BASEDIR/FASTQ/` locations. Run the Module II preparation steps from the repository root with:

```bash
SRA_RUNS_FILE=course_materials/Module_II/metadata/selected_runs_m2.txt \
    ./scripts/02_download_sra.sh
./scripts/03_convert_sra_to_fastq.sh
./scripts/05b_align_all_samples_star.sh
```

The conversion script defaults to the Module II list; `SRA_RUNS_FILE` can be set to another list when needed.

The read-validation script's default mode still checks only the first 10,000 primary mapped alignments from one sample, not the complete rMATS cohort. Large BAMs, references, and rMATS output files are local and are not tracked in the repository.

## Re-running the analysis

After confirming the sample list, FASTQ/BAM paths, read layout, read lengths, reference GTF, and local `configs/paths.sh`, run from the repository root in an environment containing SAMtools and rMATS:

```bash
./scripts/17a_check_read_lengths.sh
./scripts/18_run_rmats.sh
./scripts/19_summarize_rmats.sh
```

The first command follows the script's current one-sample defaults. It is not a substitute for checking every sample. The rMATS script uses the configured six-run groups; the summary script reads its output and writes the filtered skipped-exon table.

## Next steps

- Reconcile the nine-run Module II list with the six samples used for this comparison.
- Validate the read layout and read-length assumptions for all six BAMs.
- Review event-level results with gene IDs and genomic coordinates; do not infer gene-level counts from event totals.
- Consider the sex imbalance before making biological claims.
