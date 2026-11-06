# Environment and execution

The workflow was developed on Ubuntu/WSL2 using Bash, Python, R,
Jupyter, and Conda.

Several scripts use absolute local paths. Review these paths before
running the workflow on another machine.

## Dependencies

| Component | Purpose |
|---|---|
| SRA Toolkit | SRA download and FASTQ conversion |
| STAR | Reference indexing and alignment |
| SAMtools | BAM inspection and sorting |
| PsiCLASS | Transcript assembly |
| Cuffdiff2 | Differential expression |
| numpy, pandas, matplotlib, IPython | Notebook processing and display |
| DESeq2 | Differential expression |
| GenomicFeatures, GenomicAlignments, Rsamtools | Annotation and counting |
| SummarizedExperiment, txdbmaker | Count objects and transcript databases |
| EnhancedVolcano, ggplot2 | Figures in the whole-genome scripts |

The notebook expects PsiCLASS in a Conda environment named `psiclass`
and Rscript in an environment named `transcriptomics`.

Some tools may be installed outside Conda. Course utilities such as
`add_gn_name` are referenced in the local `Module_I/docs/` directory
and must be obtained from the course resources.

## Record the analysis environment

Run these commands on the machine used for the analysis, from the
repository root:

```bash
mkdir -p Module_I/environment

conda env export -n transcriptomics --no-builds > Module_I/environment/transcriptomics.yml
conda env export -n psiclass --no-builds > Module_I/environment/psiclass.yml

STAR --version > Module_I/environment/STAR-version.txt
samtools --version > Module_I/environment/SAMtools-version.txt