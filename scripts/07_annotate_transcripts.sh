
#!/usr/bin/env bash
set -euo pipefail

# Module I - Annotation of PsiCLASS transcripts
# Adapted from the Johns Hopkins University course

BASEDIR="/mnt/d/Bioinformatics/transcriptomics-data"
WORKDIR="${BASEDIR}/PsiCLASS"
OUTDIR="${WORKDIR}/WithGeneNames"

PROJECT_DIR="/mnt/c/Users/ianca/Documentos/Bioinformatics/transcriptomics-bioinformatics"

# Reference annotation
ANNOT="${BASEDIR}/reference/annotation/gencode.vM30.annotation.gtf"

# Annotation tools
ADDGENENAME="$(command -v add-genename || true)"
ADDSTRAND="${PROJECT_DIR}/course_materials/Module_I/docs/add-strand"

# Verify required tools and files
if [[ -z "$ADDGENENAME" ]]; then
    echo "ERROR: add-genename not found. Activate the psiclass environment." >&2
    exit 1
fi

for file in \
    "$ANNOT" \
    "${WORKDIR}/psiclass_vote.gtf" \
    "$ADDSTRAND"
do
    if [[ ! -s "$file" ]]; then
        echo "ERROR: Missing or empty file: $file" >&2
        exit 1
    fi
done

# Use temporary files to avoid overwriting valid outputs
# if an annotation step fails.
TMPDIR_ANNOT=$(mktemp -d "${WORKDIR}/.annotation_tmp.XXXXXX")
trap 'rm -rf "$TMPDIR_ANNOT"' EXIT

mkdir -p "${TMPDIR_ANNOT}/output" "$OUTDIR"

# Prepare GTF input list
printf '%s\n' \
    "${WORKDIR}/psiclass_vote.gtf" \
    > "${TMPDIR_ANNOT}/psiclass_gtf.list"

# Step 1 - Add gene names
echo "Adding gene names..."

"$ADDGENENAME" \
    "$ANNOT" \
    "${TMPDIR_ANNOT}/psiclass_gtf.list" \
    -o "${TMPDIR_ANNOT}/output"

ANNOTATED="${TMPDIR_ANNOT}/output/psiclass_vote.gtf"

if [[ ! -s "$ANNOTATED" ]]; then
    echo "ERROR: Gene annotation output is missing." >&2
    exit 1
fi

# Step 2 - Add strand information
echo "Processing strand information..."

FINAL="${TMPDIR_ANNOT}/psiclass_vote.withStrand.gtf"

"$ADDSTRAND" \
    "$ANNOT" \
    -r \
    < "$ANNOTATED" \
    > "$FINAL"

if [[ ! -s "$FINAL" ]]; then
    echo "ERROR: Strand annotation output is empty." >&2
    exit 1
fi

# Step 3 - Validate chromosome names and GTF structure
awk -F '\t' '
!/^#/ {
    if (NF != 9 || $1 ~ /^"",/ || $1 == "") {
        print "ERROR: Invalid GTF record at line", NR > "/dev/stderr"
        exit 1
    }
}' "$FINAL"

# Step 4 - Replace previous annotation outputs
mv -f "$ANNOTATED" "${OUTDIR}/psiclass_vote.gtf"
mv -f "$FINAL" "${OUTDIR}/psiclass_vote.withStrand.gtf"

echo "PsiCLASS transcript annotation completed!"
echo "Final GTF: ${OUTDIR}/psiclass_vote.withStrand.gtf"
