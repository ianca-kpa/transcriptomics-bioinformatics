
#!/usr/bin/env bash
set -euo pipefail

# ------------------------------------------------------------
# Module I - Download GENCODE reference annotation
# Johns Hopkins University - Bioinformatics Methods
# ------------------------------------------------------------

# 1. Project directories
BASEDIR="/mnt/d/Bioinformatics/transcriptomics-data/Module_I"
ANNOTDIR="${BASEDIR}/reference/annotation"

# 2. Reference annotation
RELEASE="M30"
ANNOTATION="gencode.vM30.annotation.gtf"

URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_${RELEASE}/${ANNOTATION}.gz"

# 3. Prepare directory
mkdir -p "$ANNOTDIR"
cd "$ANNOTDIR"

echo "Reference annotation: GENCODE mouse ${RELEASE}"
echo "Output directory: $ANNOTDIR"

# ------------------------------------------------------------
# 4. Download annotation
# ------------------------------------------------------------

if [[ -s "$ANNOTATION" ]]; then
    echo "Annotation already exists. Skipping download."

else
    echo "Downloading reference annotation..."

    wget -c "$URL" -O "${ANNOTATION}.gz"

    # Verify compressed file integrity
    echo "Checking gzip integrity..."
    gzip -t "${ANNOTATION}.gz"

    # Decompress annotation
    echo "Decompressing GTF..."
    gunzip "${ANNOTATION}.gz"
fi

# ------------------------------------------------------------
# 5. Verify annotation
# ------------------------------------------------------------

if [[ ! -s "$ANNOTATION" ]]; then
    echo "ERROR: Reference annotation is missing or empty." >&2
    exit 1
fi

echo "Reference annotation ready!"
echo "File: ${ANNOTDIR}/${ANNOTATION}"

ls -lh "$ANNOTATION"

