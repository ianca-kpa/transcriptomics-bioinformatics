
#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

source "${PROJECT_DIR}/configs/paths.sh"

CMDDIR="${PROJECT_DIR}/course_materials/Module_I/docs"

P_ANNOT="${PSICLASS_DIR}/WithGeneNames/psiclass_vote.withStrand.gtf"

ADDGN=${CMDDIR}/add_gn_name

WORKDIR="$DESEQ2_RESULTS_DIR"

mkdir -p ${WORKDIR}/WithGeneNames

for i in results_1mo_4mo.psiclass.csv results_1mo_4mo.pval05.psiclass.csv results_1mo_4mo.qval05.psiclass.csv raw_counts.1mo_4mo.psiclass.csv normalized_counts.1mo_4mo.psiclass.csv
do
    ${ADDGN} ${P_ANNOT} < ${WORKDIR}/${i} > ${WORKDIR}/WithGeneNames/gn_${i}
done
