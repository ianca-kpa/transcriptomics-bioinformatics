
#!/usr/bin/env bash
set -e

BASEDIR=/mnt/d/Bioinformatics/transcriptomics-data/Module_I
CMDDIR=/mnt/c/Users/ianca/Documentos/Bioinformatics/transcriptomics-bioinformatics/Module_I/docs

P_ANNOT=${BASEDIR}/results/PsiCLASS/WithGeneNames/psiclass_vote.withStrand.gtf

ADDGN=${CMDDIR}/add_gn_name

WORKDIR=${BASEDIR}/results/DESeq2/Results

mkdir -p ${WORKDIR}/WithGeneNames

for i in results_1mo_4mo.psiclass.csv results_1mo_4mo.pval05.psiclass.csv results_1mo_4mo.qval05.psiclass.csv raw_counts.1mo_4mo.psiclass.csv normalized_counts.1mo_4mo.psiclass.csv
do
    ${ADDGN} ${P_ANNOT} < ${WORKDIR}/${i} > ${WORKDIR}/WithGeneNames/gn_${i}
done
