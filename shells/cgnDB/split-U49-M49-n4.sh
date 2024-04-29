#!/bin/sh

TOP_DIR="/Users/bretbarnes/Documents"
DAT_DIR="${TOP_DIR}/data/improbe/scratch/cgnDB/dbSNP_Core4/design-output/prbs-p49"
OUT_DIR="${TOP_DIR}/data/improbe/scratch/cgnDB/dbSNP_Core4/design-output/prbs-P49-split"

U49_TSV="${DAT_DIR}/probe_U49_cgn-table.tsv.gz"
M49_TSV="${DAT_DIR}/probe_M49_cgn-table.tsv.gz"

NUCS="AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT"

echo "DAT_DIR = ${DAT_DIR}"
echo "OUT_DIR = ${OUT_DIR}"
echo "DAT_DIR = ${DAT_DIR}"
echo ""
echo "U49_TSV = ${U49_TSV}"
echo "M49_TSV = ${M49_TSV}"
echo ""

CMD="mkdir -p ${OUT_DIR}"
eval ${CMD}

for NUC in $NUCS
do

    OUT_TSV="${OUT_DIR}/${NUC}-probe_U49_cgn-table.tsv.gz"
    CMD="gzip -dc ${U49_TSV} | grep '^${NUC}' | gzip -c -> ${OUT_TSV}"
    
    echo "CMD = '$CMD'"
    eval ${CMD}
    echo "done"
    echo ""

    OUT_TSV="${OUT_DIR}/${NUC}-probe_M49_cgn-table.tsv.gz"
    CMD="gzip -dc ${M49_TSV} | grep '^${NUC}' | gzip -c -> ${OUT_TSV}"

    echo "CMD = '$CMD'"
    eval ${CMD}
    echo "done"
    echo ""

done

echo "done."

# End of file
