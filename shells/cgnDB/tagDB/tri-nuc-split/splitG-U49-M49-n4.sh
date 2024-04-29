#!/bin/sh

EXE_NAME=Infinium_Methylation_Workhorse/scripts/R/swifthoof/swifthoof_main.R
EXE_A=/repo/${EXE_NAME}
EXE_B=/Users/bretbarnes/Documents/tools/${EXE_NAME}

RSCRIPT_A=/usr/local/bin/Rscript
RSCRIPT_B="/Users/bretbarnes/Documents/tools/${EXE_NAME}"
RSCRIPT_C="/illumina/scratch/darkmatter/tools/Infinium_Methylation_Workhorse/scripts/R/swifthoof/swifthoof_main.R"

if [ -e ${RSCRIPT_A} ]; then
    RSCRIPT=${RSCRIPT_A}
    TOP_DIR="/repo"
elif [ -e ${RSCRIPT_B} ]; then
    RSCRIPT=${RSCRIPT_B}
    TOP_DIR="/Users/bretbarnes/Documents"

    DAT_DIR="${TOP_DIR}/data/improbe/scratch/cgnDB/dbSNP_Core4/design-output/prbs-p49"
    OUT_DIR="${TOP_DIR}/data/improbe/scratch/cgnDB/dbSNP_Core4/design-output/prbs-P49-split3.64"

elif [ -e ${RSCRIPT_C} ]; then
    RSCRIPT=${RSCRIPT_C}
    TOP_DIR="/illumina/scratch/darkmatter"
    DAT_DIR="${TOP_DIR}/scratch/improbe/cgnDB/dbSNP_Core4/prbs-p49"
    OUT_DIR="${TOP_DIR}/scratch/improbe/cgnDB/dbSNP_Core4/prbs-P49-split3.64"

else
    echo "Unrecognized Rscript EXE!"
    exit
fi

U49_TSV="${DAT_DIR}/probe_U49_cgn-table.csv.gz"
M49_TSV="${DAT_DIR}/probe_M49_cgn-table.csv.gz"

NUCS="GAA GAC GAG GAT GCA GCC GCG GCT GGA GGC GGG GGT GTA GTC GTG GTT"

for NUC in $NUCS
do
    echo "NUC=${NUC}"
done

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
