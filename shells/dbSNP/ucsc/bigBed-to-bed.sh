#!/bin/sh

TOP="/Users/bretbarnes/Documents"
EXE="${TOP}/tools/ucsc/bigBedToBed"

DAT_DIR="${TOP}/data/workhorse/annotation/Homo_sapiens/hg19/snp"
OUT_DIR="${TOP}/data/workhorse/annotation/Homo_sapiens/hg19/snp/Chromosomes"

NAM="dbSnp153"
BIG="${DAT_DIR}/${NAM}.bb"

CHROM="22"
CHROM="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT"
CHROM="2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 X Y MT"

echo "Starting..."
echo "  TOP='${TOP}'"
echo "  EXE='${EXE}'"
echo
echo "  DAT_DIR='${DAT_DIR}'"
echo "  OUT_DIR='${OUT_DIR}'"
echo
echo "  NAM='${NAM}'"
echo "  BIG='${BIG}'"
echo "  CHROM='${CHROM}'"
echo

mkdir -p ${OUT_DIR}

for CHR in $CHROM
do

    BED="${OUT_DIR}/chr${CHR}.${NAM}.bed"

    echo "  CHR='${CHR}'"
    echo "  BED='${BED}'"
    echo

    CMD="${EXE} -chrom=chr${CHR} ${BIG} ${BED}"
    echo "  CMD='${CMD}'"
    echo "..."
    echo
    eval $CMD

    CMD="bgzip -f ${BED}"
    echo "  CMD='${CMD}'"
    echo "..."
    echo
    eval $CMD
    
    CMD="tabix -p bed -f ${BED}.gz"
    echo "  CMD='${CMD}'"
    echo "..."
    echo
    eval $CMD
    
done

echo "done."

# End of file
