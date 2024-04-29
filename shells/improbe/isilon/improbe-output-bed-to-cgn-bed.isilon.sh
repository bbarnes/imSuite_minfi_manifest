#!/bin/bash

if [ "$#" -lt 1 ]; then
    echo "Usage: $0 genomeBuild"
    exit 1
fi
BUILD=$1

# TOP_DIR="/Users/bretbarnes/Documents"
# DAT_DIR="${TOP_DIR}/data/improbe/GRCh36-38.mm10.FWD_SEQ.21022021/designInput"
# OUT_DIR="${TOP_DIR}/data/improbe/GRCh36-38.mm10.FWD_SEQ.21022021/designOutput/bed"

TOP_DIR="/illumina/scratch/darkmatter/Projects/dbCGN"
DAT_DIR="${TOP_DIR}/data/cgnDB/GRCh36-38.mm10.FWD_SEQ.21022021/designInput"
OUT_DIR="${TOP_DIR}/data/cgnDB/GRCh36-38.mm10.FWD_SEQ.21022021/designOutput/bed"

DES="${DAT_DIR}/${BUILD}.cgn.sorted.bed.gz"
OUT_BED="${OUT_DIR}/${BUILD}-21022021.cgn.sorted.bed.gz"

echo "BUILD=${BUILD}"
echo "DAT_DIR=${DAT_DIR}"
echo "OUT_DIR=${OUT_DIR}"
echo ""
echo "DES=${DES}"
echo "OUT_BED=${OUT_BED}"
echo ""

mkdir -p ${OUT_DIR}

gzip -dc ${DES} | cut -f 1-4,6 | \
    # head -n 20 | \
    perl -pe 's/\n$//; @d=split(/\t/,$_); s/^.*$//; $d[3] =~ s/_.*$//; print "$d[0]\t$d[1]\t$d[2]\t$d[3]\t.\t$d[4]\n"' | \
    gzip -c -> ${OUT_BED}

echo "done."

# End of file
