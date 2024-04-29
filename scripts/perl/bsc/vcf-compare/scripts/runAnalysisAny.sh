## Example of how to run simluation analysis pipeline

BIN=/home/bbarnes/vcf-compare/scripts

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 candidate.vcf reference.vcf outputDir"
    exit 1
fi

CAN_VCF=$1
REF_VCF=$2
OUT_DIR=$3

if [ ! -d $OUT_DIR ]; then
    mkdir $OUT_DIR
fi

zcat ${CAN_VCF} | ${BIN}/vcfToBedSplit.pl -o ${OUT_DIR}/candidate
zcat ${REF_VCF} | ${BIN}/vcfToBedSplit.pl -o ${OUT_DIR}/reference

${BIN}/windowBed -w 10 \
    -a ${OUT_DIR}/candidate.bed \
    -b ${OUT_DIR}/reference.bed \
    >  ${OUT_DIR}/intersect.txt

${BIN}/vcfCompare.pl \
    -c ${OUT_DIR}/candidate.vcf \
    -r ${OUT_DIR}/reference.vcf \
    -i ${OUT_DIR}/intersect.txt \
    -o ${OUT_DIR}/results \
    >  ${OUT_DIR}/results.txt

## End of file
