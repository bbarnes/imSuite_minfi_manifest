## Example of how to run simluation analysis pipeline

BIN=/home/bbarnes/vcf-compare/scripts

if [ "$#" -ne 4 ]; then
    echo "Usage: $0 candidate.vcf reference1.vcf reference2.vcf outputDir"
    exit 1
fi 

CAN_VCF=$1
REF1_VCF=$2
REF2_VCF=$3
OUT_DIR=$4

if [ ! -d $OUT_DIR ]; then
    mkdir $OUT_DIR
fi

cat ${CAN_VCF} | ${BIN}/vcfToBedSplit.pl -o ${OUT_DIR}/candidate
cat ${REF1_VCF} | ${BIN}/vcfToBedSplit.pl -o ${OUT_DIR}/reference1
cat ${REF2_VCF} | ${BIN}/vcfToBedSplit.pl -o ${OUT_DIR}/reference2

${BIN}/windowBed -w 10 \
    -a ${OUT_DIR}/candidate.bed \
    -b ${OUT_DIR}/reference1.bed \
    >  ${OUT_DIR}/intersect1.txt

${BIN}/windowBed -w 10 \
    -a ${OUT_DIR}/candidate.bed \
    -b ${OUT_DIR}/reference2.bed \
    >  ${OUT_DIR}/intersect2.txt

${BIN}/vcfCompare.pl \
    -c ${OUT_DIR}/candidate.vcf \
    -r ${OUT_DIR}/reference1.vcf \
    -i ${OUT_DIR}/intersect1.txt \
    -r ${OUT_DIR}/reference2.vcf \
    -i ${OUT_DIR}/intersect2.txt \
    -o ${OUT_DIR}/results \
    >  ${OUT_DIR}/results.txt

## End of file
