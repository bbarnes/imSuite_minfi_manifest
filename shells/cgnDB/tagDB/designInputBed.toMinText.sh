

TOP="/Users/bretbarnes/Documents/data/improbe/scratch/cgnDB/dbSNP_Core4/design-input"
TOP="/illumina/scratch/darkmatter/scratch/cgnDB/dbSNP_Core4/design-input"

BLD="GRCm10"
BLD="GRCh37"

DIR="${TOP}/min"
INP="${TOP}/${BLD}.cgn.bed.gz"
OUT="${DIR}/${BLD}.cgn.min.txt.gz"

echo "TOP = ${TOP}"
echo "DIR = ${DIR}"
echo "INP = ${INP}"
echo "OUT = ${OUT}"
echo ""

mkdir -p ${DIR}

gzip -dc ${INP} | \
    cut -f 1,3,4,6 | \
    perl -pe 's/^chr//' | \
    gzip -c -> ${OUT}


echo "done."

# End of file
