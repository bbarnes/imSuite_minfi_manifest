

VER=10
SRC="GRCh37"
# NAM=".dbSNP-151.iupac"
NAM=""

TOP="/Users/bretbarnes/Documents"
EXE="${TOP}/tools/scripts/cgnDB/latest/bsc-fasta.pl"
DAT="${TOP}/data/iGenomes/Homo_sapiens/NCBI/GRCh37/Sequence/WholeGenomeFasta"

# FAS="${DAT}/${SRC}.dbSNP151-genome.fa.gz"
FAS="${DAT}/${SRC}.genome.fa.gz"
FAS="${DAT}/${SRC}.genome.dbSNP-151.iupac.fa.gz"

FAS="${DAT}/${SRC}.genome${NAM}.fa.gz"


echo "SRC=${SRC}"
echo "NAM=${NAM}"
echo ""
echo "EXE=${EXE}"
echo "DAT=${DAT}"
echo "FAS=${FAS}"
echo ""

# CMD="${EXE} -v ${VER} -o ${DAT} -s ${SRC} -n ${NAM} -fas ${FAS}"
CMD="${EXE} -v ${VER} -o ${DAT} -s ${SRC} -fas ${FAS}"

echo "CMD=${CMD}"
eval $CMD

echo ""
echo "done converting fasta."

# End of file
