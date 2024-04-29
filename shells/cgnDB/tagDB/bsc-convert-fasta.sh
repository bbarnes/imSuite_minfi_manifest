

VER=10
BLD="GRCh37"
SPE="Homo_sapiens"

# This allows the switch::
# NAM=".dbSNP-151.iupac"
NAM=""

# TOP="/Users/bretbarnes/Documents"
TOP="/illumina/scratch/darkmatter"
EXE="${TOP}/tools/scripts/cgnDB/latest/bsc-fasta.pl"
DAT="${TOP}/data/imGenomes/${SPE}/NCBI/${BLD}/Sequence/WholeGenomeFasta"

# FAS="${DAT}/${BLD}.dbSNP151-genome.fa.gz"
# FAS="${DAT}/${BLD}.genome.fa.gz"
# FAS="${DAT}/${BLD}.genome.dbSNP-151.iupac.fa.gz"

FAS="${DAT}/${BLD}.genome${NAM}.fa.gz"

echo "BLD=${BLD}"
echo "NAM=${NAM}"
echo ""
echo "EXE=${EXE}"
echo "DAT=${DAT}"
echo "FAS=${FAS}"
echo ""

# Need to swith these commands and the NAM defintion above...
#
# CMD="${EXE} -v ${VER} -o ${DAT} -s ${BLD} -n ${NAM} -fas ${FAS}"
CMD="${EXE} -v ${VER} -o ${DAT} -s ${BLD} -fas ${FAS}"

echo "CMD=${CMD}"
eval $CMD

echo ""
echo "done converting fasta."

# End of file
