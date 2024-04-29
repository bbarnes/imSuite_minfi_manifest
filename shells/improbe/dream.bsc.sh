#!/bin/bash

TOP_DIR="/Users/bbarnes/mylocal/Documents"
OUT_DIR="${TOP_DIR}/scratch/dream-designs"
DAT_CSV="${TOP_DIR}/data/manifests/methylation/GenomeStudio/EX100K_GRCh37_1F.body-section.genomic-sorted.csv.gz"

SPE="Homo_Sapiens"
MAX=1
EXE="${TOP_DIR}/tools/Workhorse-Unstained/scripts/perl/bsc/bsc.dream.pl"

mkdir -p ${OUT_DIR}

# Old Testing Command::
# CMD="gzip -dc ${DAT_CSV} | cut -d, -f 21 | ${EXE} -max ${MAX}"

# New Testing Command::
# CMD="${EXE} -s manifest -m ${DAT_CSV} -o ${OUT_DIR} -species ${SPE} -max ${MAX}"
# CMD="${EXE} -s manifest -m ${DAT_CSV} -o ${OUT_DIR} -species ${SPE} -max ${MAX}"
CMD="${EXE} -s manifest -m ${DAT_CSV} -o ${OUT_DIR} -species ${SPE} -v 10 "

echo "CMD='${CMD}'"
echo

$CMD

# exec($CMD)
# exec "${CMD}"
# exec "${CMD}"

echo
echo "Done"
echo

# End of file
