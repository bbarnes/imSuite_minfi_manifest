#!/bin/bash

if [ "$#" -lt 2 ]; then
    echo "Usage: $0 datPath outPath"
    exit 1
fi

DAT=$1
OUT=$2
# DAT="/Users/bretbarnes/Documents/scratch/cgnDB/dbSNP_update_test/design-input/GRCh38_GRCh37_GRCh36_GRCm10.cgn-set.csv.gz"

TOP_LIX="/illumina/scratch/darkmatter"
TOP_MAC="/Users/bretbarnes/Documents"

if [ -d ${TOP_MAC} ]
then
    TOP="${TOP_MAC}"
    GEN="${TOP_MAC}/data"
elif [ -d ${TOP_LIX} ]
then
    TOP="${TOP_LIX}"
    GEN="${TOP_LIX}"
    
    EXE="${TOP_LIX}/bin/improbe"
    TAN="${TOP_LIX}/dat/Tango_A_or_B_11mer_s1.dat"
    MER="${TOP_LIX}/dat/human-36.1-methyl-and-unmethyl-13mer-s3-for-infinium-methylation.dat"
    
else
    echo "Failed to locate top dir!"
    exit
fi

mkdir -p ${OUT}

INP="${OUT}/improbe-cgn-template-seqs.tsv"
DES="${OUT}/improbe-cgn-probe-designs.tsv.gz"

echo "OUT=${OUT}"
echo "DAT=${DAT}"
echo "INP=${INP}"

# First echo the header to the output::
HEAD="${TOP}/tools/Infinium_Methylation_Workhorse/dat/manifest/cgnDB/improbe-input-header.tsv"
cp ${HEAD} ${INP}

gzip -dc ${DAT} | \
    perl -pe 's/\n$//; @d=split(/,/,$_); s/^(.*)$//; print "$d[0]_$d[2]_$d[3]_$d[4]_$d[5]\t".substr($d[1],0,60)."[",substr($d[1],60,2)."]".substr($d[1],62)."\tA\t1\t$d[0]\tFALSE\n"' \
	 head \
	 >> ${INP}

# Run improbe::
${EXE} -t Both -c Both -n ${MER} -a ${TAN} -V | gzip -c - > ${DES}

echo "done."

# End of file
