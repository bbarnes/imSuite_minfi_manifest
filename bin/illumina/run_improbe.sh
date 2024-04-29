#!/bin/sh

if [ "$#" -lt 2 ]; then
    echo "Usage: $0 input.tsv.gz outName"
    exit 1
fi

CGN=$1
NAM=$2

# mkdir -p ${OUT}
# LOG=${OUT}/${NAM}.log
# DES=${OUT}/${NAM}.tsv.gz

LOG=/output/${NAM}.improbe-designOutput.log
DES=/output/${NAM}.improbe-designOutput.tsv.gz

echo "CGN="${CGN}
# echo "OUT="${OUT}
echo "NAM="${NAM}

echo "LOG="${LOG}
echo "DES="${DES}

gzip -dc /work/${CGN} | \
    /repo/bin/improbe \
	-tBoth -cBoth \
	-n/repo/bin/human-36.1-methyl-and-unmethyl-13mer-s3-for-infinium-methylation.dat \
	-a/repo/bin/Tango_A_or_B_11mer_s1.dat \
	-V - 2> ${LOG} | \
    gzip -c -> ${DES}

echo "done."

exit

gzip -dc . | \
    /repo/bin/improbe \
	-tBoth -cBoth \
	-n/repo/bin/human-36.1-methyl-and-unmethyl-13mer-s3-for-infinium-methylation.dat \
	-a/repo/bin/Tango_A_or_B_11mer_s1.dat \
	-V - 2> /output/improbe-designOutput.log | \
    gzip -c -> /output/improbe-designOutput.tsv.gz


# End of file
