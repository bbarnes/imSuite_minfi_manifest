
HEAD="/Users/bretbarnes/Documents/data/improbe/scratch/cgnDB/dbSNP_Core4/design-output/improbe-cgn-probe-designs.sort-unique.header.tsv"
DATA="/Users/bretbarnes/Documents/data/improbe/scratch/cgnDB/dbSNP_Core4/design-output/improbe-cgn-probe-designs.sort-unique.tsv"
OUT="/Users/bretbarnes/Documents/data/improbe/scratch/cgnDB/dbSNP_Core4/design-output/improbe-cgn-probe-designs.sort-unique.wHeader.tsv.gz"

cat ${HEAD} ${DATA} | perl -pe 's/\[//; s/\]//;' | gzip -c - > ${OUT}

# End of file
