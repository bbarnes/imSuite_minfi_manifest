#!/bin/sh

TOP="/Users/bretbarnes/Documents"
NAM="GRCh37.12032021_improbe-designOutput"
DAT="${TOP}/data/improbe/scratch/cgnDB/dbSNP_Core4/design-output/${NAM}.tsv.gz"

OUT_DIR="${TOP}/data/improbe/scratch/cgnDB/dbSNP_Core4/design-output/subset"
OUT_NAM="${OUT_DIR}/${NAM}.subset"
COL_TSV="${OUT_NAM}.cols.tsv"
POS_BED="${OUT_NAM}.pos-sorted.bed"
CGN_BED="${OUT_NAM}.cgn-sorted.bed"

echo "Starting..."
echo "  TOP='${TOP}'"
echo "  DAT='${DAT}'"
echo
echo "  OUT_DIR='${OUT_DIR}'"
echo "  OUT_NAM='${OUT_NAM}'"
echo "  COL_TSV='${COL_TSV}'"
echo "  POS_BED='${POS_BED}'"
echo "  CGN_BED='${CGN_BED}'"
echo

mkdir -p ${OUT_DIR}

# Write Subset Column Headers::
# cut -f 1,2,4,5,10,16,20,21,22,23,25,40,45,48,57 > ${COL_TSV}

# Skip::
#
# gzip -dc ${DAT} | head -n 1 |
#     perl -pe '@d=split(/\t/,$_); print "$d[3]\t$d[4]\t$d[4]\t$d[0]\tProbe_Type\t$d[1]\t$d[9]\t$d[15]\t$d[19]\t$d[20]\t$d[21]\t$d[22]\t$d[24]\t$d[39]\t$d[44]\t$d[47]\t$d[56]" ' > ${COL_TSV}

# Skip::
#
# gzip -dc ${DAT} | \
#     # head | \
#     grep -v "^Seq_ID" | \
#     perl -pe '@d=split(/\t/,$_); s/^.*$//; $d[0]=~ s/^cg//; $d[0]=~ s/^0+//; $p=$d[4]-1; $d[21]=substr($d[21],0,1); print "chr$d[3]\t$p\t$d[4]\t$d[0]\tcg\t$d[1]\t$d[9]\t$d[15]\t$d[19]\t$d[20]\t$d[21]\t$d[22]\t$d[24]\t$d[39]\t$d[44]\t$d[47]\t$d[56]" ' > ${POS_BED}

echo "Done Subsetting!!!"
echo

# Skip::
#
# bgzip -f ${POS_BED}
POS_BED="${POS_BED}.gz"

echo "Done Bgzip Compression:: POS!!!"
echo

# Skip::
#
# tabix -p bed -f ${POS_BED}

echo "Done Tabix-Bed Indexing!!!"
echo

gzip -dc ${POS_BED} | sort -n -k 4,4 > ${CGN_BED}

echo "Done Sorting by CGN!!!"
echo

bgzip -f ${CGN_BED}
CGN_BED="${CGN_BED}.gz"

echo "Done Bgzip Compression:: CGN!!!"
echo

tabix -b 4 -e 4 -f -s 5 ${CGN_BED}

echo "Done Tabix-Bed Indexing!!!"
echo

# tabix /Users/bretbarnes/Documents/data/improbe/scratch/cgnDB/dbSNP_Core4/design-output/subset/GRCh37.12032021_improbe-designOutput.subset.pos-sorted.bed.gz chr1:10471-10471
# tabix /Users/bretbarnes/Documents/data/improbe/scratch/cgnDB/dbSNP_Core4/design-output/subset/GRCh37.12032021_improbe-designOutput.subset.cgn-sorted.bed.gz cg:2-2

# End of file
