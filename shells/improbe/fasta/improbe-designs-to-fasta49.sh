#!/bin/sh

TOP="/Users/bretbarnes/Documents"
NAM="GRCh37.12032021_improbe-designOutput"
DAT="${TOP}/data/improbe/scratch/cgnDB/dbSNP_Core4/design-output/subset/${NAM}.subset.cgn-sorted.bed.gz"

OUT_DIR="${TOP}/data/improbe/scratch/cgnDB/dbSNP_Core4/design-output/fasta"
OUT_NAM="${OUT_DIR}/${NAM}"

U49_FAS="${OUT_NAM}.U49.fa"
M49_FAS="${OUT_NAM}.M49.fa"
TMP_TSV="${OUT_NAM}.tmp.tsv"
UNQ_TSV="${OUT_NAM}.unq-sorted.tsv"

echo "Starting..."
echo "  TOP='${TOP}'"
echo "  DAT='${DAT}'"
echo
echo "  OUT_DIR='${OUT_DIR}'"
echo "  OUT_NAM='${OUT_NAM}'"
echo
echo "  U49_FAS='${U49_FAS}'"
echo "  M49_FAS='${M49_FAS}'"
echo "  TMP_TSV='${TMP_TSV}'"
echo

mkdir -p ${OUT_DIR}

# Write Subset Column Headers::
gzip -dc ${DAT} | \
    # head -n 20000 | \
    cut -f 4,8,11,12,16 > ${TMP_TSV}

sort -u ${TMP_TSV} | sort -n -k 1,1 > ${UNQ_TSV}

echo "done."
echo

exit

cat /Users/bretbarnes/Documents/data/improbe/scratch/cgnDB/dbSNP_Core4/design-output/fasta/GRCh37.12032021_improbe-designOutput.unq-sorted.tsv | \
    head | \
    perl -pe '@d=split(/\t/,$_); s/^\s+$//; $d[4]=~ s/\s+$//; $e1=substr($d[1],-1); s/^.*$//; print ">$d[0]"."_$d[2]"."$d[3]"."$e1_M\n$d[1]' | less
 

    # perl -pe '@d=split(/\t/,$_); $e1=substr($d[1],-1); $d[1] =~ s/[A-Z]$//; $e2=substr($d[4],-1); $d[4] =~ s/[A-Z]$//; s/^.*$//; s/\s+$//; print join("\t", $d[0],$d[2],$d[3],"M\t$e1",$d[1])."\n".join("\t", $d[0],$d[2],$d[3],"U\t$e2",$d[4]);' | \
    # sort -n -k 1,5 | \
    # sort -u | \
    gzip -c - > ${TMP_TSV}.gz


    cut -f 1,2,4,5,10,16,20,21,22,23,25,40,45,48,57


gzip -dc ${DAT} | \
    head | \
    cut -f 4,8,9,10,11,12,16

#    cut -f 4,10,11,13

exit

    sort -n -u -k 1,1 | \
    perl '@d=split(/\t/,$_); s/^.*$//; $d[3]=substr($d[3],-1); print ">$d[0]"."_$d[1]"."_$d[2]_U\n$d[3]"' | \
    gzip -c - > ${U49_FAS}

echo "Done Writing Fasta:: U!!!"
echo

#
# Run bsmap::
#


# End of file
