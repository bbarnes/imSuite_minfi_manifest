
EXE="/Users/bretbarnes/Documents/tools/ucsc/bigBedToBed"
TOP="/Users/bretbarnes/Documents/data/workhorse/annotation/Homo_sapiens/hg19/snp"

TAR="dbSnp153"
TAR="dbSnp153Common"

BIG="${TOP}/${TAR}.bb"
BED="${TOP}/${TAR}.bed"

echo "Starting..."
echo "  TOP='${TOP}'"
echo "  EXE='${EXE}'"
echo
echo "  TAR='${TAR}'"
echo "  BIG='${BIG}'"
echo "  BED='${BED}'"
echo

#
# Run bigBedToBed Conversion::
#
CMD="${EXE} ${BIG} ${BED}"
echo "  CMD='${CMD}' ... "
eval $CMD
echo

#
# Run Bgzip Compression::
#
CMD="bgzip -f ${BED}"
echo "  CMD='${CMD}' ... "
eval $CMD
echo

#
# Run Tabix Indexing::
#
CMD="tabix -p bed -f ${BED}.gz"
echo "  CMD='${CMD}' ... "
eval $CMD
echo

echo "done."
echo

# End of file
