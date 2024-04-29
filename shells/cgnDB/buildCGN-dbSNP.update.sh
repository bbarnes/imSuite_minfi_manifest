#!/bin/sh

# Original build as of 20022021 with Top Sequence Output::
#  RUN=GRCh36-38.mm10.FWD_SEQ

VER=10

SRC_MAX=4
CHR_MAX=2
LEN_MAX=10000000

# New build with Forward Sequences::
# RUN=GRCh36-GRCh37-GRCh38-GRCm10.12032021
# RUN=GRCh36-GRCh37-GRCh38-GRCm10.13032021
# RUN=GRCh36-GRCh37-GRCh38-GRCm10.14032021_test
RUN=dbSNP_update_test

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
else
    echo "Failed to locate top dir!"
    exit
fi

# TOP=/illumina/scratch/darkmatter/Projects/dbCGN
# TOP="/Users/bretbarnes/Documents"
# EXE=${TOP}/tools/scripts/cgnDB/latest/buildCGN-BSCdb.cluster.pl
# EXE=${TOP}/tools/scripts/cgnDB/latest/buildCGN-dbSNP.update.pl

EXE="${TOP}/tools/Infinium_Methylation_Workhorse/scripts/perl/cgnDB/buildCGN-dbSNP.update.pl"
OUT="${TOP}/scratch/cgnDB/${RUN}/design-input"

#
# Need an option for reading in csv version...
#
# CGN=${TOP}/manifests/manifest.cgn.hash
# CGN=/illumina/scratch/darkmatter/dbCGN/manifests/manifest.cgn.hash
# CGN=${TOP}/data/cgnDB/GRCh36-38.mm10/pre-assigned.cgnTop.hash.csv.gz
#
# CAN="${TOP}/data/cgnDB/GRCh36-38.mm10/canonical-12032021.cgn-top-src.hash.csv.gz"
#
# CAN="${TOP}/data/improbe/dbCGN/pre-assigned/pre-assigned.cgnTop.hash.csv.gz"
# CAN="${TOP}/data/improbe/dbCGN/pre-assigned/canonical-assignment.cgn-top-grp.csv.gz"
CAN="${TOP}/tools/Infinium_Methylation_Workhorse/dat/manifest/cgnDB/canonical-assignment.cgn-top-grp.csv.gz"

# GEN38=/illumina/scratch/darkmatter/iGenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/GRCh38.genome.fa.gz
# GEN37=/illumina/scratch/darkmatter/iGenomes/Homo_sapiens/NCBI/GRCh37/Sequence/WholeGenomeFasta/GRCh37.genome.fa.gz
# GEN36=/illumina/scratch/darkmatter/iGenomes/Homo_sapiens/NCBI/GRCh36/Sequence/WholeGenomeFasta/GRCh36.genome.fa.gz
# GEN10=/illumina/scratch/darkmatter/iGenomes/Mus_musculus/NCBI/GRCm10/Sequence/WholeGenomeFasta/GRCm10.genome.fa.gz

GEN38="${GEN}/iGenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/GRCh38.genome.fa.gz"
GEN37="${GEN}/iGenomes/Homo_sapiens/NCBI/GRCh37/Sequence/WholeGenomeFasta/GRCh37.genome.fa.gz"
GEN36="${GEN}/iGenomes/Homo_sapiens/NCBI/GRCh36/Sequence/WholeGenomeFasta/GRCh36.genome.fa.gz"
GEN10="${GEN}/iGenomes/Mus_musculus/NCBI/GRCm10/Sequence/WholeGenomeFasta/GRCm10.genome.fa.gz"

mkdir -p ${OUT}

${EXE} -v ${VER} \
       -out ${OUT} \
       -can ${CAN} \
       -mod a -name ${RUN} \
       -src GRCh38 -f ${GEN38} \
       -src GRCh37 -f ${GEN37} \
       -src GRCh36 -f ${GEN36} \
       -src GRCm10 -f ${GEN10} \
       -srcMax ${SRC_MAX} \
       -chrMax ${CHR_MAX} \
       -lenMax ${LEN_MAX}

echo "done."

## end of file
