#!/bin/bash

#$ -cwd
#$ -N sesame_testing
#$ -q whole
#$ -pe threaded 28


PROGRAM="Workhorse-Unstained/scripts/R/triplecrown/triplecrown.R"

DOC_VER="v.2.10"
DOC_NAME="workhorse-unstained"
DOC_IMAGE="bbarnesimdocker/im_workhorse:${DOC_NAME}.${DOC_VER}"
DOC_SHELL="run_cmd.sh"

VERBOSE="100"

TOP_DIR_A="/illumina/scratch/darkmatter"
TOP_DIR_B="/Users/bretbarnes/Documents"

# REF_SOURCE="UCSC"
REF_SOURCE="NCBI"

REF_SPECIES="Homo_sapiens,Homo_sapiens,Homo_sapiens,Mus_musculus"

if [ ${REF_SOURCE} == "NCBI" ]; then
    
    REF_FILE="hg38.fa.gz,hg19.fa.gz,hg18.fa.gz,mm10.fa.gz"
    REF_BUILD="hg38,hg19,hg18,mm10"

elif [ ${REF_SOURCE} == "UCSC" ]; then

    REF_FILE="GRCh38.genome.fa.gz,GRCh37.genome.fa.gz,GRCh36.genome.fa.gz,GRCm38.genome.fa.gz"
    REF_BUILD="GRCh38,GRCh37,GRCh36,GRCm38"
    
else
    echo "Unrecognized Reference Source ${REF_SOURCE}!"
    exit
fi

echo "Reference Source = ${REF_SOURCE}"
echo

RUN_NAME=${REF_SOURCE}

if [ -e ${TOP_DIR_A} ]; then
    TOP_DIR=${TOP_DIR_A}

    SIG_IMAGE="/illumina/scratch/darkmatter/docker/software/${DOC_NAME}.${DOC_VER}.sif"
    SCRATCH_ROOT="/scratch"

    echo "Singularity Image"
    echo "  SIG_IMAGE = ${SIG_IMAGE}"

elif [ -e ${TOP_DIR_B} ]; then
    TOP_DIR=${TOP_DIR_B}

    echo "Docker Image"
    echo "  DOC_IMAGE = ${DOC_IMAGE}"
    
else
    echo "Unrecognized Rscript EXE!"
    exit
fi


DAT_PATH="/repo/Workhorse-Unstained/dat"
OUT_PATH="${TOP_DIR}/data/scratch/docker-${DOC_VER}"
GEN_PATH="${TOP_DIR}/data/imGenomes"

REF_PATH="/genome"
CHR_PATH="/genome"

mkdir -p ${OUT_PATH}
mkdir -p ${GEN_PATH}
    
CANONICAL_CSV="${DAT_PATH}/auxiliary/canonical_cgn_top_grp.csv.gz"
TEMP_OFF_CSV="${DAT_PATH}/auxiliary/triplecrown/template_offset.csv.gz"

RELOAD="1"

PARALLEL="--parallel"
PARALLEL=""

# This isn't an option yet...
# RCPP="--rcpp"
RCPP=""

TIME_TRACK=""
TIME_TRACK="--track_time"

echo "Local Mounted Directoreis::"
echo "  gen_path = ${GEN_PATH}"
echo "  out_path = ${OUT_PATH}"
echo
echo "  chr_path = ${CHR_PATH}"
echo "  ref_path = ${REF_PATH}"
echo "  ref_file = ${REF_FILE}"
echo "     build = ${REF_BUILD}"
echo "   species = ${REF_SPECIES}"
echo 
echo "  canonical = ${CANONICAL_CSV}"
echo "  temp_off  = ${TEMP_OFF_CSV}"
echo
echo "  reload     = ${RELOAD}"
echo "  parallele  = ${PARALLEL}"
echo "  rcpp       = ${RCPP}"
echo "  time_track = ${TIME_TRACK}"
echo "  verbose    = ${VERBOSE}"
echo

#
# Build Run Command::
#
OPT_STR="${PROGRAM} \
--run_name=${RUN_NAME} \
--out_path=/output \
--chr_path=${CHR_PATH} \
--ref_path=${REF_PATH} \
--ref_file=${REF_FILE} \
--ref_build=${REF_BUILD} \
--ref_species=${REF_SPECIES} \
--canonical_csv=${CANONICAL_CSV} \
--temp_off_csv=${TEMP_OFF_CSV} \
--reload=${RELOAD} \
${PARALLEL} ${RCPP} ${TIME_TRACK} \
--verbose=${VERBOSE} "

# If you just want to verify docker works::
#  print the help command::
#
# OPT_STR="${PROGRAM} --help "

if [ -e ${TOP_DIR_A} ]; then

    echo "Singluarity (Cluster) Definitions::"
    echo "  SIG_IMAGE = ${SIG_IMAGE}"
    echo "  DOC_SHELL = ${DOC_SHELL}"
    echo "    OPT_STR = ${OPT_STR}"
    echo

    RUN_CMD="singularity exec \
-B ${OUT_PATH}:/output \
-B ${GEN_PATH}:/genome \
-B ${SCRATCH_ROOT}:${SCRATCH_ROOT} \
--workdir ${SCRATCH_ROOT} \
${SIG_IMAGE} ${DOC_SHELL} ${OPT_STR}"

elif [ -e ${TOP_DIR_B} ]; then

    echo "Local Docker (Mac) Definitions::"
    echo "  DOC_IMAGE = ${DOC_IMAGE}"
    echo "  DOC_SHELL = ${DOC_SHELL}"
    echo "    OPT_STR = ${OPT_STR}"
    echo

    RUN_CMD="docker run -i --rm -w /work \
-v ${OUT_PATH}:/output \
-v ${GEN_PATH}:/genome \
${DOC_IMAGE} ${DOC_SHELL} ${OPT_STR}"

else
    echo "Unrecognized Rscript EXE!"
    exit
fi

echo
echo "RUN_CMD = ${RUN_CMD}"
echo

# ${RUN_CMD}

echo "Done: $0"

# -v ${CHR_H38_PATH}:/chr_h38 
# -v ${CHR_H37_PATH}:/chr_h37 
# -v ${CHR_H36_PATH}:/chr_h36 
# -v ${CHR_M38_PATH}:/chr_m38 
# -v ${REF_H38_PATH}:/ref_h38 
# -v ${REF_H37_PATH}:/ref_h37 
# -v ${REF_H36_PATH}:/ref_h36 
# -v ${REF_M38_PATH}:/ref_m38 

# End of file
