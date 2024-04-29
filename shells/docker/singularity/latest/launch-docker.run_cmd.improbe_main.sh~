#!/bin/bash

#$ -cwd
#$ -N sesame_testing
#$ -q whole
#$ -pe threaded 28

if [ "$#" -lt 1 ]; then
    echo "Usage: $0 run_name [ -l ] [ -p program_name ] [options]"
    exit 1
fi

echo "Args = $@"

RUN_NAME=$1
shift

LOCAL=false

PRG_NAME="ewas_loci_variation.R"
PRG_NAME="ewas_swifthoof.R"
RUN_LOCAL=false
USR_OPTS=""

while [ -n "$1" ]; do

    # $1 = `echo $1 | sed -e 's/^[[:space:]]*//'`

    case "$1" in
        -l) RUN_LOCAL=true ;;
        -p) PRR_NAME="$2"
            shift ;;
        *) USR_OPTS="${USR_OPTS} $1" ;;

        esac
        shift
done

PROGRAM="Workhorse-Unstained/scripts/R/swifthoof/${PRG_NAME}"

DOC_VER="v.2.21"
DOC_NAME="workhorse-unstained"
DOC_IMAGE="bbarnesimdocker/im_workhorse:${DOC_NAME}.${DOC_VER}"
DOC_SHELL="run_cmd.sh"

VERBOSE="100"
RELOAD="1"

SINGLE="--single"
SINGLE=""

PARALLEL="--parallel"
PARALLEL=""

# This isn't an option yet...
# RCPP="--rcpp"
RCPP=""

TIME_TRACK=""
TIME_TRACK="--track_time"

PLATFORM="EPIC"
SAMPLE_MAX="0"

# REF_SOURCE="UCSC"
REF_SOURCE="NCBI"

REF_SPECIES="Homo_sapiens,Homo_sapiens,Homo_sapiens,Mus_musculus"

TOP_DIR_A="/illumina/scratch/darkmatter"
TOP_DIR_B="/Users/bretbarnes/Documents"

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

# RUN_NAME=${REF_SOURCE}

RSCRIPT="Rscript"
if [ -e ${TOP_DIR_A} ]; then
    TOP_DIR=${TOP_DIR_A}
    
    RSCRIPT="Rscript"
    SIG_IMAGE="/illumina/scratch/darkmatter/docker/software/${DOC_NAME}.${DOC_VER}.sif"
    SCRATCH_ROOT="/scratch"

    echo "Singularity Image"
    echo "  SIG_IMAGE = ${SIG_IMAGE}"

elif [ -e ${TOP_DIR_B} ]; then
    TOP_DIR=${TOP_DIR_B}
    LOCAL_SDF_PATH="scratch/stable_ewas_swifthoof_scratch/sesame-UCSC-v3/prefix_to_sdf"

    RSCRIPT="Rscript"

    echo "Docker Image"
    echo "  DOC_IMAGE = ${DOC_IMAGE}"
    
else
    echo "Unrecognized Rscript EXE!"
    exit
fi

# Found in Workhorse Source Code (git repository)
DAT_PATH="/repo/Workhorse-Unstained/dat"
CANONICAL_CSV="${DAT_PATH}/auxiliary/canonical_cgn_top_grp.csv.gz"
TEMP_OFF_CSV="${DAT_PATH}/auxiliary/triplecrown/template_offset.csv.gz"

# Define all local directory paths, then mount them and
#  set mounted file locations
#
IDAT_PATH="${TOP_DIR}/data/idats"
OUT_PATH="${TOP_DIR}/scratch/docker-${DOC_VER}"
GEN_PATH="${TOP_DIR}/data/imGenomes"

SDF_PATH="${TOP_DIR}/${LOCAL_SDF_PATH}"
MAN_PATH="${TOP_DIR}/Projects/EWAS/data/manifests"
SAM_PATH="${TOP_DIR}/Projects/EWAS/data/sample_sheets"

REF_PATH="/genome"
CHR_PATH="/genome"

MAN_CSV="/manifests/EWAS_PQC122021-NA-NA-GRCh37_sesame.beta.historic-beadPool.csv.gz"
SAM_CSV="/samplesheets/EWAS-alpha-cross-product-testing.07012022.sampleSheet.csv.gz"
    
if [ "$RUN_LOCAL" = true ] ; then
    OUT_PATH="${TOP_DIR}/scratch/docker-local"
    MAN_CSV="${MAN_PATH}/EWAS_PQC122021-NA-NA-GRCh37_sesame.beta.historic-beadPool.csv.gz"
    SAM_CSV="${SAM_PATH}/EWAS-alpha-cross-product-testing.07012022.sampleSheet.csv.gz"
fi


mkdir -p ${OUT_PATH}
mkdir -p ${GEN_PATH}
    
echo "Command Line Parameters::"
echo "  run name   = ${RUN_NAME}"
echo "  program    = ${PRG_NAME}"
echo "  local run  = ${RUN_LOCAL}"
echo "  user opts  = ${USR_OPTS}"
echo
echo "Local Mounted Directoreis::"
echo "  gen_path = ${GEN_PATH}"
echo "  top_path = ${TOP_PATH}"
echo "  out_path = ${OUT_PATH}"
echo " idat_path = ${IDAT_PATH}"
echo
echo "  platform   = ${PLATFORM}"
echo "  sample_max = ${SAMPLE_MAX}"
echo
echo "  chr_path = ${CHR_PATH}"
echo "  ref_path = ${REF_PATH}"
echo "  ref_file = ${REF_FILE}"
echo "     build = ${REF_BUILD}"
echo "   species = ${REF_SPECIES}"
echo 
echo "Hard Coded Parameters::"
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
--top_path=${TOP_DIR} \
--sdf_path=${SDF_PATH} \
--manifest=${MAN_CSV} \
--sample_sheet=${SAM_CSV} \
--platform=${PLATFORM} \
--sample_max=${SAMPLE_MAX} \
--reload=${RELOAD} \
${PARALLEL} ${RCPP} ${TIME_TRACK} \
${SINGLE} \
--verbose=${VERBOSE} \
${USR_OPTS} "

LOCAL_DIR="--out_path=${OUT_PATH} --idat_path=${IDAT_PATH}"
MOUNT_DIR="--out_path=/output --idat_path=/idats"

# If you just want to verify docker works::
#  print the help command::
#
# OPT_STR="${PROGRAM} --help "

if [ -e ${TOP_DIR_A} ]; then

    echo "Singluarity (Cluster) Definitions::"
    echo "  SIG_IMAGE = ${SIG_IMAGE}"
    echo "  DOC_SHELL = ${DOC_SHELL}"
    echo "  MOUNT_DIR = ${MOUNT_DIR}"
    echo "    OPT_STR = ${OPT_STR}"
    echo

    RUN_CMD="singularity exec \
-B ${OUT_PATH}:/output \
-B ${GEN_PATH}:/genome \
-B ${MAN_PATH}:/manifests \
-B ${SAM_PATH}:/samplesheets \
-B ${IDAT_PATH}:/idats \
-B ${SCRATCH_ROOT}:${SCRATCH_ROOT} \
--workdir ${SCRATCH_ROOT} \
${SIG_IMAGE} ${DOC_SHELL} ${OPT_STR} ${MOUNT_DIR}"

elif [ -e ${TOP_DIR_B} ]; then

    echo "Local Docker (Mac) Definitions::"
    echo "  DOC_IMAGE = ${DOC_IMAGE}"
    echo "  DOC_SHELL = ${DOC_SHELL}"
    echo "  MOUNT_STR = ${MOUNT_DIR}"
    echo "    OPT_STR = ${OPT_STR}"
    echo

    RUN_CMD="docker run -i --rm -w /work \
-v ${OUT_PATH}:/output \
-v ${GEN_PATH}:/genome \
-v ${MAN_PATH}:/manifests \
-v ${SAM_PATH}:/samplesheets \
-v ${IDAT_PATH}:/idats \
${DOC_IMAGE} ${DOC_SHELL} ${OPT_STR} ${MOUNT_DIR}"

else
    echo "Unrecognized Rscript EXE!"
    exit
fi

if [ "$RUN_LOCAL" = true ]; then
    LOCAL_CMD="${RSCRIPT} ${OPT_STR} ${LOCAL_DIR}"
    echo
    echo "RUN_LOCAL_CMD = ${LOCAL_CMD}"
    echo

    eval $LOCAL_CMD

else
    
    echo
    echo "RUN_CMD = ${RUN_CMD}"
    echo

    eval $RUN_CMD
    
fi

echo "Done: $0"

# --chr_path=${CHR_PATH} \
# --ref_path=${REF_PATH} \
# --ref_file=${REF_FILE} \
# --ref_build=${REF_BUILD} \
# --ref_species=${REF_SPECIES} \
# --canonical_csv=${CANONICAL_CSV} \
# --temp_off_csv=${TEMP_OFF_CSV} \

# -v ${CHR_H38_PATH}:/chr_h38 
# -v ${CHR_H37_PATH}:/chr_h37 
# -v ${CHR_H36_PATH}:/chr_h36 
# -v ${CHR_M38_PATH}:/chr_m38 
# -v ${REF_H38_PATH}:/ref_h38 
# -v ${REF_H37_PATH}:/ref_h37 
# -v ${REF_H36_PATH}:/ref_h36 
# -v ${REF_M38_PATH}:/ref_m38 

# End of file
