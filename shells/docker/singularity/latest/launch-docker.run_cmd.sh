#!/bin/bash

#$ -cwd
#$ -N sesame_testing
#$ -q whole
#$ -pe threaded 28

PROGRAM="Infinium_Methylation_Workhorse/scripts/R/workhorse/workhorse_main.R"

DOC_VER="v.1.49"
DOC_NAME="Infinium_Methylation_Workhorse_Centos"
DOC_IMAGE="bbarnesimdocker/im_workhorse:${DOC_NAME}.${DOC_VER}"
DOC_SHELL="run_cmd.sh"

VERBOSE="5"
SPECIES="Homo_sapiens"

TOP_DIR_A="/illumina/scratch/darkmatter"
TOP_DIR_B="/Users/bretbarnes/Documents"

ORD_CSV="McMaster_CpG_DesignFile_v4.csv.gz"
MAT_TSV="20532820_probes1.match.gz,20532820_probes2.match.gz"
AQP_TSV="20051339_A_ProductQC.txt.gz"

RELOAD="--reload"
RELOAD=""

PARALLEL="--parallel"
PARALLEL=""

OPT_STR="${PROGRAM} \
  --Rscript='Rscript' \
  --bsmap_exe='${BSP_EXE}' \
  --genome_build=GRCh37 \
  --platform=MCM \
  --version=v1 \
  --run_name=McMaster10Kselection \
  --ord_csv=${ORD_CSV} \
  --mat_tsv=${MAT_TSV} \
  --aqp_tsv=${AQP_TSV} \
  --canonical_cgn_csv=canonical.cgn-top-grp.csv.gz \
  ${RELOAD} ${PARALLEL} --verbose=${VERBOSE} "

# --memory-swap="[memory_limit]"
# --memory=${MEM}
MEM="16g"

if [ -e ${TOP_DIR_A} ]; then
    TOP_DIR=${TOP_DIR_A}

    SIG_IMAGE="/illumina/scratch/darkmatter/docker/software/${DOC_NAME}.${DOC_VER}.sif"

    SCRATCH_ROOT="/scratch"

    # MAN_LDIR="/illumina/scratch/darkmatter/tools/Infinium_Methylation_Workhorse/dat/manifest/core"
    # MAN_SDIR="-B ${MAN_LDIR:/tmp}"

    IMP_DIR="${TOP_DIR}/data/cgnDB"
    SEQ_DIR="${TOP_DIR}/scratch/improbe/cgnDB/dbSNP_Core4/prbs-P49-split-3.64"
    POS_DIR="${TOP_DIR}/scratch/cgnDB/dbSNP_Core4/design-input/min"
    CAN_DIR="${TOP_DIR}/data/cgnDB/canonical"

    OUT_DIR="${TOP_DIR}/data/scratch/docker-${DOC_VER}"

    echo "Singularity Image"
    echo "  SIG_IMAGE = ${SIG_IMAGE}"

elif [ -e ${TOP_DIR_B} ]; then
    TOP_DIR=${TOP_DIR_B}

    IMP_DIR="${TOP_DIR}/data/improbe"
    SEQ_DIR="${IMP_DIR}/scratch/cgnDB/dbSNP_Core4/design-output/prbs-p49-split"
    POS_DIR="${IMP_DIR}/scratch/cgnDB/dbSNP_Core4/design-input/min"
    CAN_DIR="${IMP_DIR}/data/improbe/scratch/cgnDB/dbSNP_Core4/design-input"

    OUT_DIR="${TOP_DIR}/scratch/docker-${DOC_VER}"

else
    echo "Unrecognized Rscript EXE!"
    exit
fi

ORD_DIR="${TOP_DIR}/data/CustomContent/McMaster/McMaster10Kselection/AQP.v2/order"
MAT_DIR="${TOP_DIR}/data/CustomContent/McMaster/McMaster10Kselection/AQP.v2/match"
AQP_DIR="${TOP_DIR}/data/CustomContent/McMaster/McMaster10Kselection/AQP.v2/aqp"

GEN_DIR="${TOP_DIR}/data/imGenomes/${SPECIES}/NCBI"
MAN_DIR="${TOP_DIR}/data/manifests"
ANN_DIR="${TOP_DIR}/data/annotation"

# If you just want to verify docker works::
#  print the help command::
#
# OPT_STR="${PROGRAM} --help "

mkdir -p ${OUT_DIR}

echo "Local Mounted Directoreis::"
echo "  OUT_DIR = ${OUT_DIR}"
echo
echo "  ORD_DIR = ${ORD_DIR}"
echo "  MAT_DIR = ${MAT_DIR}"
echo "  AQP_DIR = ${AQP_DIR}"
echo
echo "  GEN_DIR = ${GEN_DIR}"
echo "  MAN_DIR = ${MAN_DIR}"
echo "  IMP_DIR = ${IMP_DIR}"
echo "  ANN_DIR = ${ANN_DIR}"
echo
echo "  SEQ_DIR = ${SEQ_DIR}"
echo "  POS_DIR = ${POS_DIR}"
echo "  CAN_DIR = ${CAN_DIR}"
echo

echo "Input file names (order, match, aqp/pqc)::"
echo "  ORD_CSV = ${ORD_CSV}"
echo "  MAT_TSV = ${MAT_TSV}"
echo "  AQP_TSV = ${AQP_TSV}"
echo

if [ -e ${TOP_DIR_A} ]; then

    echo "Singluarity Definitions::"
    echo "  SIG_IMAGE = ${SIG_IMAGE}"
    echo "  DOC_SHELL = ${DOC_SHELL}"
    echo

    RUN_CMD="singularity exec \
	 -B ${ORD_DIR}:/order \
	 -B ${MAT_DIR}:/match \
	 -B ${AQP_DIR}:/aqp \
	 -B ${OUT_DIR}:/output \
	 -B ${SEQ_DIR}:/sequence \
	 -B ${POS_DIR}:/coordinate \
	 -B ${CAN_DIR}:/canonical \
	 -B ${GEN_DIR}:/genome \
	 -B ${IMP_DIR}:/improbe \
	 -B ${ANN_DIR}:/annotation \
	 -B ${MAN_DIR}:/manifest \
      	 -B ${SCRATCH_ROOT}:${SCRATCH_ROOT} \
      	 --workdir ${SCRATCH_ROOT} \
      	 ${SIG_IMAGE} ${DOC_SHELL} ${OPT_STR}"

elif [ -e ${TOP_DIR_B} ]; then

    RUN_CMD="docker run -i --rm -w /work \
             -v ${ORD_DIR}:/order \
     	     -v ${MAT_DIR}:/match \
     	     -v ${AQP_DIR}:/aqp \
     	     -v ${OUT_DIR}:/output \
     	     -v ${SEQ_DIR}:/sequence \
     	     -v ${POS_DIR}:/coordinate \
     	     -v ${CAN_DIR}:/canonical \
     	     -v ${GEN_DIR}:/genome \
     	     -v ${IMP_DIR}:/improbe \
     	     -v ${ANN_DIR}:/annotation \
     	     -v ${MAN_DIR}:/manifest \
     	     ${DOC_IMAGE} ${DOC_SHELL} ${OPT_STR}"

else
    echo "Unrecognized Rscript EXE!"
    exit
fi

echo
echo "RUN_CMD = ${RUN_CMD}"
echo

${RUN_CMD}

echo "Done: $0"

# End of file
