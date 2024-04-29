#!/bin/sh

EXE_NAME=Infinium_Methylation_Workhorse/scripts/R/swifthoof/swifthoof_main.R
EXE_A=/repo/${EXE_NAME}
EXE_B=/Users/bretbarnes/Documents/tools/${EXE_NAME}

RSCRIPT_A=/usr/local/bin/Rscript
RSCRIPT_B=/usr/bin/Rscript

if [ -e ${RSCRIPT_A} ]; then
    RSCRIPT=${RSCRIPT_A}
elif [ -e ${RSCRIPT_B} ]; then
    RSCRIPT=${RSCRIPT_B}
else
    echo "Unrecognized Rscript EXE!"
    exit
fi

if [ -e ${EXE_A} ]; then
    EXE=${EXE_A}

    INP=/input
    OUT=/output
    CMD="${RSCRIPT} ${EXE} --Rscript ${RSCRIPT} --outDir=${OUT} $@"

    echo "Docker Run..."
    
elif [ -e ${EXE_B} ]; then
    EXE=${EXE_B}

    # INP=$1
    # OUT=$2

    CMD="${RSCRIPT} ${EXE} --Rscript ${RSCRIPT} $@"
    echo "Local Run..."
    
else
    echo "Unrecognized EXE directory!"
    exit
fi

echo "RSC="${RSCRIPT}
echo "EXE="${EXE}
echo "INP="${INP}
echo "OUT="${OUT}
echo "CMD="${CMD}
echo ""

# echo ${RSCRIPT} ${EXE} --Rscript ${RSCRIPT} --idatsDir=${INP} --outDir=${OUT} $@
${CMD}

echo "done."

# End of file
