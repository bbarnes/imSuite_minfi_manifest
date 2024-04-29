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

echo "RSC="${RSCRIPT}

if [ -e ${EXE_A} ]; then
    EXE=${EXE_A}

    echo "Docker Run..."
    echo "EXE="${EXE}

    INP=/input
    OUT=/output
    CMD="${RSCRIPT} ${EXE} --Rscript ${RSCRIPT} --idatsDir=${INP} --outDir=${OUT} $@"
    
elif [ -e ${EXE_B} ]; then
    EXE=${EXE_B}
    
    echo "Local Run..."
    echo "EXE="${EXE}
    
    INP=$1
    OUT=$2

    if [ "$#" -lt 2 ]; then
	echo "Usage: $0 idatsDir outDir [options]"
	exit 1
    fi
    
    CMD="${RSCRIPT} ${EXE} --Rscript ${RSCRIPT} $@"
    
else
    echo "Unrecognized EXE directory!"
    exit
fi

echo "INP="${INP}
echo "OUT="${OUT}
echo ""

echo "CMD="${CMD}
echo ""

# echo ${RSCRIPT} ${EXE} --Rscript ${RSCRIPT} --idatsDir=${INP} --outDir=${OUT} $@
${CMD}

echo "done."

exit

${RSCRIPT} ${EXE} \
    --Rscript ${RSCRIPT} \
    -i ${INP} -o ${OUT} \
    --workflows="i,ind" \
    --writeCalls \
    --writeSsheet


    
#    --minNegPval=0.02 \
#    --minOobPval=0.1 \
#    --minNegPerc=98 \
#    --minOobPerc=90 \
#    --minDeltaBeta=0.2 \
#    --percisionBeta=4 \
#    --percisionPval=6 \
#    --verbose=${VER}

# End of file
