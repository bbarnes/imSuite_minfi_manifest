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
    echo "EXE=${EXE}"

    INP=/input
    OUT=/output
    WRK=/tmp
    
    # CMD="${RSCRIPT} ${EXE} --Rscript ${RSCRIPT} --datDir=${INP} --outDir=${OUT} $@"
    # CMD="${RSCRIPT} ${EXE} --Rscript ${RSCRIPT} --datDir=${INP} --outDir=${OUT} --manifestPath=/manDirPath $@"
    CMD="${RSCRIPT} ${EXE} --Rscript ${RSCRIPT} --datDir=${INP} --outDir=${OUT} --manDirPath=${WRK} $@"
    
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

${CMD}

echo "done."
    
#    --minNegPval=0.02 \
#    --minOobPval=0.1 \
#    --minNegPerc=98 \
#    --minOobPerc=90 \
#    --minDeltaBeta=0.2 \
#    --percisionBeta=4 \
#    --percisionPval=6 \
#    --verbose=${VER}

# End of file
