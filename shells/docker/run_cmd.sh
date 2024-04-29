#!/bin/sh

if [ "$#" -lt 1 ]; then
    echo "Usage: $0 program_path [options]"
    exit 1
fi

# echo "Args = $@"

# EXE_NAME=Infinium_Methylation_Workhorse/scripts/R/swifthoof/swifthoof_main.R
EXE_NAME=$1
shift

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
    
    CMD="${RSCRIPT} ${EXE} --Rscript ${RSCRIPT} $@"
    
elif [ -e ${EXE_B} ]; then
    EXE=${EXE_B}
    
    echo "Local Run..."
    echo "EXE="${EXE}

    if [ "$#" -lt 2 ]; then
	echo "Usage: $0 idatsDir outDir [options]"
	exit 1
    fi
    
    CMD="${RSCRIPT} ${EXE} --Rscript ${RSCRIPT} $@"
    
else
    echo "Unrecognized EXE directory!"
    exit
fi

echo "OUT="${OUT}
echo ""

echo "CMD="${CMD}
echo ""

${CMD}

echo "done."

# End of file
