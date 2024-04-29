#!/bin/sh

if [ "$#" -lt 1 ]; then
    echo "Usage: $0 program_path [options]"
    exit 1
fi

# echo "Args = $@"

# EXE_NAME=Infinium_Methylation_Workhorse/scripts/R/swifthoof/swifthoof_main.R
EXE_NAME=$1
shift

# val2=$1
# shift

# echo "Args = $@"
# echo "Path = ${path}"
# echo "Val2 = ${val2}"

# EXE_PATH="Infinium_Methylation_Workhorse/scripts/R"

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


    ORD=/order
    MAT=/match
    AQP=/aqp
    OUT=/output

    SEQ=/sequence
    POS=/coordinate
    CAN=/canonical
    
    GEN=/genome
    MAN=/manifest
    IMP=/improbe
    ANN=/annotation

    # cgn_seq_dir
    # cgn_bed_dir
    # canonical_cgn_dir
    # gen_dir
    # man_dir
    # imp_dir
    # ann_dir
    
    CMD="${RSCRIPT} ${EXE} --Rscript ${RSCRIPT} \
    		    --ord_dir=${ORD} \
    		    --mat_dir=${MAT} \
    		    --aqp_dir=${AQP} \
		    --out_dir=${OUT} \
		    --cgn_seq_dir=${SEQ} \
		    --cgn_bed_dir=${POS} \
		    --canonical_cgn_dir=${CAN} \
		    --gen_dir=${GEN} \
		    --man_dir=${MAN} \
		    --imp_dir=${IMP} \
		    --ann_dir=${ANN} $@"
    
elif [ -e ${EXE_B} ]; then
    EXE=${EXE_B}
    
    echo "Local Run..."
    echo "EXE="${EXE}
    
    # INP=$1
    # OUT=$2

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
    
#    --minNegPval=0.02 \
#    --minOobPval=0.1 \
#    --minNegPerc=98 \
#    --minOobPerc=90 \
#    --minDeltaBeta=0.2 \
#    --percisionBeta=4 \
#    --percisionPval=6 \
#    --verbose=${VER}

# End of file
