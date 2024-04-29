#!/bin/sh

RSCRIPT=/usr/bin/Rscript
EXE_NAME=Infinium_Methylation_Workhorse/scripts/R/swifthoof/swifthoof_main.R
EXE=/repo/${EXE_NAME}
INP=/input
OUT=/output
CMD="${RSCRIPT} ${EXE} --Rscript ${RSCRIPT} --idatsDir=${INP} --outDir=${OUT} $@"

echo "Docker Run..."
    
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
