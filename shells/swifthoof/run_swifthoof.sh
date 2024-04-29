#!/bin/bash


# EXP="Customer-Facing"

# Running below...

# Done below...
#
# EXP="Excalibur-New-1609202"
# EXP="Excalibur-Old-1609202"
#
# EXP="BETA-8x1-EPIC-Bad"
#
# EXP="BETA-8x1-EPIC-Core"
# EXP="DELTA-8x1-EPIC-Core"
#
# EXP="COVIC-Set1-15052020"
# EXP="COVIC-Set2-31052020"
# EXP="COVIC-Set3-05062020"
# EXP="COVIC-Set4-09062020"
# EXP="COVIC-Set5-10062020"
# EXP="COVIC-Set7-06082020"
#
# EXP="COVIC-Set8-26182020"
#
# EXP="CNTL-Samples_VendB_10092020"
# EXP="CNTL-Samples_VendA_10092020"


SRC="/Users/bretbarnes/Documents/tools/Infinium_Methylation_Workhorse"
TOP="/Users/bretbarnes/Documents"

EXE="${SRC}/scripts/R/swifthoof/swifthoof_main.R"
# AUTO=${SRC}/dat/ref/AutoSampleDetection_EPIC-B4_8x1_pneg98_Median_beta_noPval_BETA-Zymo_Mean-COVIC-280-NP-ind_negs-0.02.csv.gz

DAT="${TOP}/data/idats/idats_${EXP}"

OUT="${TOP}/scratch"

# WORK="nd,ind"
WORK="ind"
PVAL_KEY="pOOBAH,PnegEcdf"
PVAL_MIN="0.1,0.02"
PVAL_CUT="90,98"

IS_SINGLE="--single"
IS_SINGLE=""
VERBOSE=4

RUN="Rscript ${EXE} \
	--Rscript=Rscript \
	--outDir=${OUT} \
	--idatsDir=${DAT} \
	--workflow=${WORK} \
	--runName=${EXP} \
	--pval=${PVAL_KEY} \
	--minPval=${PVAL_MIN} \
	--minPerc=${PVAL_CUT} \

	--dpi=120 \
	--manDirName=core \
	--manifestPath=auto \
	--minDeltaBeta=0.2 \
	--percision_beta=4 \
	--percision_pval=6 \
	--percision_sigs=1 \
	--plotFormat=png \
	--plotMax=10000 \
	--plotSub=5000 \

	--load_idat \
	--load_sset \
	--save_sset \

	--write_call \
	--write_sigs \
	--write_ssum \
	--write_csum \

        --verbose=${VERBOSE} 
        --auto_detect \
	--parallel \
        ${IS_SINGLE}"

echo "RUN=${RUN}"
eval $RUN

echo "done swifthoof..."
exit

# 	--save_idat \

echo "done swifthoof..."
exit

# 	--auto_detect \
#	--auto_sam_csv=${AUTO} \
#	--single


# End of file
