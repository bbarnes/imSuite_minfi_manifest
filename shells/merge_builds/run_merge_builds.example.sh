

# EXE="/Users/bretbarnes/Documents/tools/shells/merge_builds/merge_builds.sh"
EXE="/Users/bretbarnes/Documents/tools/Infinium_Methylation_Workhorse/shells/merge_builds/merge_builds.sh"
OUT="/Users/bretbarnes/Documents/scratch"

SAMP="/Users/bretbarnes/Documents/data/human_sample_sheets/NA12878_SampleSheet.csv"
BLD="/Users/bretbarnes/Documents/scratch/RStudio/swifthoof_main/NA12878.v0"
EXP="NA12878.v0"

SAMP="/Users/bretbarnes/Documents/data/CustomContent/AstraZeneca_EM_SamplePrep/AstraZeneca_20201007_Human_SampleSheet.csv.gz"
SAMP="/Users/bretbarnes/Documents/data/CustomContent/AstraZeneca_EM_SamplePrep/AstraZeneca_30042021_Human_SampleSheet.csv.gz"
BLD="/Users/bretbarnes/Documents/scratch/swifthoof_main/EPIC-8x1-EM-Sample-Prep"
EXP="EPIC-8x1-EM-Sample-Prep"

BLD="/Users/bretbarnes/Documents/scratch/swifthoof_main/EPIC-8x1-EM-Sample-Prep.v0"
EXP="EPIC-8x1-EM-Sample-Prep.v0"

WORK="ind"
PLAT="EPIC"
VERS="B4"

${EXE} \
    --buildDir=${BLD} \
    --outDir=${OUT} \
    --runName=${EXP} \
    --workflow=${WORK} \
    --platform=${PLAT} \
    --version=${VERS} \
    --sampleCsv=${SAMP} \
    --clean

# End of file
