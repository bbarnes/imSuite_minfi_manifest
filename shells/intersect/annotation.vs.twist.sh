
GEN="GRCh37"

TOP="/Users/bretbarnes/Documents"
EXE="${TOP}/tools/Workhorse-Unstained/scripts/perl/bsc/intersect/intersectAll.pl"
OUT="${TOP}/Projects/Twist/intersectAll/${GEN}/ann_vs_man"
DES="${TOP}/Projects/Twist/manifest/covered_targets_Twist_Methylome_hg19_annotated_collapsed_uniq_final.bed.gz"
ANN="${TOP}/data/workhorse/annotation/Homo_sapiens/${GEN}"

mkdir -p ${OUT}

REG1="${ANN}/EPIC_CORE_CLEAN/ENCODE/DNaseI_FAIRE_ChIP_Synthesis.bed.gz"
REG2="${ANN}/EPIC_CORE_CLEAN/ENCODE/DNaseI_Hypersensitivity_Clusters.bed.gz"
REG3="${ANN}/EPIC_CORE_CLEAN/ENCODE/TFBS_PeakSeq_based_Peaks.bed.gz"
REG4="${ANN}/EPIC_CORE_CLEAN/ENCODE/wgEncodeGencodeBasicV12.gp.bed.gz"
REG5="${ANN}/EPIC_CORE_CLEAN/ENCODE/wgEncodeGencodeCompV12.gp.bed.gz"
REG6="${ANN}/EPIC_CORE_CLEAN/PHANTOM5/PHANTOM5_Enhacners.bed.gz"
REG7="${ANN}/EPIC_CORE_CLEAN/RRBS/RRBS_DMRs.bed.gz"
REG8="${ANN}/EPIC_CORE_CLEAN/UCSC/cpgIslands.bed.gz"
REG9="${ANN}/EPIC_CORE_CLEAN/UCSC/ucscRefGene.bed.gz"
REG10="${ANN}/EPIC_CORE_CLEAN/WGBS/WGBS_DMRs.bed.gz"
REG11="${ANN}/wgEncodeBroadHmm/GRCh36-GRCh37-liftover/ucsc_liftover_main.time-tracker.csv.gz"
REG12="${ANN}/wgEncodeBroadHmm/GRCh36-GRCh37-liftover/wgEncodeBroadHmmGm12878HMM-lift-hg18ToHg19.map.bed.gz"
REG13="${ANN}/wgEncodeBroadHmm/GRCh36-GRCh37-liftover/wgEncodeBroadHmmGm12878HMM-lift-hg18ToHg19.non.bed.gz"
REG14="${ANN}/wgEncodeBroadHmm/GRCh36-GRCh37-liftover/wgEncodeBroadHmmH1hescHMM-lift-hg18ToHg19.map.bed.gz"
REG15="${ANN}/wgEncodeBroadHmm/GRCh36-GRCh37-liftover/wgEncodeBroadHmmH1hescHMM-lift-hg18ToHg19.non.bed.gz"
REG16="${ANN}/wgEncodeBroadHmm/GRCh36-GRCh37-liftover/wgEncodeBroadHmmHepg2HMM-lift-hg18ToHg19.map.bed.gz"
REG17="${ANN}/wgEncodeBroadHmm/GRCh36-GRCh37-liftover/wgEncodeBroadHmmHepg2HMM-lift-hg18ToHg19.non.bed.gz"
REG18="${ANN}/wgEncodeBroadHmm/GRCh36-GRCh37-liftover/wgEncodeBroadHmmHmecHMM-lift-hg18ToHg19.map.bed.gz"
REG19="${ANN}/wgEncodeBroadHmm/GRCh36-GRCh37-liftover/wgEncodeBroadHmmHmecHMM-lift-hg18ToHg19.non.bed.gz"
REG20="${ANN}/wgEncodeBroadHmm/GRCh36-GRCh37-liftover/wgEncodeBroadHmmHsmmHMM-lift-hg18ToHg19.map.bed.gz"
REG21="${ANN}/wgEncodeBroadHmm/GRCh36-GRCh37-liftover/wgEncodeBroadHmmHsmmHMM-lift-hg18ToHg19.non.bed.gz"
REG22="${ANN}/wgEncodeBroadHmm/GRCh36-GRCh37-liftover/wgEncodeBroadHmmHuvecHMM-lift-hg18ToHg19.map.bed.gz"
REG23="${ANN}/wgEncodeBroadHmm/GRCh36-GRCh37-liftover/wgEncodeBroadHmmHuvecHMM-lift-hg18ToHg19.non.bed.gz"
REG24="${ANN}/wgEncodeBroadHmm/GRCh36-GRCh37-liftover/wgEncodeBroadHmmK562HMM-lift-hg18ToHg19.map.bed.gz"
REG25="${ANN}/wgEncodeBroadHmm/GRCh36-GRCh37-liftover/wgEncodeBroadHmmK562HMM-lift-hg18ToHg19.non.bed.gz"
REG26="${ANN}/wgEncodeBroadHmm/GRCh36-GRCh37-liftover/wgEncodeBroadHmmNhekHMM-lift-hg18ToHg19.map.bed.gz"
REG27="${ANN}/wgEncodeBroadHmm/GRCh36-GRCh37-liftover/wgEncodeBroadHmmNhekHMM-lift-hg18ToHg19.non.bed.gz"
REG28="${ANN}/wgEncodeBroadHmm/GRCh36-GRCh37-liftover/wgEncodeBroadHmmNhlfHMM-lift-hg18ToHg19.map.bed.gz"
REG29="${ANN}/wgEncodeBroadHmm/GRCh36-GRCh37-liftover/wgEncodeBroadHmmNhlfHMM-lift-hg18ToHg19.non.bed.gz"

${EXE} -o ${OUT} -r ${DES} \
       -d ${REG1} \
       -d ${REG2} \
       -d ${REG3} \
       -d ${REG4} \
       -d ${REG5} \
       -d ${REG6} \
       -d ${REG7} \
       -d ${REG8} \
       -d ${REG9} \
       -d ${REG10} \
       -d ${REG12} \
       -d ${REG13} \
       -d ${REG14} \
       -d ${REG15} \
       -d ${REG16} \
       -d ${REG17} \
       -d ${REG18} \
       -d ${REG19} \
       -d ${REG20} \
       -d ${REG21} \
       -d ${REG22} \
       -d ${REG23} \
       -d ${REG24} \
       -d ${REG25} \
       -d ${REG26} \
       -d ${REG27} \
       -d ${REG28} \
       -d ${REG29}
       
#        -d ${REG11} \

# End of file
