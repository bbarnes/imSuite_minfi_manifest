#!/bin/bash
#! -N sesame_build
#! -pe threaded 28
#! -q whole
#! -V

# image_id="Infinium_Methylation_Workhorse"
# VER="v.4.3"

image_id="Infinium_Methylation_Workhorse_Centos"
# VER="v.1.13"
VER="v.1.11"
VER="v.1.12"
VER="v.1.21"
VER="v.1.24"
VER="v.1.25"
VER="v.1.27"
VER="v.1.29"
VER="v.1.30"
VER="v.1.31"
VER="v.1.32"
VER="v.1.33"
VER="v.1.34"
VER="v.1.35"
VER="v.1.36"
VER="v.1.37"
VER="v.1.38"
VER="v.1.39"
VER="v.1.40"
VER="v.1.41"
VER="v.1.42"
VER="v.1.43"
VER="v.1.44"

VER="v.1.49"

outdir="/illumina/scratch/darkmatter/docker/software"
mkdir -p ${outdir}

docker_source="docker://bbarnesimdocker/im_workhorse:${image_id}.${VER}"

# output_sif="/illumina-services/scratch/ils_methylation/sesame_testing/software/${image_id}.${VER}.sif"
output_sif="${outdir}/${image_id}.${VER}.sif"

echo "docker_source=${docker_source}"
echo "output_sif=${output_sif}"

cmd="singularity build ${output_sif} ${docker_source}"
echo ${cmd}
eval ${cmd}

echo "done singularity."

# End of file
