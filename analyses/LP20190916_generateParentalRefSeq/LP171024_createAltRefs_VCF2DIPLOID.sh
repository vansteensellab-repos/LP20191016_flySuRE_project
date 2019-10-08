#!/bin/bash

# set some paths to data files and runfiles
JAR="/home/ludo/projects/LP140430_SureSeq_JvArensbergen/analyses/LP171023_createAltRefs/VCF2DIPLOID/vcf2diploid_v0.2.6a/vcf2diploid.jar"
CHR="/home/ludo/projects/LP140430_SureSeq_JvArensbergen/analyses/LP171023_createAltRefs/fasta/hs37d5.fa"
VCF="/home/ludo/projects/LP140430_SureSeq_JvArensbergen/analyses/LP171023_createAltRefs/VCF/1000Genomes_1-22-X_HG02601-NA18983-HG01241-HG03464_strucVarFltr_LP171107.vcf"
DATE=$(date +"%y%m%d")
OUTDIR="${PWD}/fasta/altRef_${DATE}"

# define samples from 1000genome project to process
declare -a SAMPLES=("HG02601" "NA18983" "HG01241" "HG03464")

# loop pver samples and create all ref sequences
for sample in "${SAMPLES[@]}"; do
  echo "$sample"
  echo "${OUTDIR}/${sample}"
  CMD=( java -jar "${JAR}" -id "${sample}" -pass -chr "${CHR}" -vcf "${VCF}" -outDir "${OUTDIR}/${sample}" )
  echo "CMD = ${CMD[@]}"
  mkdir -p "${OUTDIR}/${sample}"
  eval "${CMD[@]}"
done


