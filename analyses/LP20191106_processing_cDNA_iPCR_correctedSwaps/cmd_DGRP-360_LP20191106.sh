DRYRUN=""
# DRYRUN="--dryrun "

DATETAG="LP$( date +"%Y%m%d_%H%M" )"
SNAKEFILE=../../code/pipeline/SuRE-snakemake
CONFIG=config-Dm04_DGRP-360_LP20191106.yml
LOG="${CONFIG%.yml}_run-${DATETAG}.log"
NCORES=35
RAM=200
TARGET="bedpe_merged_smpls"
TARGET="merged_ipcr_cdna"
TARGET="reversed_liftover"
TARGET="trim_iPCR"
TARGET="split_bam"
TARGET="trim_iPCR"
TARGET="sorted_cnt_tbls"
#TARGET="bedpe_BC"
# TARGET="bed2coverage_done"

CMD="/usr/bin/time -v nice -19 snakemake ${DRYRUN}-prs ${SNAKEFILE} --use-conda --resources ram=${RAM} --configfile ${CONFIG} --cores $NCORES ${TARGET} &> ${LOG}"
echo "${CMD}"
eval ${CMD}


