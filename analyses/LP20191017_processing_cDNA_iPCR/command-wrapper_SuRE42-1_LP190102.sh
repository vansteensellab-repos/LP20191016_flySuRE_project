DRYRUN=""
# DRYRUN="--dryrun "

DATETAG="LP$( date +"%y%m%d_%H%M" )"
SNAKEFILE=pipeline/SuRE-snakemake
CONFIG=config-SuRE42_1_LP181224.yml
LOG="${CONFIG%.yml}_run-${DATETAG}.log"
NCORES=30
TARGET="bedpe_merged_smpls"
TARGET="merged_ipcr_cdna"
TARGET="reversed_liftover"
TARGET="sorted_cnt_tbls"
# TARGET="bed2coverage_done"

CMD="/usr/bin/time -v nice -19 snakemake "${DRYRUN}"-prs "${SNAKEFILE}" --use-conda --resources ram=150 --configfile "${CONFIG}" --cores $NCORES ${TARGET} &> "${LOG}""
echo "${CMD}"
eval ${CMD}


