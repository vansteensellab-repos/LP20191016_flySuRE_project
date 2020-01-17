#!/bin/bash

# AUTHOR / DATE
#   Ludo Pagie; January 14, 2020; getTotalSuREcounts.sh

# INTRO / BACKGROUND
#   script to perform the first steps of the SuREcounts2SNPcalls pipeline:
#   1. go through all input data (all bedpe files for a single cell-line,
#     single bio-rep, whole genome)
#   2. determine total counts for iPCR and cDNA
#   3. write total counts to some intermediate file
#   The input data is (a collection of) tabular text files. The counts are in 2
#   or more columns:
#   1. column with iPCR counts (called 'iPCR' or 'counts')
#   2. One or more collumns with cDNA counts, called anything but generally
#      combinations of sample name and bio-rep nr
#   The count columns should be specified by the user as follows:
#   - name of column with iPCR counts: '-i column-name'
#   - name(s) of the column(s) with cDNA count(s): '"-c col-name1 col-name2 ..."'
#   The means that the latter (possibly long) string should be constructed in
#     the Snakemake file.
#
# USAGE / INPUT / ARGUMENTS / OUTPUT
# USAGE:
#   required:
#   - trailing arguments: (multiple) input data files (columns: BC chrom start end strand SNP_ABS_POS SNP_REL_POS SNP_ID SNP_SEQ SNP_VAR SNP_PARENT SNP_TYPE SNP_SUBTYPE count I33_B1 start_hg19 end_hg19 SNP_ABS_POS_hg19)
#   -o: output file for total counts
#   -i: name of coilumn with iPCR counts
#   -c: name(s) of column(s) with cDNA counts, in single, double quoted string separated by spaces
#   optional:
#   -l: write to logfile instead of stdout
#   -n: number of cores used in parallel processes (10)
# INPUT:
#   trimmed fastq files
# OUTPUT:
#   bam files with aligned reads
#
# VERSIONS:
#   -170206: initial version, VERSION set to 0.0.1
#
# TODO
#   - parameterize filter criteria (min MAPQ score, BC length, etc)

SCRIPTNAME=getTotalSuREcounts.sh

# EXTERNAL SOFTWARE
AWK=gawk

# GLOBAL VARIABLES
# LOG="false"

# ERROR_EXIT FUNCTION
error_exit()
{
#	----------------------------------------------------------------
#	Function for exit due to fatal program error
#		Accepts 2 arguments:
#     1. exit code of last run program
#			2. string containing descriptive error message
#	----------------------------------------------------------------

	echo -e "${SCRIPTNAME}: ${2:-"Unknown Error"}" 1>&2 # print $2 or "Unkno..." if $2 not exists
	exit $1
}

# # Example call of the error_exit function.  Note the inclusion
# # of the LINENO environment variable.  It contains the current
# # line number.
# echo "Example of error with line number and message"
# error_exit "$LINENO: An error has occurred."
# 

# FUNCTION TO MAKE FILE PATHS ABSOLUTE
make_path_absolute() {
# ----------------------------------------------------------------
# Function to make a file-/dir-path absolute
# accepts 1 argument: path
# ----------------------------------------------------------------

  D=`dirname "$1"` || error_exit $? "$LINENO: error dirname $1"
  B=`basename "$1"` || error_exit $? "$LINENO: error; basename $1"
  DD="`cd $D 2>/dev/null && pwd || echo $D`" || error_exit $? "$LINENO: error; making filename absolute ($1)"
  echo "$DD/$B"
}

# PARSE OPTIONS
OPTIND=1         # Reset in case getopts has been used previously in the shell.
USAGE=
usage() {
  echo >&2 "usage: ${SCRIPTNAME} -o:i:c:[l:h] input-bedpe-data-file(s)"
  echo >&2 "OPTIONS:"
  echo >&2 "  -o: filename; output file for total counts [required]"
  echo >&2 "  -i: string; name of column with iPCR counts [required]"
  echo >&2 "  -c: string; names(s) of column(s) with cDNA counts in double quoted string [required]"
  echo >&2 "  -l: filename; write messages to file [default: stdout]"
  echo >&2 "  -h: flag; print this message"
  echo >&2 ""
}

while getopts "h?o:l:n:i:c:" opt; do
  case $opt in
    l)
      LOG=$OPTARG;
      ;;
    o)
      OUT=$OPTARG;
      ;;
    i)
      IPCR=$OPTARG;
      ;;
    c)
      CDNA=$OPTARG;
      ;;
    h)
      usage;
      ;;
    \?)
      # echo "option not recognized: "$opt
      usage
      error_exit 1 "$LINENO: unrecognized option: $opt"
      ;;
  esac
done
shift $(( OPTIND - 1 ))

# the remaining CLI arguments are the input data files in (kinda) bedpe format
# check we have at least 1 remaining arguments
if [ $# -lt 1 ]; then
  echo -e "\nerror: too few arguments left after options are parsed (should be at least 1 filename).\nThe remaining args are:"
  while test $# -gt 0; do
    echo $1
    shift
  done
  # echo -e "Aborting\n\n"
  usage
  error_exit 1 "$LINENO: error; script called with incorrect arguments"
fi

# retrieve input bedpe files from command line
declare -a BEDPE=( "$@" );

# check input file(s) exists
abort_flag="false"
for f in ${BEDPE[@]}; do
  if [ ! -f ${f} ]; then
    echo -e "error; input file (${f}) doesn't exist.\n" 
    abort_flag=true
  fi
done
if [ $abort_flag == 'true' ]; then
  error_exit 1 "error; input file(s) do not exist.\n" 
fi
unset abort_flag

# Make name of bedpe file absolute
for (( i=0; i<${#BEDPE[@]}; i++ )); do
  D=`dirname "${BEDPE[$i]}"` || error_exit $? "$LINENO: error dirname ${BEDPE[$i]}"
  B=`basename "${BEDPE[$i]}"` || error_exit $? "$LINENO: error; basename ${BEDPE[$i]}"
  DD="`cd $D 2>/dev/null && pwd || echo $D`" || error_exit $? "$LINENO: error; making filename absolute (${BEDPE[$i]})"
  BEDPE[$i]=$( make_path_absolute "${BEDPE[$i]}" )
done

# check all required options are set
if [ -z ${OUT+x} ]; then error_exit 1 "$LINENO: option -o not set (directory for output files)"; fi
OUTDIR=$( dirname $OUT ) || error_exit 1 "$LINENO: dirname ($OUT)"
# check required subdirectories existte if not or error_exit
[[ -d $( dirname ${OUT} ) ]] || mkdir -p $( dirname "${OUT}" ) || error_exit $? "$LINENO: can't create directory for output (${OUT})"
# make path to OUT absolute
OUT=$( make_path_absolute "${OUT}" )
if [ -z ${IPCR+x} ]; then error_exit 1 "$LINENO: option -i not set (name of column in input file with iPCR counts)"; fi
if [ -z ${CDNA+x} ]; then error_exit 1 "$LINENO: option -c not set (name(s) of column(s) in input file with cDNA counts)"; fi

######################################
# write stdout to stdout or a log file
######################################
[[ -z ${LOG+x} ]] || exec 1>>${LOG}

# print values of variables and CLI args for log
# print header for log
######################
LINE="running "${SCRIPTNAME}
SEPARATOR=$(head -c ${#LINE} </dev/zero | tr '\0' '=')
echo $SEPARATOR; 
echo $SEPARATOR
echo $LINE; 
echo $SEPARATOR
echo $SEPARATOR
echo "script context"
echo "=============="
starttime=$(date +%c)
echo "starting date/time = "${starttime}
echo "User set variables:"
echo "==================="
echo "file for output=${OUT}"
echo "name of column with iPCR counts=${IPCR}"
echo "name(s) of column(s) with cDNA counts=${CDNA}"
echo "file for output=${OUT}"
[[ -z ${LOG+x} ]] && echo "LOG=stdout" || echo "LOG=${LOG}"
echo "NCORES=${NCORES}"
echo ""
echo "bedpe files for input:"
echo "================================="
for f in ${BEDPE[@]}; do echo $f; done
echo ""
# print some software version info
echo "Used software:"
echo "=============="
echo "unix/host"; uname -a; echo "---------------";
echo "bash:"; bash --version 2>&1 head -3; echo "---------------";
echo "awk:"; echo "executable used: ${AWK}"; ${AWK} --version; echo "---------------";
echo "=============="
echo ""

echo -e "finished prepping for processing"
echo -e "================================\n\n"

echo "====================================="
echo "====================================="
echo "MAIN: starting to process bedpe files" 
echo "====================================="
echo "====================================="
echo ""

#################################
#######  MAIN  ##################
#################################

INPUT="<( gzip -dc ${BEDPE[0]} )"
for ((i = 1; i < ${#BEDPE[@]}; i++)); do
  INPUT+=" <( gzip -dc ${BEDPE[$i]} )"
done
echo "INPUT = ${INPUT}"

read -r -d '' CMD<<-'EOS'
awk -F"\t" -v ipcr=${IPCR} -v cdna="${CDNA}" '
  BEGIN{
    counts["ipcr"]=0
    split(cdna, cdnanames)
    for (c in cdnanames) {
      counts[c]=0
    }
    split("",cdnacols," ")
  }
  NR==1 { # read header line first time; set correct column indices
    for (i=1; i<=NF; i++) {
      if ($i==ipcr) {
        ipcr=i
        continue
      }
      for (c=1; c<=length(cdnanames); c++) {
        if ($i==cdnanames[c]) {
          cdnacols[cdnanames[c]]=i
          continue
        }
      }
    }
    next
  }
  FNR==1 { # read header lines in subsequent files, no need to set column indices again
    next
  }
  { # read data, add counts to counters
    counts["ipcr"]+=$ipcr
    for (c in cdnanames) {
      counts[cdnanames[c]]+=$cdnacols[cdnanames[c]]
    }
  }
  END {
    # print header
    header = "FRAGM_COUNT\tiPCR"
    values = (NR-ARGC+1)"\t"counts["ipcr"] # FRAGM_COUNT is number of records read, minus number of input files to account for headers
    PROCINFO["sorted_in"] = "@ind_str_asc"
    for (c in cdnanames) {
      header = header"\t"cdnanames[c]
      values = values"\t"counts[cdnanames[c]]
    }
    print header
    print values
  }'

EOS

eval "${CMD}" "${INPUT}" > "${OUT}"
ec=$?
[ "${ec}" -eq 0 ] || error_exit "${ec}" "line $LINENO: error in main command"

##############################
########## DONE ##############
##############################
LINE="finished "${SCRIPTNAME}
SEPARATOR=$(head -c ${#LINE} </dev/zero | tr '\0' '=')
echo $SEPARATOR; 
echo $SEPARATOR; 
echo $LINE; 
echo $SEPARATOR
echo $SEPARATOR; 
endtime=$(date +%c)
echo "end date/time = "${endtime}
echo $SEPARATOR; 
echo $SEPARATOR; 
