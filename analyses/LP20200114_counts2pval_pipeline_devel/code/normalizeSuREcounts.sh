#!/bin/bash

# AUTHOR / DATE
#   Ludo Pagie; January 23, 2020; normalizeSuREcounts.sh

# INTRO / BACKGROUND
#   script to perform the second step of the SuREcounts2SNPcalls pipeline:
#   1. go through single input data (bedpe file for a single cell-line,
#     single bio-rep)
#   2. computed SuRE score, 
#      ie (cDNA-count/iPCR-count)*(total-iPCR-count/total-cDNA-count)
#   3. write new tabular text file with annotation columns (seqname, start,
#      end, SNPs, etc) and computed SuRE-score
#   The input data is a tabular text files. The counts are in 2 or more
#      columns:
#   1. column with iPCR counts (mostly called 'iPCR' or 'counts')
#   2. One or more collumns with cDNA counts, called anything but generally
#      combinations of sample name and bio-rep nr
#   The count columns should be specified by the user as follows:
#   - name of column with iPCR counts: '-i column-name'
#   - name(s) of the column(s) with cDNA count(s): '"-c col-name1 col-name2 ..."'
#   This means that the latter (possibly long) string should be constructed in
#     the Snakemake file.
#   Additional input is the file generated in pipeline-step1, with the total
#     cDNA/iPCR-counts
#
# USAGE / INPUT / ARGUMENTS / OUTPUT
# USAGE:
#   required:
#   -d: input data file (columns: BC chrom start end strand SNP_ABS_POS SNP_REL_POS SNP_ID SNP_SEQ SNP_VAR SNP_PARENT SNP_TYPE SNP_SUBTYPE count I33_B1 start_hg19 end_hg19 SNP_ABS_POS_hg19)
#   -o: output file
#   -i: name of column with iPCR counts
#   -c: name(s) of column(s) with cDNA counts, in single, double quoted string separated by spaces
#   -t: input data file with total iPCR/cDNA counts
#   optional:
#   -l: write to logfile instead of stdout
#   example usage: bash /DATA/usr/ludo/projects/LP190425_flySuRE/analyses/LP20200114_counts2pval_pipeline_devel/getTotalSuREcounts.sh -o /tmp/out -i count -l /tmp/log -c DGRP-324_B1 chr3L.bedpe.gz
# INPUT:
#   tabular text files, compressed, with iPCR/cDNA counts and SuRE fragment
#     annotation. Including columns with names set by user (options -i/-c)
#   tabular text file with total iPCR/cDNA counts for samples specified by '-c'
# OUTPUT:
#   tabular text file, with the same data as the input data file, except
#     columns are added with normalized SuRE-score for each cDNA sample
#     NOTE: the column names with the cDNA-counts are changed by appending
#     "_count". The new columns now are named by the corresponding cDNA sample
#     name
#
# TODO
#   - 

SCRIPTNAME=normalizeSuREcounts.sh

# EXTERNAL SOFTWARE
AWK=gawk

# GLOBAL VARIABLES
DEBUG=0

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
  echo >&2 "usage: ${SCRIPTNAME} -d:o:i:c:t:[l:h]"
  echo >&2 "OPTIONS:"
  echo >&2 "  -d: filename; input file with iPCR/cDNA counts [required]"
  echo >&2 "  -o: filename; output file for normalized SuRE-scores [required]"
  echo >&2 "  -i: string; name of column with iPCR counts [required]"
  echo >&2 "  -c: string; names(s) of column(s) with cDNA counts in double quoted string [required]"
  echo >&2 "  -c: string; name of column with cDNA counts [required];"
  echo >&2 "      multiple names can be specified using the option multiple times"
  echo >&2 "      eg. '-c colname1 -c colname2 ..'"
  echo >&2 "  -t: filename; input file with total iPCR/cDNA counts [required]"
  echo >&2 "  -l: filename; write messages to file [default: stdout]"
  echo >&2 "  -h: flag; print this message"
  echo >&2 ""
}

while getopts "h?d:o:i:c:t:l:h" opt; do
  case $opt in
    d)
      IN=$OPTARG;
      ;;
    o)
      OUT=$OPTARG;
      ;;
    i)
      IPCR=$OPTARG;
      ;;
    c)
      # this option may be used repeatedly, therefor built an array
      CDNA+=("$OPTARG");
      ;;
    t)
      TOTALS=$OPTARG;
      ;;
    l)
      LOG=$OPTARG;
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

# # the remaining CLI arguments are the input data files in (kinda) bedpe format
# # check we have at least 1 remaining arguments
# if [ $# -lt 1 ]; then
#   echo -e "\nerror: too few arguments left after options are parsed (should be at least 1 filename).\nThe remaining args are:"
#   while test $# -gt 0; do
#     echo $1
#     shift
#   done
#   # echo -e "Aborting\n\n"
#   usage
#   error_exit 1 "$LINENO: error; script called with incorrect arguments"
# fi

# no remaining CLI arguments should be left at this point
if [ $# -gt 0 ]; then
  usage
  error_exit 1 "$LINENO: error; too many user arguments"
fi

# check input file(s) exists ($IN, $TOTALS)
if [ -z ${IN+x} ]; then error_exit 1 "$LINENO: option -d not set (name of file with total counts)"; fi
if [ -z ${TOTALS+x} ]; then error_exit 1 "$LINENO: option -t not set (name of file with total counts)"; fi
abort_flag="false"
for f in "${IN}" "${TOTALS}"; do
  if [ ! -f ${f} ]; then
    echo -e "error; input file (${f}) doesn't exist.\n" 
    abort_flag=true
  fi
done
if [ $abort_flag == 'true' ]; then
  error_exit 1 "error; input file(s) do not exist.\n" 
fi
unset abort_flag

# Make name of input file absolute
IN=$( make_path_absolute "${IN}" )
TOTALS=$( make_path_absolute "${TOTALS}" )

# check all required options are set
if [ -z ${OUT+x} ]; then error_exit 1 "$LINENO: option -o not set (directory for output files)"; fi
OUTDIR=$( dirname $OUT ) || error_exit 1 "$LINENO: dirname ($OUT)"
# check required subdirectories existte if not or error_exit
[[ -d $( dirname ${OUT} ) ]] || mkdir -p $( dirname "${OUT}" ) || error_exit $? "$LINENO: can't create directory for output (${OUT})"
# make path to OUT absolute
OUT=$( make_path_absolute "${OUT}" )

# check column names are defined and present in input
if [ -z ${IPCR+x} ]; then error_exit 1 "$LINENO: option -i not set (name of column in input file with iPCR counts)"; fi
if [ -z ${CDNA+x} ]; then error_exit 1 "$LINENO: option -c not set (name(s) of column(s) in input file with cDNA counts)"; fi
# check all required column names (iPCR and cDNA columns) are present in all input files
HEADER=$( zcat "${IN}" | head -1 )
for col in ${IPCR} ${CDNA[@]}; do
#   if ( ! echo $HEADER | grep -wq "${col}" ); then
  if [[ ${HEADER} != *"$col"* ]]; then # test whether sub-string '$c' is contained in HEADER (https://stackoverflow.com/a/229606)
#     abort_flag=true
    error_exit $? "$LINENO: error; column name '$col' not found in input (${IN})"
#     echo "error: column-name \"${col}\" does not exist in input (${IN})"
  fi
done
HEADER=$( cat "${TOTALS}" | head -1 )
for col in ${IPCR} ${CDNA[@]}; do
  if [[ ${HEADER} != *"${col}"* ]]; then
    error_exit $? "$LINENO: error; column name '$col' not found in input (${TOTALS})"
  fi
done


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
echo "reading input from datafile=${IN}"
echo "reading total cDNA/iPCR counts from datafile=${TOTALS}"
echo "file for output=${OUT}"
echo "name of column with iPCR counts=${IPCR}"
echo "name(s) of column(s) with cDNA counts=${CDNA[@]}"
[[ -z ${LOG+x} ]] && echo "LOG=stdout" || echo "LOG=${LOG}"
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

## READ cDNA/iPCR TOTALS FROM FILE ##
#####################################

# read total counts from total-counts file

join_arr() {   local IFS="$1";   shift;   echo "$*"; }
CDNA_JOINED=$( join_arr , "${CDNA[@]}" )

echo "CDNA jioned = "$CDNA_JOINED

read -r -d '' CMD<<-'EOS'
awk -F"\t" -v ipcrname="${IPCR}" -v cdnanames="${CDNA_JOINED}" -v DEBUG="${DEBUG}" '
  BEGIN {
    if(DEBUG!=0) { # Set debug to a value awk understands (so in bash it can be set to either 0 or anything else, like 'true')
      DEBUG=1
    } else {
      DEBUG=0
    }
    # fix array scanning order: use ninteger index and scan in ascending order
    PROCINFO["sorted_in"] = "@ind_num_asc"
    split("", cols)  # array: column-name -> column-index
    split("",scale)  # array: column-name -> scaling-factor (ie total_iPCR/total_cDNA)
    split(cdnanames, cdnaname, ",") # cdnanames are user-supplied string of comma separated cDNA column-names
    # print for debug (to stderr): int-index, column-names
    if (DEBUG) {
      for (c in cdnaname) {
        print c, cdnaname[c] > "/dev/stderr"
      }
    }
  }
  # read header of the total-counts file (1st file argument)
  # iterate over the header, ie the column names,
  #   and set the column name/indices in array 'cols'
  NR==FNR && NR==1 { # header of totals file
    if (DEBUG) {print $0 > "/dev/stderr"}
    for (i=1; i<=NF; i++) {
      if ($i==ipcrname) { # if column-name == ipcr-name
        cols[$i]=i
        continue
      }
      for (c in cdnaname) { # loop over all cdna-names and compare to column name
        if(DEBUG) {
          print "c2: "c, "#"cdnaname[c]"#", "@"$i"@" > "/dev/stderr"
          print "bool: "($i == cdnaname[c]) > "/dev/stderr"
        }
        if ($i==cdnaname[c]) {  # if column-name = cdna-name[c]
          if (DEBUG) { print "catch" > "/dev/stderr" }
          cols[$i]=i
          continue
        }
      }
      if (DEBUG) {print "i = "i, $i > "/dev/stderr"}
    }
    if (DEBUG) {
      for (c in cols) {
        print "cols: "c, cols[c] > "/dev/stderr"
      }
    }
    next
  }
  # data of totals file: 
  # read total counts and compute scaling factors, stored in array 'scale'
  NR==FNR { 
    for (n in cdnaname) {
      scale[cdnaname[n]]=$cols[ipcrname]/$cols[cdnaname[n]]
    }
    next
  }
  # 2nd file, bedpe file with fragments and raw counts
  # The header line is changed by adding columns for the normalized counts, and
  # adapting the column names with the raw counts by appending "_count"
  FNR==1 { 
    # 1st line; read header and set column indices 
    split("",cols) # array: column-name -> column-index
    for (i=1; i<=NF; i++) {
      if ($i==ipcrname) {
        cols[$i]=i
        continue
      }
      for (c in cdnaname) {
        if ($i==cdnaname[c]) {
          cols[$i]=i
          $i=$i"_count"  # change the cdna raw count column name by appending "_count"
          continue
        }
      }
    }
    # create header for writing to output
    header=$1
    for (i=2; i<=NF; i++) {
      header=header"\t"$i
    }
    for (c in cdnaname) {
      header=header"\t"cdnaname[c]
    }
    print header
    next
  }
  { # second file (main data); 
    # read data and print scaled output
    line=$0
    for (c in cdnaname) {
      line = line"\t"($cols[cdnaname[c]]/$cols[ipcrname]) * scale[cdnaname[c]]
    }
    print line
  }
  '
EOS

echo "running following awk command:"
echo "${CMD} $TOTALS $IN > $OUT"
eval "${CMD}" "$TOTALS" <( gzip -cd "$IN" ) | gzip -c > "$OUT"
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
