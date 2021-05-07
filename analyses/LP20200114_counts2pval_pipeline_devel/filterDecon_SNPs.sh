#!/bin/bash

# AUTHOR / DATE
#   Ludo Pagie; January 29, 2020; filterDecon_SNPs.sh

# INTRO / BACKGROUND
#   script to perform the third step of the SuREcounts2SNPcalls pipeline:
#   1. go through single input data (bedpe file for a single cell-line,
#     single bio-rep, with normalized SuRE-scores)
#   2. Filter:
#      - Discard data without SNPs
#      - Discard secondary SNPs, if SNP is read twice in single fragment (in forw/rev read direction)
#      - Remove NAs
#      - Discard ambiguous SNP data
#   3. Deconvolute; for each fragment containing multiple SNPs; split the
#      SNP-data and associate each SNP with the same SuRE annotation data.
#   The input data is a tabular text files. The normalized SuRE-scores are in 1 or more
#      columns:
#      - One or more columns with normalized SuRE-scores, called anything, but generally
#        combinations of sample name and bio-rep nr
#   The SuRE-score columns should be specified by the user as follows:
#      - name(s) of the column(s) with cDNA count(s): '"-c col-name1 col-name2 ..."'
#   This means that the latter (possibly long) string should be constructed in
#      the Snakemake file.
#
# USAGE / INPUT / ARGUMENTS / OUTPUT
# USAGE:
#   required:
#   -d: input data file (columns: BC chrom start end strand SNP_ABS_POS SNP_REL_POS SNP_ID SNP_SEQ SNP_VAR SNP_PARENT SNP_TYPE SNP_SUBTYPE count I33_B1_count start_hg19 end_hg19 SNP_ABS_POS_hg19 I33_B1)
#   -o: output file
#   -c: name(s) of column(s) with SuRE-scores, in single, double quoted string separated by spaces
#   optional:
#   -l: write to logfile instead of stdout
#   example usage: bash /DATA/usr/ludo/projects/LP190425_flySuRE/analyses/LP20200114_counts2pval_pipeline_devel/filterDecon_SNPs.sh -o /tmp/out -l /tmp/log -c DGRP-324_B1 chr3L.bedpe.gz
# INPUT:
#   tabular text files, compressed, with normalized SuRE-scores and SuRE fragment
#     annotation. Including columns with names set by user (options -c)
# OUTPUT:
#   tabular text file, with the same data as the input data file, with some
#     data discarded and other data split in multiple lines
#
# TODO
#   - 

SCRIPTNAME=filterDecon_SNPs.sh

# EXTERNAL SOFTWARE
AWK=gawk

# GLOBAL VARIABLES
# 

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
  echo >&2 "usage: ${SCRIPTNAME} -d:o:c:[l:h]"
  echo >&2 "OPTIONS:"
  echo >&2 "  -d: filename; input file with normalized SuRE-scores [required]"
  echo >&2 "  -o: filename; output file for filtered and deconvoluted data [required]"
  echo >&2 "  -c: string; names(s) of column(s) with SuRE-scores in double quoted string [required]"
  echo >&2 "  -l: filename; write messages to file [default: stdout]"
  echo >&2 "  -h: flag; print this message"
  echo >&2 ""
}

while getopts "h?d:o:c:l:h" opt; do
  case $opt in
    d)
      IN=$OPTARG;
      ;;
    o)
      OUT=$OPTARG;
      ;;
    c)
      SCORE_COLS=$OPTARG;
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
if [ ! -f ${IN} ]; then
  error_exit 1 "error; input file (${IN}) do not exist.\n" 
fi

# Make name of input file absolute
IN=$( make_path_absolute "${IN}" )

# check all required options are set
if [ -z ${OUT+x} ]; then error_exit 1 "$LINENO: option -o not set (directory for output files)"; fi
OUTDIR=$( dirname $OUT ) || error_exit 1 "$LINENO: dirname ($OUT)"
# check required subdirectories exist, create if not or error_exit
[[ -d $( dirname ${OUT} ) ]] || \
  mkdir -p $( dirname "${OUT}" ) || \
  error_exit $? "$LINENO: can't create directory for output (${OUT})"
# make path to OUT absolute
OUT=$( make_path_absolute "${OUT}" )

# check column names are defined and present in input
if [ -z ${SCORE_COLS+x} ]; then error_exit 1 "$LINENO: option -i not set (name of column in input file with SuRE-scores)"; fi
abort_flag="false"
header=$( zcat "${IN}" | head -1 )
for col in ${SCORE_COLS}; do
  if ( ! echo $header | grep -wq "${col}" ); then
    abort_flag=true
    echo "error: column-name \"${col}\" does not exist in input (${IN})"
  fi
done
if [ $abort_flag == 'true' ]; then
  error_exit 1 "$LINENO; error: some column-name(s) does not exist in input\n" 
fi
unset abort_flag


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
echo "file for output=${OUT}"
echo "name(s) of column(s) with nrmalized SuRE-scores=${SCORE_COLS}"
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

#   -d: input data file (columns: BC chrom start end strand SNP_ABS_POS SNP_REL_POS SNP_ID SNP_SEQ SNP_VAR SNP_PARENT SNP_TYPE SNP_SUBTYPE count I33_B1_count start_hg19 end_hg19 SNP_ABS_POS_hg19 I33_B1)

#   2. Filter:
#      - Discard data without SNPs
#      - Discard secondary SNPs, if SNP is read twice in single fragment (in forw/rev read direction)
#      - Remove NAs
#      - Discard ambiguous SNP data
#   3. Deconvolute; for each fragment containing multiple SNPs; split the
#      SNP-data and associate each SNP with the same SuRE annotation data.


read -r -d '' CMD<<-'EOS'
awk -F"\t" -v scorecols="${SCORE_COLS}" '
  BEGIN {
    split("", cols)
    split("",totals)
    split("",scale)
    totals[ipcrname] = 0
    split(scorecols, scorecol)
    for (c in scorecol) {
      totals[scorecol[c]] = 0
    }
  }
  NR==FNR && NR==1 { # header of totals file
    for (i=1; i<=NF; i++) {
      if ($i==ipcrname) {
        cols[$i]=i
        continue
      }
      for (c in cdnaname) {
        if ($i==cdnaname[c]) {
          cols[$i]=i
          continue
        }
      }
    }
    next
  }
  NR==FNR { # data of totals file
    for (n in cols) {
      totals[n]=$cols[n]
    }
    for (n in cdnaname) {
      scale[cdnaname[n]]=totals[ipcrname]/totals[cdnaname[n]]
    }
    next
  }
  FNR==1 { # 2nd file, 1st line; read header and set column indices
    split("",cols) 
    for (i=1; i<=NF; i++) {
      if ($i==ipcrname) {
        cols[$i]=i
        continue
      }
      for (c in cdnaname) {
        if ($i==cdnaname[c]) {
          cols[$i]=i
          $i=$i"_count"
          continue
        }
      }
    }
    header=$0
    print "adapting header"
    for (c in cdnaname) {
      header=header"\t"cdnaname[c]
    }
    print header
    print "header printed"
    next
  }
  { # second file (main data), past header; read data and print scaled output
    line=$0
    for (c in cdnaname) {
      line = line"\t"($cols[cdnaname[c]]/$cols[ipcrname]) * scale[cdnaname[c]]
    }
    print line
  }
  '
EOS

# read -r -d '' CMD<<-'EOS'
# awk -F"\t" -v ipcrcol=${IPCR} -v cdnacols="${CDNA}" -v ipcrtotal="" -v cdnatotals="" '
#   BEGIN{
#     counts["ipcr"]=0
#     split(cdna, cdnanames)
#     for (c in cdnanames) {
#       counts[c]=0
#     }
#     split("",cdnacols," ")
#   }
#   NR==1 { # read header line first time; set correct column indices
#     for (i=1; i<=NF; i++) {
#       if ($i==ipcr) {
#         ipcr=i
#         continue
#       }
#       for (c=1; c<=length(cdnanames); c++) {
#         if ($i==cdnanames[c]) {
#           cdnacols[cdnanames[c]]=i
#           continue
#         }
#       }
#     }
#     next
#   }
#   FNR==1 { # read header lines in subsequent files, no need to set column indices again
#     next
#   }
#   { # read data, add counts to counters
#     counts["ipcr"]+=$ipcr
#     for (c in cdnanames) {
#       counts[cdnanames[c]]+=$cdnacols[cdnanames[c]]
#     }
#   }
#   END {
#     # print header
#     header = "FRAGM_COUNT\tiPCR"
#     values = (NR-ARGC+1)"\t"counts["ipcr"] # FRAGM_COUNT is number of records read, minus number of input files to account for headers
#     PROCINFO["sorted_in"] = "@ind_str_asc"
#     for (c in cdnanames) {
#       header = header"\t"cdnanames[c]
#       values = values"\t"counts[cdnanames[c]]
#     }
#     print header
#     print values
#   }'
# 
# EOS

# eval "${CMD}" "${INPUT}" > "${OUT}"
echo "${CMD} $TOTALS $IN > $OUT"
eval "${CMD}" "$TOTALS" <( gzip -cd "$IN" ) > "$OUT"
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
