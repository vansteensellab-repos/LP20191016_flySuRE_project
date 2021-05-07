# AUTHOR / DATE
#   Ludo Pagie; January 30, 2020; filterDecon_SNPs.py

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
#   -input: input data file (columns: BC chrom start end strand SNP_ABS_POS SNP_REL_POS SNP_ID SNP_SEQ SNP_VAR SNP_PARENT SNP_TYPE SNP_SUBTYPE count I33_B1_count start_hg19 end_hg19 SNP_ABS_POS_hg19 I33_B1)
#   -output: output file
#   -columns: name(s) of column(s) with SuRE-scores, in single, double quoted string separated by spaces
#   optional:
#   -log: write to logfile instead of stdout
#   example usage: python /DATA/usr/ludo/projects/LP190425_flySuRE/analyses/LP20200114_counts2pval_pipeline_devel/filterDecon_SNPs.py -out /tmp/out -log /tmp/log -columns DGRP-324_B1 -input chr3L.bedpe.gz
# INPUT:
#   tabular text files, compressed, with normalized SuRE-scores and SuRE fragment
#     annotation. Including columns with names set by user (options -columns)
# OUTPUT:
#   tabular text file, with the same data as the input data file, with some
#     data discarded and other data split in multiple lines
#
# TODO
#   - 

SCRIPTNAME="filterDecon_SNPs.py"

import sys
import argparse
# import os.path
# import glob
import gzip
import logging

REF_POS_COL = dict( (('start',"start_hg19"), ('end',"end_hg19"), ('snp',"SNP_ABS_POS_hg19")))
SNP_COL_NAMES = ('SNP_ABS_POS','SNP_REL_POS','SNP_ID','SNP_SEQ','SNP_VAR','SNP_PARENT','SNP_TYPE','SNP_SUBTYPE',REF_POS_COL['snp'])
SNP_COL_NAMES = ('SNP_ABS_POS','SNP_ID','SNP_SEQ','SNP_VAR','SNP_PARENT','SNP_TYPE','SNP_SUBTYPE',REF_POS_COL['snp'])

def parse_options():
# parse user options:
    parser = argparse.ArgumentParser(
        description="Filter and deconvolute data with regard to SNPs")
    parser.add_argument('-i', '--input', required=True,
                        help=('input data file (compressed)'))
    parser.add_argument('-o', '--output', required=True, 
                        help=('output bedpe-like file'))
#    parser.add_argument('-c', '--columns', required=True, nargs='+',
#                        help=('name(s) of column(s) containing SuRE-score(s)'))
    parser.add_argument('-l', '--log', required=True,
                        help=('log file'))

    options = parser.parse_args()
    return options

def init_logger(log_fname):
    # setup a log 'device'
    logger = logging.getLogger(SCRIPTNAME)
    logger.setLevel(logging.DEBUG)
    # handler
    hndlr = logging.FileHandler(log_fname, mode='w')
    hndlr.setFormatter(logging.Formatter('%(name)s:%(levelname)s:%(asctime)s:%(lineno)d - %(message)s', 
        datefmt='%d-%b-%y %H:%M:%S'))
    logger.addHandler(hndlr)
    return logger
    # use as: logger.info('my message'), or: logger.error('error msg')

def discardFragment_no_SNP(vals, colindex):
    # filter data from input ('vals') if it does not contain SNP data
    # return filtered data ('vals' or an empty list [])
    if len(vals)==0: # check whether input is empty already
        return vals
    if vals[colindex['SNP_ID']] == "":
        return []
    return vals

def discardFragment_NA(vals,colindex):
    # filter data from input ('vals') if data contains empty values ("") in position columns
    # return filtered data ('vals' or an empty list [])
    if len(vals)==0: # check whether input is empty already
        return vals
    if vals[colindex[REF_POS_COL['start']]]=="" or vals[0][colindex[REF_POS_COL['end']]]=="" or vals[0][colindex[REF_POS_COL['snp']]]=="":
        return []
    return vals

def deconvolute_SNPs(vals,colindex):
    # split data if it contains multiple SNP annotations
    # return split data ('newvals' or an empty list [])
    if len(vals)==0: # check whether input is empty already
        return vals

    # if only single SNP is annotated return 'vals' directly
    if ',' not in vals[colindex['SNP_ABS_POS']]:
        return [vals]

    # split all columns ('SNP_COL_NAMES') which contain SNP annotation data
    splitvals = { n:vals[colindex[n]].split(",") for n in SNP_COL_NAMES }

    # check whether all annotations contain the same number of elements
    # by making a 'set' of all lengths: a set contains only uniq values so should have length 1 itself
    # if columns contain different number of values discard the entire data set and issue a warning
    if len(set([ len(splitvals[e]) for e in splitvals ])) != 1:
        logger.warning("varying number of SNP annotation data found in different columns, in line:\n\t%s" % "\t".join(vals))
        return []

    # check whether SNPs are duplicated in this fragment (due to overlap with both forward and reverse read)
    # remove the duplicated reads from the reverse reads as these are expected to be of worse quality
    def discard_one_of_SNPs(vals,s1,s2):
        return s2
    # loop over the vector SNP_ABS_POS in splitvals (these values should be
    #   unique for all SNP's in the fragment). Thus, duplicated values indicate
    #   duplicated SNP's
    # use dictionary 'retainIdx' to check if a abs_pos is duplicated
    #   if duplicated do not add this SNP index to the retainIdx dictionary
    retainIdx = {}
    for i in range(len(splitvals['SNP_ABS_POS'])):
        if splitvals['SNP_ABS_POS'][i] in retainIdx:
            # do not add 2nd copy of this SNP to dictionary 'retainIdx'
            continue
        else:
            retainIdx[splitvals['SNP_ABS_POS'][i]] = i
    # init a new list which'll contain the split data
    newvals = []
    # loop over split values; copy original data and replace convoluted SNP data with single values
    for i in retainIdx.values():
        v = vals.copy()
        for n in SNP_COL_NAMES:
            v[colindex[n]] = splitvals[n][i]
        newvals.append(v)

    # return split data
    return newvals

def discardSNP_ambiguous(vals,colindex):
    # discard all data for which the SNP_PARENT is *not* 'paternal' or 'maternal'
    if len(vals)==0: # check whether input is empty already
        return vals
    for i in reversed(range(len(vals))): # 'reversed(..)', so we can delete elements along the way
        if vals[i][colindex['SNP_PARENT']] != 'expected':
            del vals[i]
    return vals

def iterate_input(input_fname, output_fname):
    # open output, as output
    with gzip.open(output_fname,'wt') as output:
        with gzip.open(input_fname,'rt') as input:
            # read header
            header = input.readline().rstrip()
            output.write(header+"\n")
            colnames = header.split("\t")
            colindex = {}
            for i in range(len(colnames)):
                colindex[colnames[i]] = i
            for line in input:
                vals = line.rstrip().split("\t")
                vals = discardFragment_no_SNP(vals, colindex)
                vals = discardFragment_NA(vals, colindex)
                vals = deconvolute_SNPs(vals, colindex)
                vals = discardSNP_ambiguous(vals,colindex)
                for v in vals:
                    output.write("\t".join(v)+"\n")
                #process_line(line.rstrip(), output, colindex)

def main(options):
    # open log file
    logger = init_logger(options.log)
    logger.info("processing started")
    # process input
    iterate_input(options.input,options.output)
    logger.info("processing done")
    return


if __name__ == "__main__":
    sys.stderr.write("command: %s\n" % " ".join(sys.argv))
    options = parse_options()
    main(options)

