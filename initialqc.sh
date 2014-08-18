#!/bin/bash
#PBS -N initialqc
#PBS -r n
#PBS -V
#PBS -l nodes=1:ppn=6
#PBS -l walltime=32:00:00
#PBS -d /home/mrals/Final
#------------------------------------------------
# Title: initialqc.sh
#
# Copyright 2014 Matt Ralston
# 
# This script performs a quality check on fastq
# format files in a directory, reporting the 
# results to a separate directory.
# 
#------------------------------------------------

#------------------------------------------------
# Functions
#------------------------------------------------
source functions.sh


#------------------------------------------------
# Parameters
#------------------------------------------------
# CORES: the number of processors available
CORES=6
# INDIR: the name of a directory that contains
#        fastq or fastq.gz format reads.
export INDIR=rawdata
# OUTDIR: the name of a directory where the 
#         results will be reported
export OUTDIR=initialqc
# PHRED: the phred offset/quality encoding of the
#        reads
PHRED=33
# FILES: An array of the names of the files to process.

FILES=`/usr/bin/ls $INDIR`


#------------------------------------------------
# Quality
#------------------------------------------------
parallel -j $CORES 'quality {} $INDIR $OUTDIR' ::: $FILES



#qsub trimfiltercheck.sh
## EOF-------------------------------------------
