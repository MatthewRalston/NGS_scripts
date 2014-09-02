#!/bin/bash
#PBS -N TFC
#PBS -r n
#PBS -V
#PBS -l nodes=1:ppn=20
#PBS -l walltime=120:00:00
#PBS -d /home/mrals/Final
#------------------------------------------------
# Title: trimfiltercheck.sh
#
# Copyright 2014 Matt Ralston
# 
# This script performs trimming with trimmomatic
# and then performs a quality check with fastqc.
# 
# 
#------------------------------------------------


#------------------------------------------------
# Functions
#------------------------------------------------
source functions.sh

#------------------------------------------------
# Parameters
#------------------------------------------------

# CORES is the number of cores for parallelization of the trimming
CORES=20
# INDIR is the location of the unprocessed fastq
# reads that will be trimmed.
export INDIR=rawdata
# OUTDIR is the location where the results will
# be deposited
OUTDIR=processed
mkdir $OUTDIR/fastq $OUTDIR/qc
# QC is the location where fastqc reports will
# be generated.
export QC=$OUTDIR/qc
# OUTPUT is the location where final fastq files
# will be deposited
export OUTPUT=$OUTDIR/fastq
FILES=(`/usr/bin/ls $INDIR`)
# PHRED encoding/offset of input files (33 sanger etc.)
PHRED=33
# WINDOW:
WINDOW=4
# AVGQ:
AVGQ=20
# MINLEN:
MINLEN=30

#------------------------------------------------
# Pre-processing
#------------------------------------------------
# Each file is processed, concatenating the Casava 1.8+ header

#FILES=`/usr/bin/ls $INDIR`
#parallel -j $CORES 'preprocess {} $INDIR' ::: $FILES

#------------------------------------------------
# Trimming
#------------------------------------------------


for ((i=0;i<${#FILES[@]};i+=2))
do
	java -jar /home/mrals/pckges/Trimmomatic-0.32/trimmomatic-0.32.jar PE -threads $CORES -phred33 -trimlog logs/${FILES%_*}.trimmomatic.log $INDIR/${FILES[$i]} $INDIR/${FILES[$i+1]} $OUTPUT/${FILES[$i]} $OUTPUT/${FILES[$i]}.for.unpaired.gz $OUTPUT/${FILES[$i+1]} $OUTPUT/${FILES[$i+1]}.rev.unpaired.gz ILLUMINACLIP:/home/mrals/pckges/Trimmomatic-0.32/adapters/TruSeq3-PE.fa:1:40:12:15 LEADING:5 TRAILING:5 SLIDINGWINDOW:$WINDOW:$AVGQ MINLEN:$MINLEN
	# here, the unpaired reads from the reverse strand of an unpaired read
	# are reverse complemented, so that all unpaired reads are strand specific
	# with respect to the forward strand.
	gunzip -c $OUTPUT/${FILES[$i+1]}.rev.unpaired.gz | fastx_reverse_complement -Q33 | gzip > $OUTPUT/${FILES[$i+1]}.rev.gz
	cat $OUTPUT/${FILES[$i+1]}.rev.gz $OUTPUT/${FILES[$i]}.for.unpaired.gz > $OUTPUT/${FILES[$i]%_*}_unpaired.gz
	rm $OUTPUT/${FILES[$i+1]}.rev.gz $OUTPUT/${FILES[$i+1]}.rev.unpaired.gz $OUTPUT/${FILES[$i]}.for.unpaired.gz
done

FINALFILES=`/usr/bin/ls $OUTPUT`
parallel -j $CORES 'quality {} $OUTPUT $QC' ::: $FINALFILES

qsub bowtie.sh


## EOF-------------------------------------------
