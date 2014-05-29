#!/bin/bash
#PBS -N initialqc
#PBS -r n
#PBS -V
#PBS -l nodes=1:ppn=1
#PBS -l walltime=32:00:00
#PBS -d /home/mrals/ETP
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
# Parameters
#------------------------------------------------
# INDIR: the name of a directory that contains
#        fastq or fastq.gz format reads.
INDIR=rawdata
# OUTDIR: the name of a directory where the 
#         results will be reported
OUTDIR=initialqc
# PHRED: the phred offset/quality encoding of the
#        reads
PHRED=33

FILES=`/usr/bin/ls $INDIR`
for f in $FILES
do
	mkdir ${OUTDIR}/${f%.*.*}
	fastqc -j /usr/bin/java -f fastq -o ${OUTDIR}/${f%.*.*} $INDIR/$f
	zcat $INDIR/$f | fastx_quality_stats -Q $PHRED  > $OUTDIR/${f%.*.*}/fastx_report.txt
	zcat $INDIR/$f | prinseq -fastq stdin -stats_all > $OUTDIR/${f%.*.*}/prinseq_stats.txt
done



## EOF-------------------------------------------
