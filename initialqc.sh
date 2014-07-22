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
# CORES: The number of cores to speed the preprocessing.
CORES=6
# FILES: An array of the names of the files to process.

FILES=`/usr/bin/ls $INDIR`


# Pre-processing
# First, each fastq file in unzipped, changing its ending
# from .gz to .fastq or .fq
parallel -j $CORES 'gunzip {}' ::: $FILES
# Then, these files are processed, concatenating the Casava 1.8+ header
# an compressing to the original file names
UNZIPD=`/usr/bin/ls $INDIR`
export MOD='puts($_.chomp.split.join("_"))'
parallel -j $CORES "cat {} | ruby -ne '$MOD' | gzip > {}.gz" ::: $UNZIPD
# Then, the uncompressed files are removed.
parallel -j $CORES "rm {}" ::: $UNZIPD

for f in $FILES
do

	mkdir ${OUTDIR}/${f%.*.*}
	fastqc -j /usr/bin/java -f fastq -o ${OUTDIR}/${f%.*.*} $INDIR/$f
	zcat $INDIR/$f | fastx_quality_stats -Q $PHRED  > $OUTDIR/${f%.*.*}/fastx_report.txt
	zcat $INDIR/$f | prinseq -fastq stdin -stats_all > $OUTDIR/${f%.*.*}/prinseq_stats.txt
done


qsub trimfiltercheck.sh
## EOF-------------------------------------------
