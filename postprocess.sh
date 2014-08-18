#!/bin/bash
#PBS -N postprocessing
#PBS -r n
#PBS -V
#PBS -l nodes=1:ppn=12
#PBS -l walltime=240:00:00
#PBS -d /home/mrals/NichSandoval
#------------------------------------------------
# Title: postprocess.sh
#
# Copyright 2014 Matt Ralston
# 
# This script performs some housekeeping
# operations on alignment files.
# First, each alignment is cleaned and sorted,
# duplicates are marked, the file is validated
# and indexed. Additionally, several summary
# statistics are generated for each alignment.
# 
#------------------------------------------------


#set -e

#------------------------------------------------
# Parameters
#------------------------------------------------
#CORES: This is the number of processors available for parallelization
CORES=12
# PICARD: this is the location of the jar files
# for the picard suite
export PICARD=/usr/local/picard-tools-1.67
# REFERENCE: this is the fasta file of the
# C. ac. genome for InsertSizeMetrics
export REFERENCE=reference/CAC.txt
# INDIR: this is the location of the unprocessed
# alignment files.
export INDIR=SAM_unprocessed
# OUTDIR: this is the directory where processed
# alignment files will be produced with the
# suffix *.3.bam
export OUTDIR=SAM_processed
# TMPD: this is a directory for temporary files
export TMPD=/home/mrals/NichSandoval/tmp/
# LOGDIR: this is a directory for log files to be written
export LOGDIR=logs
FILES=`/usr/bin/ls $INDIR`


function postprocess {
    file=$1
    ############   1   ################
    # CleanSam: Reads a SAM file, soft-clips alignments that overhang the end of the reference sequence (e.g. for circular bact. chromosomes)
    #              also, sets unmapped reads MAPQ to 0, if not already.
    samtools view -h -q 5 -o /dev/stdout $INDIR/$file | java -jar $PICARD/CleanSam.jar INPUT=/dev/stdin OUTPUT=$OUTDIR/${file%.*}.1.sam

    # ValidateSamFile: Checks the validity of a SAM/BAM file, prints the results to cleansam.log
    java -jar $PICARD/ValidateSamFile.jar INPUT=$OUTDIR/${file%.*}.1.sam OUTPUT=$LOGDIR/${file%.*}.1_sam_validation.log

    ############   2   ################
    # SAM->BAM + sort: Converts SAM to BAM format and sorts, removing unmapped reads. An alternative is provided which does not filter unmapped reads.
    samtools view -hubS $file.1.sam | samtools sort - $file.2
    samtools view -hubS -F 4 $OUTDIR/${file%.*}.1.sam | samtools sort - $OUTDIR/${file%.*}.2

    ############   3   ################
    # MarkDuplicates: Marks duplicate reads as such in bam file, reports metrics to log file.
    java -jar $PICARD/MarkDuplicates.jar INPUT=$OUTDIR/${file%.*}.2.bam METRICS_FILE=$LOGDIR/${file%.*}.2.duplication.log OUTPUT=$OUTDIR/${file%.*}.3.bam ASSUME_SORTED=true TMP_DIR=$TMPD

    # CollectAlignmentSummaryMetrics
    java -jar $PICARD/CollectAlignmentSummaryMetrics.jar INPUT=$OUTDIR/${file%.*}.3.bam OUTPUT=$LOGDIR/${file%.*}.3.summary.log REFERENCE_SEQUENCE=$REFERENCE ASSUME_SORTED=true TMP_DIR=$TMPD

    # Create bam index
    samtools index $OUTDIR/${file%.*}.3.bam $OUTDIR/${file%.*}.3.bam.bai TMP_DIR=$TMPD

    # Validate final bam file
    java -jar $PICARD/ValidateSamFile.jar INPUT=$OUTDIR/${file%.*}.3.bam OUTPUT=$LOGDIR/${file%.*}.3_bam_validation.log TMP_DIR=$TMPD
	
    # PAIRED END ONLY: CollectInsertSizeMetrics
    java -jar $PICARD/CollectInsertSizeMetrics.jar INPUT=$OUTDIR/${file%.*}.3.bam  HISTOGRAM_FILE=$LOGDIR/${file%.*}.3.insert_size_histogram.pdf OUTPUT=$LOGDIR/${file%.*}.3.insert_size_stats.log REFERENCE_SEQUENCE=$REFERENCE ASSUME_SORTED=true TMP_DIR=$TMPD/
}
export -f postprocess
parallel -j $CORES 'postprocess {}' ::: $FILES



#qsub trinity.sh


## EOF-------------------------------------------
