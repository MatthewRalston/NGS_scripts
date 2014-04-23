#!/bin/bash
#PBS -N postprocessing
#PBS -r n
#PBS -V
#PBS -l nodes=1:ppn=1
#PBS -l walltime=72:00:00,cput=72:00:00
#PBS -d /home/mrals/ETP

#set -e

. ~/.bash_profile
#General
PATH=$PATH:/home/mrals/pckges/Vienna/bin:/home/mrals/home/bin/:/home/mrals/bin/
#NGS
PATH=$PATH:/home/mrals/pckges/sratoolkit.2.3.4-2-centos_linux64/bin/:/usr/local/bowtie-0.12.7/:/usr/local/bowtie2-2.1.0/:/usr/local/bwa-0.7.4/:/usr/local/cufflinks-2.0.2/:/usr/local/FastQC/:/usr/local/samtools-0.1.18/:/usr/local/tophat-2.0.4/:/home/mrals/pckges/seqtk-master/

export PATH
[[ -s "$HOME/.rvm/scripts/rvm" ]] && source "$HOME/.rvm/scripts/rvm"
rvm use 2.0.0

#########################
#  D E S C R I P T I O N
#########################
#  AUTHOR: MATT RALSTON
#  DATE: 1/22/14
#
#  This script is capable of postprocessing .sam files.
#  1. First, the SAM file is cleaned and verified
#  2. Second, the SAM file is sorted, and converted to BAM
#       By default, unmapped reads are removed during this step.
#  3. Third, duplicates are marked in the BAM file. A BAM index is created After this step and the final BAM file is validated.
#	   Additionally, summary statistics are reported for the bam file and its degree of duplication is reported.
#  4. OPTIONALLY, additional summary files can be created for
#		RNAseq metrics (strand specificity, etc.)
#		Illumina TruSeq librarys
#		Insert size metrics
#        
#	This is a Torque submission BASH script that runs in the directory described in the Torque header's directory option: -d
#   Inside this directory, Picard Tools and SAMtools are expected to run. The location of Picard-tools .jar files is specified
#   by the option PICARD. 
#
#######################################################################
#  N O T E :    This script will run with 4 mandatory arguments (below)
#######################################################################
#  Here, we distinguish a filepath as the name of a directory (without the terminal backslash "/")
#  And a filename as the relative or full filename of a file.
#   1. PICARD: The filepath of the picard install, containing the .jar files.
#   2. REFERENCE: The filename of the reference genome, in FASTA format.
#   3. FILES: The string-list of filename prefixes, each on a separate line, without the terminal ".sam" suffix.
#			e.g. Two .sam files: samfiles/fileone.sam and samfiles/filetwo.sam, FILES and SAMDIR would be:
#
#  					FILES="fileone
#				    filetwo"
#					SAMDIR=samfiles
#
#	4. SAMDIR: Ths filepath of the directory which contains FILES.
#	5. TMPD: The filpath of a directory for temporary files.
#	6. LOGDIR: The (OPTIONAL) filepath to a directory for .log files (for summary statistics).
#			If left unspecified, the resulting .log files will spill to the working directory (Torque option -d above)

#########################
# ARGUMENTS
#########################
PICARD=/usr/local/picard-tools-1.67
REFERENCE=reference/CAC.txt
# Optional for RNAseq: RefFlat file
# REFFLAT=reference/H.pylori.refFlat
INDIR=SAM_unprocessed
OUTDIR=SAM_processed
FILES=`/usr/bin/ls $INDIR`
TMPD=/home/mrals/ETP/tmp/
LOGDIR=logs
REFFLAT=reference/CAC.refflat

#########################
# S C R I P T
#########################


for file in $FILES
do
	############   1   ################
	# CleanSam: Reads a SAM file, soft-clips alignments that overhang the end of the reference sequence (e.g. for circular bact. chromosomes)
	#              also, sets unmapped reads MAPQ to 0, if not already.
	samtools view -h -o /dev/stdout $INDIR/$file | java -jar $PICARD/CleanSam.jar INPUT=/dev/stdin OUTPUT=$OUTDIR/${file%.*}.1.sam

	# ValidateSamFile: Checks the validity of a SAM/BAM file, prints the results to cleansam.log
	java -jar $PICARD/ValidateSamFile.jar INPUT=$OUTDIR/${file%.*}.1.sam OUTPUT=$LOGDIR/${file%.*}.1_sam_validation.log

	############   2   ################
	# SAM->BAM + sort: Converts SAM to BAM format and sorts, removing unmapped reads. An alternative is provided which does not filter unmapped reads.
	#samtools view -hubS $file.1.sam | samtools sort - $file.2
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
	java -jar $PICARD/CollectInsertSizeMetrics.jar INPUT=$OUTDIR/${file%.*}.3.bam  HISTOGRAM_FILE=$LOGDIR/${file%.*}.3.insert_size_histogram.jpg OUTPUT=$LOGDIR/${file%.*}.3.insert_size_stats.log REFERENCE_SEQUENCE=$REFERENCE ASSUME_SORTED=true TMP_DIR=$TMPD/
done






## EOF-------------------------------------------
