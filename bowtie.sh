#!/bin/bash
#PBS -N bwtie
#PBS -r n
#PBS -V
#PBS -l nodes=1:ppn=16
#PBS -l walltime=300:00:00
#PBS -d /home/mrals/Final

#------------------------------------------------
# Title: bowtie.sh
#
# Copyright 2014 Matt Ralston
# 
# This file performs in silico ribosomal RNA 
# removal and then performs a reference alignment
# to the C. acetobutylicum genome. If desirable,
# final fastqc checks can be performed to check
# the quality of final rRNA-free reads. If this
# is done, however, the alignmentsummary.rb
# script will fail.
# 
#------------------------------------------------

#------------------------------------------------
# Functions
#------------------------------------------------
source functions.sh


#------------------------------------------------
# Parameters
#------------------------------------------------
# CORES: this is the number of processors requested
# for this job.
CORES=16
# INDIR: this is the directory which contains the
# reads for rRNA removal and reference alignment.
INDIR=processed/fastq
# OUTDIR: this is the directory where the unprocessed
# alignment files will be produced
OUTDIR=SAM_unprocessed
# FINALFASTQ: this is the directory where the final
# rRNA-free fastq files will be produced
export FINALFASTQ=final/finalfastq
# FINALQC: this is the directory where the final
# fastqc output files will be produced.
export FINALQC=final/finalqc
# BT2BASE: This is the prefix of the bowtie2 indexes
# that will be used for the reference alignement
BT2BASE=CAC
# CHROMOSOMES: these are files of fasta sequences that can
# be used to build bowtie2 indexes, if necessary
#CHROMOSOMES='CAC.fa, CAP.fa'
# Build indexes if necessary
#cd reference; bowtie2-build -f $CHROMOSOMES $BT2BASE; cd ..
# REFERENCE: this is the location of the bowtie2 reference
# indexes for the alignment
REFERENCE=reference/$BT2BASE
# rRNA: this is the location of the ribosomal RNA indexes
# that will be used to perform the in silico rRNA removal.
rRNA=reference/rRNA
# FILTER: This is a gtf format of the rRNA sequences
FILTER=reference/mask.gtf
# CONTAM: Store rRNA alignments
CONTAM=rrna
# PHRED: this is the PHRED offset that would be used for fastx toolkit
# summary statistics.
PHRED=33
# INTERMEDIATE: This is a location to store alignment files for filtering and statistics
INTERMEDIATE=Fullsam
# SUFFIX: this is the suffix that will be used to label the output files of
# the reference alignment, in this case, bam.
SUFFIX='.bam'
FILES=(`/usr/bin/ls $INDIR`)
# This script runs bowtie on both paired and unpaired reads, aligning them first to the 
# rRNA sequences. After this in-silico rRNA-read removal, the remaining reads are then 
# aligned to the genome. FastQC checks are performed on the rRNA-filtered reads.


#------------------------------------------------
# In-silico rRNA removal + Bowtie
#------------------------------------------------
for ((i=0; i < ${#FILES[@]}; i+=3))
do
    echo ${FILES[$i]} | ruby -e "puts gets.chomp.split('_')[0]"  >&2
   # O L D    B O W T I E   P I P E L I N E
    bowtie2 -p $CORES --mp 20,5 --no-mixed --fr --dpad 0 --fast --no-discordant -x $rRNA --un-conc-gz $FINALFASTQ/${FILES[$i]%_*}.fastq.gz -1 $INDIR/${FILES[$i]} -2 $INDIR/${FILES[$i+1]} -S /dev/null
    bowtie2 -p $CORES --mp 20,5 --dpad 0 --fast -x $rRNA --un-gz $FINALFASTQ/${FILES[$i+2]} -U $INDIR/${FILES[$i+2]} -S /dev/null
    bowtie2 -p $CORES --very-sensitive -x $REFERENCE -1 $FINALFASTQ/${FILES[$i]%_*}.fastq.1.gz -2 $FINALFASTQ/${FILES[$i]%_*}.fastq.2.gz -U $FINALFASTQ/${FILES[$i+2]} -S /dev/stdout | samtools view -bhS - > $OUTDIR/${FILES[$i]%_*}$SUFFIX


    #bowtie2 -p $CORES -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 --dpad 0 -x $REFERENCE -1 $INDIR/${FILES[$i]} -2 $INDIR/${FILES[$i+1]} -U $INDIR/${FILES[$i+2]} -S /dev/stdout | samtools view -bhS - > $INTERMEDIATE/${FILES[$i]%_*}$SUFFIX
    #bedtools intersect -v -abam $INTERMEDIATE/${FILES[$i]%_*}$SUFFIX -b $FILTER > $OUTDIR/${FILES[$i]%_*}$SUFFIX
    #bedtools intersect -wa -abam $INTERMEDIATE/${FILES[$i]%_*}$SUFFIX -b $FILTER > $CONTAM/${FILES[$i]%_*}$SUFFIX
    
    
done

#------------------------------------------------
# Quality
#------------------------------------------------
FINALFILES=`/usr/bin/ls $FINALFASTQ`

qsub postprocess.sh
parallel -j $CORES 'quality {} $FINALFASTQ $FINALQC' ::: $FINALFILES


## EOF-------------------------------------------
