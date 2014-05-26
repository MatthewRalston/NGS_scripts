#!/bin/bash
#PBS -N bwtie
#PBS -r n
#PBS -V
#PBS -l nodes=1:ppn=8
#PBS -l walltime=300:00:00
#PBS -d /home/mrals/ETP

#------------------------------------------------
# Title: bowtie.sh
#
# Matt Ralston
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
# Parameters
#------------------------------------------------
# CORES: this is the number of processors requested
# for this job.
CORES=8
# INDIR: this is the directory which contains the
# reads for rRNA removal and reference alignment.
INDIR=processed_fastq/processed
# OUTDIR: this is the directory where the unprocessed
# alignment files will be produced
OUTDIR=SAM_unprocessed
# FINALFASTQ: this is the directory where the final
# rRNA-free fastq files will be produced
FINALFASTQ=final/finalfastq
# FINALQC: this is the directory where the final
# fastqc output files will be produced.
FINALQC=final/finalqc
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
# PHRED: this is the PHRED offset that would be used for fastx toolkit
# summary statistics.
PHRED=33
# SUFFIX: this is the suffix that will be used to label the output files of
# the reference alignment, in this case, bam.
SUFFIX='.bam'
FILES=(`/usr/bin/ls $INDIR`)
# This script runs bowtie on both paired and unpaired reads, aligning them first to the 
# rRNA sequences. After this in-silico rRNA-read removal, the remaining reads are then 
# aligned to the genome. FastQC checks are performed on the rRNA-filtered reads.


for ((i=0;i<${#FILES[@]}; i+=3))
do
    mkdir $FINALQC/${FILES[$i]%_*} tmp/${FILES[$i]%_*} $FINALQC/${FILES[$i+1]%_*} $FINALQC/${FILES[$i+2]%.gz}
    echo ${FILES[$i]*.*.*} >&2
    bowtie2 -p $CORES --mp 20,5 --no-mixed --fr --dpad 0 --fast --no-discordant -x $rRNA --un-conc-gz $FINALFASTQ/${FILES[$i]%_*}.fastq.gz -1 $INDIR/${FILES[$i]} -2 $INDIR/${FILES[$i+1]} -S /dev/null
    bowtie2 -p $CORES --mp 20,5 --dpad 0 --fast -x $rRNA --un-gz $FINALFASTQ/${FILES[$i+2]} -U $INDIR/${FILES[$i+2]} -S /dev/null
    #bowtie2 -p $CORES --very-sensitive -x $REFERENCE -1 $FINALFASTQ/${FILES[$i]%_*}.fastq.1.gz -2 $FINALFASTQ/${FILES[$i]%_*}.fastq.2.gz -U $FINALFASTQ/${FILES[$i+2]} -S /dev/stdout | samtools view -bhS - > $OUTDIR/${FILES[$i]%_*}$SUFFIX
    #fastqc -j /usr/bin/java -f fastq -o $FINALQC/${FILES[$i]%_*} $FINALFASTQ/${FILES[$i]%_*}.fastq.1.gz
    #fastqc -j /usr/bin/java -f fastq -o $FINALQC/${FILES[$i+1]%_*} $FINALFASTQ/${FILES[$i]%_*}.fastq.2.gz
    #fastqc -j /usr/bin/java -f fastq -o $FINALQC/${FILES[$i+1]%_*} $FINALFASTQ/${FILES[$i+2]}
	#zcat $FINALFASTQ/$f | fastx_quality_stats -Q $PHRED > $FINALQC/${f%.*.*}/fastx_report.txt
	#zcat $FINALFASTQ/$f | prinseq -fastq stdin -stats_all > $FINALQC/${f%.*.*}/prinseq_stats.txt
done

#qsub postprocess.sh

## EOF-------------------------------------------
