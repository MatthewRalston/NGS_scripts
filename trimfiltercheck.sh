#!/bin/bash
#PBS -N TFC
#PBS -r n
#PBS -V
#PBS -l nodes=1:ppn=2
#PBS -l walltime=120:00:00
#PBS -d /home/mrals/ETP

. ~/.bash_profile
#General
PATH=$PATH:/home/mrals/pckges/Vienna/bin:/home/mrals/home/bin/:/home/mrals/bin/
#NGS
PATH=$PATH:/home/mrals/pckges/sratoolkit.2.3.4-2-centos_linux64/bin/:/usr/local/bowtie-0.12.7/:/usr/local/bowtie2-2.1.0/:/usr/local/bwa-0.7.4/:/usr/local/cufflinks-2.0.2/:/usr/local/FastQC/:/usr/local/samtools-0.1.18/:/usr/local/tophat-2.0.4/:/home/mrals/pckges/seqtk-master/

export PATH
[[ -s "$HOME/.rvm/scripts/rvm" ]] && source "$HOME/.rvm/scripts/rvm"
rvm use 2.0.0

INDIR=rawdata
OUTDIR=processed_fastq
mkdir $OUTDIR/processed $OUTDIR/qc
QC=$OUTDIR/qc
OUTPUT=$OUTDIR/processed
FILES=(`/usr/bin/ls $INDIR`)
# PHRED encoding/offset of input files (33 sanger etc.)
PHRED=33

# This script contains a rather lengthy call to trimmomatic to trim quality 
# and clip adapters so that there is good quality in paired end reads. 
#The TruSeq v3 adapters used in the ScriptSeq (Epicentre) v2 kit are clipped from the reads
# Resulting unpaired reads are combined into a single file, after reverse complementing
# reads from the second mate.
# FastQC quality checks are then run on the resulting files.



for ((i=0;i<${#FILES[@]};i+=2))
do
	java -jar /home/mrals/pckges/Trimmomatic-0.32/trimmomatic-0.32.jar PE -threads 6 -phred33 -trimlog logs/${FILES%_*}.trimmomatic.log $INDIR/${FILES[$i]} $INDIR/${FILES[$i+1]} $OUTPUT/${FILES[$i]} $OUTPUT/${FILES[$i]}.for.unpaired.gz $OUTPUT/${FILES[$i+1]} $OUTPUT/${FILES[$i+1]}.rev.unpaired.gz ILLUMINACLIP:/home/mrals/pckges/Trimmomatic-0.32/adapters/TruSeq3-PE.fa:1:50:12:15 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:30

	gunzip -c $OUTPUT/${FILES[$i+1]}.rev.unpaired.gz | fastx_reverse_complement -Q33 | gzip > $OUTPUT/${FILES[$i+1]}.rev.gz
	cat $OUTPUT/${FILES[$i+1]}.rev.gz $OUTPUT/${FILES[$i]}.for.unpaired.gz > $OUTPUT/${FILES[$i]%_*}.unpaired.gz
	rm $OUTPUT/${FILES[$i+1]}.rev.gz $OUTPUT/${FILES[$i+1]}.rev.unpaired.gz $OUTPUT/${FILES[$i]}.for.unpaired.gz
		
	fastqc -j /usr/bin/java -f fastq -o $QC/${FILES[$i]%_*} $OUTPUT/${FILES[$i]}
	fastqc -j /usr/bin/java -f fastq -o $QC/${FILES[$i+1]%_*} $OUTPUT/${FILES[$i+1]}
	fastqc -j /usr/bin/java -f fastq -o $QC/${FILES[$i+1]%_*}.unpaired $OUTPUT/${FILES[$i]%_*}.unpaired.gz


#zcat $OUTPUT/${FILES[$i]} | fastx_quality_stats -Q $PHRED > $QC/${FILES[$i]%_*}/${FILES[$i]%.*}.fastx
	#zcat $OUTPUT/${FILES[$i]} | prinseq -fastq stdin -stats_all > $QC/${FILES[$i]%_*}/${FILES[$i]%.*}.prinseq
	#zcat $OUTPUT/${FILES[$i+1]} | fastx_quality_stats -Q $PHRED > $QC/${FILES[$i+1]%_*}/${FILES[$i+1]%.*}.fastx
	#zcat $OUTPUT/${FILES[$i+1]} | prinseq -fastq stdin -stats_all > $QC/${FILES[$i+1]%_*}/${FILES[$i+1]%.*}.prinseq
done




## EOF-------------------------------------------
