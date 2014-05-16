#!/bin/bash
#PBS -N TFC
#PBS -r n
#PBS -V
#PBS -l nodes=1:ppn=2
#PBS -l walltime=120:00:00
#PBS -d /home/mrals/ETP
#------------------------------------------------
# Title: trimfiltercheck.sh
#
# Matt Ralston
# 
# This script performs trimming with trimmomatic
# and then performs a quality check with fastqc.
# 
# 
#------------------------------------------------




#------------------------------------------------
# Parameters
#------------------------------------------------

# INDIR is the location of the unprocessed fastq
# reads that will be trimmed.
INDIR=rawdata
# OUTDIR is the location where the results will
# be deposited
OUTDIR=processed_fastq
mkdir $OUTDIR/processed $OUTDIR/qc
# QC is the location where fastqc reports will
# be generated.
QC=$OUTDIR/qc
# OUTPUT is the location where final fastq files
# will be deposited
OUTPUT=$OUTDIR/processed
FILES=(`/usr/bin/ls $INDIR`)
# PHRED encoding/offset of input files (33 sanger etc.)
PHRED=33
# WINDOW:
WINDOW=4
# AVGQ:
20
# MINLEN:
30


for ((i=0;i<${#FILES[@]};i+=2))
do
	java -jar /home/mrals/pckges/Trimmomatic-0.32/trimmomatic-0.32.jar PE -threads 6 -phred33 -trimlog logs/${FILES%_*}.trimmomatic.log $INDIR/${FILES[$i]} $INDIR/${FILES[$i+1]} $OUTPUT/${FILES[$i]} $OUTPUT/${FILES[$i]}.for.unpaired.gz $OUTPUT/${FILES[$i+1]} $OUTPUT/${FILES[$i+1]}.rev.unpaired.gz ILLUMINACLIP:/home/mrals/pckges/Trimmomatic-0.32/adapters/TruSeq3-PE.fa:1:40:12:15 LEADING:5 TRAILING:5 SLIDINGWINDOW:$WINDOW:$AVGQ MINLEN:$MINLEN
	# here, the unpaired reads from the reverse strand of an unpaired read
	# are reverse complemented, so that all unpaired reads are strand specific
	# with respect to the forward strand.
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
