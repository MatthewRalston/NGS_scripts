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
#########################
# T R I M M I N G
#########################
# Minimum quality to keep from sliding window trimming
q=20
# Minimum read size to keep
min=25

#########################
# F I L T E R
#########################
# Minimum Quality score to keep across 100% of the reads  
q100=20
# The fastq_quality filter can be retooled to filter differently, e.g. by setting the percentage of bases to carry that base. 



for ((i=0;i<${#FILES[@]};i+=2))
do
        #gunzip $INDIR/${FILES[$i]}; gunzip $INDIR/${FILES[$i+1]}


        #sickle pe -t sanger -q $q -l $min -f $INDIR/${FILES[$i]%.*}.clipped -r $INDIR/${FILES[$i+1]%.*}.clipped -o $OUTPUT/${FILES[$i]%.*} -p $OUTPUT/${FILES[$i+1]%.*}

	java -jar /home/mrals/pckges/Trimmomatic-0.32/trimmomatic-0.32.jar PE -threads 2 --phred33 -trimlog logs/${FILES%_*}.trimmomatic.log $INDIR/${FILES[$i]} $INDIR/${FILES[$i+1]} $OUTPUT/${FILES[$i]} /dev/null $OUTPUT/${FILES[$i+1]}.clipped /dev/null ILLUMINACLIP:/home/mrals/pckges/trimmomatic-0.32/adapters/ScriptSeq.fa:1:29:14 LEADING:20 TRAILING:15 SLIDINGWINDOW:4:20 MINLENGTH:30

	#cutadapt -a ACACTCTTTCCCTACACGACGCTCTTCCGATCT $INDIR/${FILES[$i]} > $INDIR/${FILES[$i]%.*}.clipped 
	#cutadapt -a GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT $INDIR/${FILES[$i+1]} > $INDIR/${FILES[$i+1]%.*}.clipped


	#fastq_quality_filter -Q $PHRED -q $q100 -p 100 -z > $OUTPUT/$file
	#gzip $INDIR/${FILES[$i]%.*}.clipped; gzip $INDIR/${FILES[$i+1]%.*}.clipped; gzip $OUTPUT/${FILES[$i]%.*}; gzip $OUTPUT/${FILES[$i+1]%.*}
	mkdir $QC/${FILES[$i]%_*} 
	fastqc -j /usr/bin/java -f fastq -o $QC/${FILES[$i]%_*} $OUTPUT/${FILES[$i]%.*}
	fastqc -j /usr/bin/java -f fastq -o $QC/${FILES[$i+1]%_*} $OUTPUT/${FILES[$i+1]%.*}
	zcat $OUTPUT/${FILES[$i]} | fastx_quality_stats -Q $PHRED > $QC/${FILES[$i]%_*}/${FILES[$i]%.*}.fastx
	zcat $OUTPUT/${FILES[$i]} | prinseq -fastq stdin -stats_all > $QC/${FILES[$i]%_*}/${FILES[$i]%.*}.prinseq
	zcat $OUTPUT/${FILES[$i+1]} | fastx_quality_stats -Q $PHRED > $QC/${FILES[$i+1]%_*}/${FILES[$i+1]%.*}.fastx
	zcat $OUTPUT/${FILES[$i+1]} | prinseq -fastq stdin -stats_all > $QC/${FILES[$i+1]%_*}/${FILES[$i+1]%.*}.prinseq
done




## EOF-------------------------------------------
