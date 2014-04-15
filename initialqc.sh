#!/bin/bash
#PBS -N initialqc
#PBS -r n
#PBS -V
#PBS -l nodes=1:ppn=1
#PBS -l walltime=32:00:00
#PBS -d /home/mrals/ETP

. ~/.bash_profile
#General
PATH=$PATH:/home/mrals/pckges/Vienna/bin:/home/mrals/home/bin/
#NGS
PATH=$PATH:/home/mrals/pckges/sratoolkit.2.3.4-2-centos_linux64/bin/:/usr/local/bowtie-0.12.7/:/usr/local/bowtie2-2.1.0/:/usr/local/bwa-0.7.4/:/usr/local/cufflinks-2.0.2/:/usr/local/FastQC/:/usr/local/samtools-0.1.18/:/usr/local/tophat-2.0.4/

export PATH
[[ -s "$HOME/.rvm/scripts/rvm" ]] && source "$HOME/.rvm/scripts/rvm"
rvm use 2.0.0


INDIR=rawdata
OUTDIR=initialqc
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
