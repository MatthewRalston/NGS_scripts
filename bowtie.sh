#!/bin/bash
#PBS -N alignment
#PBS -r n
#PBS -V
#PBS -l nodes=1:ppn=4
#PBS -l walltime=40:00:00,cput=160:00:00
#PBS -d /home/mrals/test

. ~/.bash_profile
#General
PATH=$PATH:/home/mrals/pckges/Vienna/bin:/home/mrals/home/bin/:/home/mrals/bin/
#NGS
PATH=$PATH:/home/mrals/pckges/sratoolkit.2.3.4-2-centos_linux64/bin/:/usr/local/bowtie-0.12.7/:/usr/local/bowtie2-2.1.0/:/usr/local/bwa-0.7.4/:/usr/local/cufflinks-2.0.2/:/usr/local/FastQC/:/usr/local/samtools-0.1.18/:/usr/local/tophat-2.0.4/:/home/mrals/pckges/seqtk-master/

export PATH
[[ -s "$HOME/.rvm/scripts/rvm" ]] && source "$HOME/.rvm/scripts/rvm"
rvm use 2.0.0

INDIR=one/processed/
OUTDIR=SAM_unprocessed
FINALFASTQ=finalfastq
FINALQC=finalqc
PICARD=/usr/local/picard-tools-1.67
REFERENCE=reference/CAC
rRNA=reference/rRNA
PHRED=33
SUFFIX='.bam'
FILES=`/usr/bin/ls $INDIR`
for f in $FILES
do
	mkdir $FINALQC/${f%.*.*}
	bowtie2 -p 4 --very-sensitive -x $rRNA --un-gz $FINALFASTQ/$f -U $INDIR/$f -S /dev/null
	bowtie2 -p 4 -N 1 --very-sensitive -x $REFERENCE --un-gz tmp/$f -U $FINALFASTQ/$f -S /dev/stdout | samtools view -bhS - > $OUTDIR/${f%.*.*}$SUFFIX
	
	fastqc -j /usr/bin/java -f fastq -o $FINALQC/${f%.*.*} $FINALFASTQ/$f
	zcat $FINALFASTQ/$f | fastx_quality_stats -Q $PHRED > $FINALQC/${f%.*.*}/fastx_report.txt
	zcat $FINALFASTQ/$f | prinseq -fastq stdin -stats_all > $FINALQC/${f%.*.*}/prinseq_stats.txt
done

#qsub postprocess.sh

## EOF-------------------------------------------
