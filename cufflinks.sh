#!/bin/bash
#PBS -N assembly
#PBS -r n
#PBS -V
#PBS -l nodes=1:ppn=2
#PBS -l walltime=40:00:00,cput=80:00:00
#PBS -d /home/mrals/test

. ~/.bash_profile
#General
PATH=$PATH:/home/mrals/pckges/Vienna/bin:/home/mrals/home/bin/:/home/mrals/bin/
#NGS
PATH=$PATH:/home/mrals/pckges/sratoolkit.2.3.4-2-centos_linux64/bin/:/usr/local/bowtie-0.12.7/:/usr/local/bowtie2-2.1.0/:/usr/local/bwa-0.7.4/:/usr/local/cufflinks-2.0.2/:/usr/local/FastQC/:/usr/local/samtools-0.1.18/:/usr/local/tophat-2.0.4/:/home/mrals/pckges/seqtk-master/

export PATH
[[ -s "$HOME/.rvm/scripts/rvm" ]] && source "$HOME/.rvm/scripts/rvm"
rvm use 2.0.0

INDIR=SAM_processed
CUFFLINKS=Cufflinks_assemblies
CUFFCOMPARE
REFERENCE=reference/CAC.gff
FILES=`/usr/bin/ls $INDIR/*.3.bam`
for f in $FILES
do
	mkdir $OUTDIR/${f%.*.*}
	cufflinks -o $CUFFLINKS/${f%.*.*} -p 2 -g $REFERENCE --3-overhang-tolerance 50 $f
	cuffcompare -r $REFERENCE -o ${f%.*.*} $CUFFLINKS/${f%.*.*}/transcripts.gtf
done



## EOF-------------------------------------------
