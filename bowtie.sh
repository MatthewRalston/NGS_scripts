#!/bin/bash
#PBS -N bwtie
#PBS -r n
#PBS -V
#PBS -l nodes=1:ppn=8
#PBS -l walltime=300:00:00
#PBS -d /home/mrals/ETP

set -e

. ~/.bash_profile
#General
PATH=$PATH:/home/mrals/pckges/Vienna/bin:/home/mrals/home/bin/:/home/mrals/bin/
#NGS
PATH=$PATH:/home/mrals/pckges/sratoolkit.2.3.4-2-centos_linux64/bin/:/usr/local/bowtie-0.12.7/:/usr/local/bowtie2-2.1.0/:/usr/local/bwa-0.7.4/:/usr/local/cufflinks-2.0.2/:/usr/local/FastQC/:/usr/local/samtools-0.1.18/:/usr/local/tophat-2.0.4/:/home/mrals/pckges/seqtk-master/

export PATH
[[ -s "$HOME/.rvm/scripts/rvm" ]] && source "$HOME/.rvm/scripts/rvm"
rvm use 2.0.0
CORES=8
INDIR=processed_fastq/processed
OUTDIR=SAM_unprocessed
FINALFASTQ=final/finalfastq
FINALQC=final/finalqc
PICARD=/usr/local/picard-tools-1.67
CHROMOSOMES='CAC.fa, CAP.fa'
BT2BASE=CAC
# Build indexes if necessary
#bowtie2-build -f $CHROMOSOMES $BT2BASE
REFERENCE=reference/$BT2BASE
rRNA=reference/rRNA
PHRED=33
SUFFIX='.bam'
FILES=(`/usr/bin/ls $INDIR`)
# This script runs bowtie on both paired and unpaired reads, aligning them first to the 
# rRNA sequences. After this in-silico rRNA-read removal, the remaining reads are then 
# aligned to the genome. FastQC checks are performed on the rRNA-filtered reads.


for ((i=0;i<${#FILES[@]}; i+=3))
do
    #mkdir $FINALQC/${FILES[$i]%_*} tmp/${FILES[$i]%_*}
    echo "$i paired rRNA $INDIR/${FILES[$i]} -> $FINALFASTQ/${FILES[$i]%_*}.fastq.gz"
    bowtie2 -p $CORES --mp 20,5 --no-mixed --fr --dpad 0 --fast --no-discordant -x $rRNA --un-conc-gz $FINALFASTQ/${FILES[$i]%_*}.fastq.gz -1 $INDIR/${FILES[$i]} -2 $INDIR/${FILES[$i+1]} -S /dev/null
    echo "$i unpaired rRNA $INDIR/${FILES[$i+3]} -> $FINALFASTQ/${FILES[$i+3]}"
    bowtie2 -p $CORES --mp 20,5 --dpad 0 --fast -x $rRNA --un-gz $FINALFASTQ/${FILES[$i+2]} -U $INDIR/${FILES[$i+2]} -S /dev/null
    echo "$i alignment"
    bowtie2 -p $CORES --very-sensitive -x $REFERENCE -1 $FINALFASTQ/${FILES[$i]%_*}.fastq.1.gz -2 $FINALFASTQ/${FILES[$i]%_*}.fastq.2.gz -U $FINALFASTQ/${FILES[$i+2]} -S /dev/stdout | samtools view -bhS - > $OUTDIR/${FILES[$i]%_*}$SUFFIX
    #fastqc -j /usr/bin/java -f fastq -o $FINALQC/${FILES[$i]%_*} $FINALFASTQ/${FILES[$i]%_*}.fastq.1.gz
    #fastqc -j /usr/bin/java -f fastq -o $FINALQC/${FILES[$i+1]%_*} $FINALFASTQ/${FILES[$i]%_*}.fastq.2.gz
    #fastqc -j /usr/bin/java -f fastq -o $FINALQC/${FILES[$i+3]} $FINALFASTQ/${FILES[$i+3]}
	#zcat $FINALFASTQ/$f | fastx_quality_stats -Q $PHRED > $FINALQC/${f%.*.*}/fastx_report.txt
	#zcat $FINALFASTQ/$f | prinseq -fastq stdin -stats_all > $FINALQC/${f%.*.*}/prinseq_stats.txt
done

#qsub postprocess.sh

## EOF-------------------------------------------
