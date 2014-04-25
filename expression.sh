#!/bin/bash
#PBS -N gene_express
#PBS -r n
#PBS -V
#PBS -l nodes=1:ppn=6
#PBS -l walltime=240:00:00
#PBS -d /home/mrals/ETP

#set -e


. ~/.bash_profile
#General
PATH=$PATH:/home/mrals/pckges/Vienna/bin:/home/mrals/home/bin/:/home/mrals/bin/
#NGS
PATH=$PATH:/home/mrals/pckges/sratoolkit.2.3.4-2-centos_linux64/bin/:/usr/local/bowtie-0.12.7/:/usr/local/bowtie2-2.1.0/:/usr/local/bwa-0.7.4/:/usr/local/cufflinks-2.2.0/:/usr/local/FastQC/:/usr/local/samtools-0.1.18/:/usr/local/tophat-2.0.4/:/home/mrals/pckges/seqtk-master/

export PATH
[[ -s "$HOME/.rvm/scripts/rvm" ]] && source "$HOME/.rvm/scripts/rvm"
rvm use 2.0.0
# The first goal of this script is to perform a standard reference guided transcriptome assembly and comparison to the ref. annotation.
# Next, the new assembly will be used to guide the detection of strong signal near the transcription start sites with MACS2. This is accomplished as follows:
# First, the transcripts.gtf file is parsed into a dictionary.
# Next, the information is used to extract the reads found near the transcription start site in a strand specific manner with SAMtools view
# These reads are then passed to MACS for its model to evaluate the coverage near the transcription start site.
# This is mostly to just check the transcriptome assembly, showing a significant increase in signal near the transcription start site, indicative of multiple
# unique reads mapping near the start site.


CORES=6
INDIR=SAM_processed
LOGDIR=logs
CUFFQUANT=Cuffquant
REFERENCE=reference/merged.gtf
MASK=reference/mask.gtf
export REFERENCE
REFFASTA=reference/CAC.txt
EXPRDIR=counts
BAM=`/usr/bin/ls $INDIR/*.3.bam`
NS30='$CUFFQUANT/NS30A/abundances.cxb,$CUFFQUANT/NS30B/abundances.cxb'
NS75='$CUFFQUANT/NS75A/abundances.cxb,$CUFFQUANT/NS75B/abundances.cxb'
NS270='$CUFFQUANT/NS270A/abundances.cxb'

PICARD=/usr/local/picard-tools-1.67
TMPD=/home/mrals/ETP/tmp/

parallel -j$CORES 'htseq-count -f bam -r pos -s yes {} $REFERENCE > counts/{/.}.counts' ::: $BAM

for file in $BAM
do
    f=${file##*/};f=${f%*.*.*}
    #mkdir $CUFFQUANT/$f
    cuffquant -o $CUFFQUANT/$f -v -p $CORES -M $MASK -b $REFFASTA -u --library-type fr-firststrand $REFERENCE $file
done


cuffdiff -p $CORES -o $CUFFQUANT -L NS30,NS75,NS270 -T -b $REFFASTA -u -M $MASK --library-type fr-firststrand --library-norm-method geometric --min-reps-for-js-test 2 $REFERENCE $NS30 $NS75 $NS270

cuffnorm -p $CORES-o $CUFFQUANT -L NS30,NS75,NS270 --library-type fr-firststrand --library-norm-method geometric $REFERENCE $NS30 $NS75 $NS270

## EOF-------------------------------------------
