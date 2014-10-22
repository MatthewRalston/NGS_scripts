#!/bin/bash
#PBS -N gene_express
#PBS -r n
#PBS -V
#PBS -l nodes=1:ppn=6
#PBS -l walltime=12:00:00
#PBS -d /home/mrals/Final
#------------------------------------------------
# Title: expression.sh
#
# Copyright 2014 Matt Ralston
# 
# This script first produces raw read counts for
# each gene with HTSeq. Next, Cuffquant and
# Cuffnorm produce normalized gene expression
# measurements. Finally, cuffdiff can be used for
# differential expression measurements. Then,
# Picard RNAseqmetrics are calculated for an
# improved gene assembly.
# 
#------------------------------------------------


#set -e

#------------------------------------------------
# Parameters
#------------------------------------------------
# cuffdiff: This is the location of the cuffdiff
# executable.
cuffdiff='/usr/local/cufflinks-2.2.0/cuffdiff'
# CORES: This is the number of processors available
# for parallelization.
CORES=6
# INDIR: This is the location of the processed
# alignment files that will be input for expression
# measurements.
INDIR=SAM_processed
# LOGDIR: This is the location of the directory
# where logs from the picard RNAseq metrics will be
# generated
LOGDIR=logs
# CUFFQUANT: This is the location where Cuffquant/
# Cuffnorm output will be produced.
CUFFQUANT=Cuffquant
# REFERENCE: This is the reference annotation that 
# will be used to calculate gene expression.
#REFERENCE=reference/ref.gtf
#REFERENCE=reference/CAC.gtf
REFERENCE=annotation/Trin-blast.gtf
#REFERENCE=summary/summary2000.gtf
export REFERENCE
# MASK: This is an annotation file of genes that will
# be omitted from normalization procedures (tRNA, rRNA)
MASK=reference/mask.gtf
# REFFASTA: This is a fasta format file
# of the reference Cac genome.
REFFASTA=reference/CAC.txt
# REFGENOME: This is a .genome file, deprecated.
REFGENOME=reference/CAC.genome
# EXPRDIR: This is the location where expression will be
# generated.
#EXPRDIR=counts
export EXPRDIR=Expression/rawcounts
# BAM: This is a list of the files for expression measurments
BAM=`/usr/bin/ls $INDIR/*.3.bam`
# NS30: This is a comma separated list of cxb cuffquant output
# files for replicates of this condition
NS30='Cuffquant/NS30A/abundances.cxb,Cuffquant/NS30B/abundances.cxb'
# NS75: This is a comma separated list of cxb cuffquant output
# files for replicates of this condition
NS75='Cuffquant/NS75A/abundances.cxb,Cuffquant/NS75B/abundances.cxb'
# NS270: This is a comma separated list of cxb cuffquant output
# files for replicates of this condition
NS270='Cuffquant/NS270A/abundances.cxb'
# PICARD: This is the location of picard jar files, used for 
# RNAseq metrics.

PICARD=/usr/local/picard-tools-1.67
export PICARD
# TYPE: This is the type of feature to count for HTSeqcount
export TYPE="transcript"
# TMPD: This is a location for temporary files.
TMPD=/home/mrals/Final/tmp/
export TMPD




#########################     HTSeq

#parallel -j$CORES 'samtools view -h {} | ./unprocess.rb | samtools view -bhS - > $TMPD/{/.}.tempbam' ::: $BAM
#TEMPBAM=`/usr/bin/ls $TMPD/*.tempbam`
#parallel -j$CORES 'java -jar $PICARD/MarkDuplicates.jar INPUT={} OUTPUT=/dev/stdout METRICS_FILE=/dev/null REMOVE_DUPLICATES=true | samtools sort -n - $TMPD/{/.}' ::: $TEMPBAM
#rm $TMPD/*.tempbam
#BAM=`/usr/bin/ls $TMPD/*.3.bam`
# paired
parallel -j$CORES 'samtools view -h -f 1 {} | htseq-count -f sam -t $TYPE -a 20 -r name -m intersection-nonempty -s yes - $REFERENCE > $EXPRDIR/{/.}.counts.paired' ::: $BAM
# unpaired
parallel -j$CORES 'samtools view -h -F 1 {} | htseq-count -f sam -t $TYPE -a 20 -r name -m intersection-nonempty -s yes - $REFERENCE > $EXPRDIR/{/.}.counts.unpaired' ::: $BAM
#rm $TMPD/*.3.bam

#./countsummary.rb


#########################     CUFFQUANT

#for file in $BAM
#do
#    f=${file##*/};f=${f%*.*.*}
    #mkdir $CIRC/$f
    #mkdir $CUFFQUANT/$f
    #cuffquant -o $CUFFQUANT/$f -v -p $CORES -M $MASK -b $REFFASTA -u --library-type fr-firststrand $REFERENCE $file
#done

#########################     CUFFNORM
#cuffnorm -p $CORES -o $CUFFQUANT -L NS30,NS75,NS270 --library-type fr-firststrand --library-norm-method geometric $REFERENCE $NS30 $NS75 $NS270
#mv $CUFFQUANT/genes.count_table $CUFFQUANT/genes.count_table.geometric
#cuffnorm -p $CORES -o $CUFFQUANT -L NS30,NS75,NS270 --library-type fr-firststrand --library-norm-method classic-fpkm $REFERENCE $NS30 $NS75 $NS270
#mv $CUFFQUANT/genes.count_table $CUFFQUANT/genes.count_table.fpkm




#########################     CUFFDIFF
#/usr/local/cufflinks-2.2.0/cuffdiff -p $CORES -o $CUFFQUANT -T -b $REFFASTA -u -M $MASK --library-type fr-firststrand --library-norm-method geometric --min-reps-for-js-test 1 $REFERENCE $NS30 $NS75 $NS270




# Additional RNAseq metrics with picard
#gtfToGenePred -genePredExt $CUFFLINKS/merged.gtf $CUFFLINKS/tmp.refflat
#awk 'BEGIN{FS="\t"};{print $12"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' $CUFFLINKS/tmp.refflat > $CUFFLINKS/merged.refflat
#REFFLAT=$CUFFLINKS/merged.refflat
#ANNOTATION=/home/mrals/ETP/$CUFFLINKS/merged.gtf
#export ANNOTATION
#for file in $FILES
#do
#	java -jar $PICARD/CollectRnaSeqMetrics.jar INPUT=$file REF_FLAT=$REFFLAT TMP_DIR=$TMPD STRAND_SPECIFICITY=FIRST_READ_TRANSCRIPTION_STRAND MINIMUM_LENGTH=50 CHART_OUTPUT=$LOGDIR/${file##*/*}.RNAseq_metrics.pdf OUTPUT=$LOGDIR/${file##*/*}.RNAseq_metrics.log REFERENCE_SEQUENCE=$REFFASTA ASSUME_SORTED=true
#done

#########################     Expression ratios of counts
# Here the counts are calculated as expression ratios
# log2(ratio) = log2 ( gene-a-t1 / gene-a-t2 )
# INDIR = Expression/counts
# OUTDIR = circos/data/DE_BuOH-30_0.txt


## EOF-------------------------------------------
