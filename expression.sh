#!/bin/bash
#PBS -N gene_express
#PBS -r n
#PBS -V
#PBS -l nodes=1:ppn=6
#PBS -l walltime=24:00:00
#PBS -d /home/mrals/NichSandoval
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
# REFERENCE: This is the reference annotation that 
# will be used to calculate gene expression.
export REFERENCE=reference/combined.gtf
export ECOREF=reference/eco_k12_mg1655.new.gtf
export LPLREF=reference/Lpl.new.gtf

# REFFASTA: This is a fasta format file
# of the reference Cac genome.
REFFASTA=reference/Lpl_eco_fos
# REFGENOME: This is a .genome file, deprecated.
REFGENOME=reference/CAC.genome
# EXPRDIR: This is the location where expression will be
# generated.
#EXPRDIR=counts
export EXPRDIR=Expression/counts/All
export ECODIR=Expression/counts/Eco
export LPLDIR=Expression/counts/Lpl
# BAM: This is a list of the files for expression measurments
BAM=`/usr/bin/ls $INDIR/*.3.bam`
# PICARD: This is the location of the picard tool jars
export PICARD=/usr/local/picard-tools-1.67
# TMPD: This is a location for temporary files.
TMPD=/home/mrals/NichSandoval/tmp/
export TMPD





#########################     HTSeq combined

#parallel -j$CORES 'samtools view -h {} | ./unprocess.rb > $TMPD/{/.}.tempbam' ::: $BAM
#TEMPBAM=`/usr/bin/ls $TMPD/*.tempbam`
#parallel -j$CORES 'java -jar $PICARD/MarkDuplicates.jar INPUT={} OUTPUT=/dev/stdout METRICS_FILE=/dev/null REMOVE_DUPLICATES=true | samtools sort -n - $TMPD/{/.}' ::: $TEMPBAM
#rm $TMPD/*.tempbam
BAM=`/usr/bin/ls $TMPD/*.3.bam`
# paired
parallel -j$CORES 'samtools view -h -f 1 {} | htseq-count -f sam -t CDS -a 20 -r name -m intersection-nonempty -s yes - $REFERENCE > $EXPRDIR/{/.}.counts.paired' ::: $BAM
# unpaired
parallel -j$CORES 'samtools view -h -F 1 {} | htseq-count -f sam -t CDS -a 20 -r name -m intersection-nonempty -s yes - $REFERENCE > $EXPRDIR/{/.}.counts.unpaired' ::: $BAM

echo $EXPRDIR | ./countsummary.rb 


#########################     HTSeq Eco


# paired
parallel -j$CORES 'samtools view -h -f 1 {} | htseq-count -f sam -t CDS -a 20 -r name -m intersection-nonempty -s yes - $ECOREF > $ECODIR/{/.}.counts.paired' ::: $BAM
# unpaired
parallel -j$CORES 'samtools view -h -F 1 {} | htseq-count -f sam -t CDS -a 20 -r name -m intersection-nonempty -s yes - $ECOREF > $ECODIR/{/.}.counts.unpaired' ::: $BAM

echo $ECODIR | ./countsummary.rb


#########################     HTSeq Lpl

# paired
parallel -j$CORES 'samtools view -h -f 1 {} | htseq-count -f sam -t CDS -a 20 -r name -m intersection-nonempty -s yes - $LPLREF > $LPLDIR/{/.}.counts.paired' ::: $BAM
# unpaired
parallel -j$CORES 'samtools view -h -F 1 {} | htseq-count -f sam -t CDS -a 20 -r name -m intersection-nonempty -s yes - $LPLREF > $LPLDIR/{/.}.counts.unpaired' ::: $BAM
#rm $TMPD/*.3.bam

echo $LPLDIR | ./countsummary.rb




## EOF-------------------------------------------
