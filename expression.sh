#!/bin/bash
#PBS -N gene_express
#PBS -r n
#PBS -V
#PBS -l nodes=1:ppn=4
#PBS -l walltime=240:00:00
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
REFERENCE=reference/combined.gtf
ECOREF=reference/eco_k12_mg1655.gtf
LPLREF=reference/Lpl.gtf
export REFERENCE
# REFFASTA: This is a fasta format file
# of the reference Cac genome.
REFFASTA=reference/Lpl_eco
# REFGENOME: This is a .genome file, deprecated.
REFGENOME=reference/CAC.genome
# EXPRDIR: This is the location where expression will be
# generated.
#EXPRDIR=counts
export EXPRDIR=Expression/counts
export ECODIR=Expression/Eco
export LPLDIR=Expression/Lpl
# BAM: This is a list of the files for expression measurments
BAM=`/usr/bin/ls $INDIR/*.3.bam`

# TMPD: This is a location for temporary files.
TMPD=/home/mrals/Final/tmp/
export TMPD




#########################     HTSeq combined

parallel -j$CORES 'java -jar $PICARD/MarkDuplicates.jar INPUT={} OUTPUT=/dev/stdout METRICS_FILE=/dev/null REMOVE_DUPLICATES=true | java -jar $PICARD/SortSam.jar INPUT=/dev/stdin OUTPUT=$TMPD/{/} SORT_ORDER=queryname' ::: $BAM
BAM=`/usr/bin/ls $TMPD/*.3.bam`
# paired
parallel -j$CORES 'samtools view -h -f 1 {} | htseq-count -f sam -a 20 -r name -m intersection-nonempty -s yes - $REFERENCE > $EXPRDIR/{/.}.counts.paired' ::: $BAM
# unpaired
parallel -j$CORES 'samtools view -h -F 1 {} | htseq-count -f sam -a 20 -r name -m intersection-nonempty -s yes - $REFERENCE > $EXPRDIR/{/.}.counts.unpaired' ::: $BAM

echo /home/mrals/NichSandoval/Expression/counts | ./countsummary.rb 


#########################     HTSeq Eco


# paired
parallel -j$CORES 'samtools view -h -f 1 {} | htseq-count -f sam -a 20 -r name -m intersection-nonempty -s yes - $ECOREF > $ECODIR/{/.}.counts.paired' ::: $BAM
# unpaired
parallel -j$CORES 'samtools view -h -F 1 {} | htseq-count -f sam -a 20 -r name -m intersection-nonempty -s yes - $ECOREF > $ECODIR/{/.}.counts.unpaired' ::: $BAM

echo /home/mrals/NichSandoval/Expression/Eco | ./countsummary.rb


#########################     HTSeq Lpl

# paired
parallel -j$CORES 'samtools view -h -f 1 {} | htseq-count -f sam -a 20 -r name -m intersection-nonempty -s yes - $LPLREF > $LPLDIR/{/.}.counts.paired' ::: $BAM
# unpaired
parallel -j$CORES 'samtools view -h -F 1 {} | htseq-count -f sam -a 20 -r name -m intersection-nonempty -s yes - $LPLREF > $LPLDIR/{/.}.counts.unpaired' ::: $BAM
rm $TMPD/*.3.bam

echo /home/mrals/NichSandoval/Expression/Lpl | ./countsummary.rb




## EOF-------------------------------------------
