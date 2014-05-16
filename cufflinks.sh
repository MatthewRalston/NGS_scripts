#!/bin/bash
#PBS -N assembly
#PBS -r n
#PBS -V
#PBS -l nodes=1:ppn=6
#PBS -l walltime=240:00:00
#PBS -d /home/mrals/ETP
#------------------------------------------------
# Title: cufflinks.sh
#
# Matt Ralston
#
# This script performs a reference transcriptome
# assembly for each condition, which are then
# merged into a consensus reference assembly for
# manual curation.
#
#------------------------------------------------


#set -e

#------------------------------------------------
# Parameters
#------------------------------------------------
# cufflinks: this is the location of the cufflinks
# executable
cufflinks='/usr/local/cufflinks-2.2.0/cufflinks'
# cuffcompare: this is the location of the cuffcompare
# executable
cuffcompare='/usr/local/cufflinks-2.2.0/cuffcompare'
# cuffmerge: this is the location of the cuffmerge
# executable
cuffmerge='/usr/local/cufflinks-2.2.0/cuffmerge'
# CORES: This is the number of cores available for
# parallelization
CORES=6
# INDIR: This is the location of the processed
# alignment files for transcriptome assembly.
INDIR=SAM_processed
# CUFFLINKS: This is the location where the 
# assemblies will be generated
CUFFLINKS=Cufflinks_assemblies
# CUFFCOMPARE: This is the location where the 
# output files from the assembly comparison will
# be produced.
CUFFCOMPARE=Cuffcompare
# MASK: This is the location of a gtf annotation
# of tRNAs and rRNAs from the C. acetobutylicum
# genome that will be masked from assembly.
MASK=reference/mask.gtf
# LOGDIR: This is the location where log files will
# be produced
LOGDIR=logs
# REFERENCE: This is the location of the reference
# CDS annotation to guide the transcriptome assembly
REFERENCE=reference/CAC.gtf
# REFFASTA: This is the location of the reference
# fasta file to aid transcriptome assembly.
REFFASTA=reference/CAC.txt
FILES=`/usr/bin/ls $INDIR/*.3.bam`

#########################     2     ##########################
# Next, the files are iterated through, performing a reference transcriptome assembly for each one.
for f in $FILES
do
        f=${f##*/};f=${f%*.*.*}
	mkdir $CUFFLINKS/$f $CUFFCOMPARE/$f
	$cufflinks -o $CUFFLINKS/${f} -p $CORES -b $REFFASTA -u -M $MASK -m 520 -s 49 -L $f -I 1 --min-intron-length 0 --max-multiread-fraction 0.2 --overlap-radius 5 --3-overhang-tolerance 10 --library-type fr-firststrand -g $REFERENCE --min-frags-per-transfrag 15 --3-overhang-tolerance 25 $INDIR/$f.3.bam
done

#########################     3     ##########################
# Here we move the location of all of the transcriptome assemblies into a file 'assemblies.txt'
rm assemblies.txt
find $CUFFLINKS/ -name 'transcripts.gtf' > assemblies.txt
#########################     4     ##########################
# Here, cuffcompare is used to find transcripts that do not contain the reference #BROKEN
$cuffcompare -r $REFERENCE -R -i assemblies.txt -s $REFFASTA -e 50 -d 10 -p Cuffcompare/novel -T -o $CUFFCOMPARE/novel
#########################     5     ##########################
# Here, cuffcompare is used to find transcripts that contain the reference #BROKEN
$cuffcompare -r $REFERENCE -Q -i assemblies.txt -s $REFFASTA -e 50 -d 10 -p Cuffcompare/canonical -T -o $CUFFCOMPARE/canonical
#########################     6     ##########################
# Finally, cuffmerge is used to merge all transcripts from the individual assemblies.
$cuffmerge -o $CUFFLINKS -g $REFERENCE -p $CORES -s $REFFASTA assemblies.txt







## EOF-------------------------------------------
