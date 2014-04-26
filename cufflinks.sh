#!/bin/bash
#PBS -N assembly
#PBS -r n
#PBS -V
#PBS -l nodes=1:ppn=6
#PBS -l walltime=240:00:00
#PBS -d /home/mrals/ETP

#set -e


#General
PATH=$PATH:/home/mrals/pckges/Vienna/bin:/home/mrals/home/bin/:/home/mrals/bin/
#NGS
PATH=$PATH:/home/mrals/pckges/sratoolkit.2.3.4-2-centos_linux64/bin/:/usr/local/bowtie-0.12.7/:/usr/local/bowtie2-2.1.0/:/usr/local/bwa-0.7.4/:/usr/local/cufflinks-2.2.0/:/usr/local/FastQC/:/usr/local/samtools-0.1.18/:/usr/local/tophat-2.0.4/:/home/mrals/pckges/seqtk-master/

export PATH
[[ -s "$HOME/.rvm/scripts/rvm" ]] && source "$HOME/.rvm/scripts/rvm"
rvm use 2.0.0

#########################     1     ##########################
# First, several environmental variables are defined, including the name of the executables that will be called,
# The number of processors for parallelization, the directory of input files, the directory of output files,
# a .gtf file of sequences to omit from assembly (e.g. tRNA, rRNA), the directory of log files to describe the assembly,
# and the reference genome and annotation to guide the assembly.

cufflinks='/usr/local/cufflinks-2.2.0/cufflinks'
cuffcompare='/usr/local/cufflinks-2.2.0/cuffcompare'
cuffmerge='/usr/local/cufflinks-2.2.0/cuffmerge'
CORES=6
INDIR=SAM_processed
CUFFLINKS=Cufflinks_assemblies
CUFFCOMPARE=Cuffcompare
MASK=reference/mask.gtf
LOGDIR=logs
REFERENCE=reference/CAC.gtf
REFFASTA=reference/CAC.txt
FILES=`/usr/bin/ls $INDIR/*.3.bam`

#########################     2     ##########################
# Next, the files are iterated through, performing a reference transcriptome assembly for each one.
for f in $FILES
do
        f=${f##*/};f=${f%*.*.*}
	mkdir $CUFFLINKS/$f $CUFFCOMPARE/$f
	$cufflinks -o $CUFFLINKS/${f} -p $CORES -b $REFFASTA -u -M $MASK -m 520 -s 49 -L $f -I 1 --min-intron-length 0 --max-multiread-fraction 0.2 --overlap-radius 5 --3-overhang-tolerance 10 --library-type fr-firststrand -g $REFERENCE --3-overhang-tolerance 25 $INDIR/$f.3.bam
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
