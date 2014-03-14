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
# The first goal of this script is to perform a standard reference guided transcriptome assembly and comparison to the ref. annotation.
# Next, the new assembly will be used to guide the detection of strong signal near the transcription start sites with MACS2. This is accomplished as follows:
# First, the transcripts.gtf file is parsed into a dictionary.
# Next, the information is used to extract the reads found near the transcription start site in a strand specific manner with SAMtools view
# These reads are then passed to MACS for its model to evaluate the coverage near the transcription start site.
# This is mostly to just check the transcriptome assembly, showing a significant increase in signal near the transcription start site, indicative of multiple
# unique reads mapping near the start site.
INDIR=SAM_processed
CUFFLINKS=Cufflinks_assemblies
CUFFCOMPARE=Cuffcompare
PEAKCALLING=Peak_call
REFERENCE=reference/CAC.gff
FILES=`/usr/bin/ls $INDIR/*.3.bam`
for f in $FILES
do
	mkdir $CUFFLINKS/${f%.*.*} $CUFFCOMPARE/${f%.*.*} $PEAKCALLING/${f%.*.*}
	# First, this script uses cufflinks to produce a reference guided transcriptome assembly
	# Cufflinks produces a transcripts.gtf which contains the new assembly
	# It produces a isoforms.fpkm_tracking, a estimated isoform expression values in FPKM Tracking Format
	# Finally, it produces a genes.fpkm_tracking, a etimated gene expression values in FPKM Tracking Format
	cufflinks -o $CUFFLINKS/${f%.*.*} -p 2 -g $REFERENCE --3-overhang-tolerance 50 $f
	# Second, it runs cuffcompare to compare the new assembly to the reference assembly.
	# This produces a .stats file that describes the sensitivity and specificity for detecting nucleotide, exons, introns, transcripts genes
	# It also produces a .combined.gtf which is mostly pointless here.
	# It produces a .tracking file which keeps track of the transcripts, and labels transfrags as novel or not based on the reference annotation
	# It produces a .refmap file which also labels the transcripts as novel or not
	# It produces a .tmap file which lists retains the novel/partial/full transcript match (to the refrence)
	#    and lists the closest matching reference transcript, the coverage, expression (FPKM), confidence intervals for the FPKM, and length of the resulting transcript
	cuffcompare -r $REFERENCE -o $CUFFCOMPARE/${f%.*.*} $CUFFLINKS/${f%.*.*}/transcripts.gtf
	for 
done




## EOF-------------------------------------------
