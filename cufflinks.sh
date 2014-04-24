#!/bin/bash
#PBS -N assembly
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
CUFFLINKS=Cufflinks_assemblies
CUFFCOMPARE=Cuffcompare
#PEAKCALLING=Peak_call
LOGDIR=logs
which cufflinks >> temp.txt
REFERENCE=reference/CAC.gtf
REFFASTA=reference/CAC.txt
EXPRDIR=counts
FILES=`/usr/bin/ls $INDIR/*.3.bam`

PICARD=/usr/local/picard-tools-1.67
TMPD=/home/mrals/ETP/tmp/



for f in $FILES
do
        f=${f##*/};f=${f%*.*.*}
	mkdir $CUFFLINKS/$f $CUFFCOMPARE/$f

	# First, this script uses cufflinks to produce a reference guided transcriptome assembly
	# Cufflinks produces a transcripts.gtf which contains the new assembly
	# It produces a isoforms.fpkm_tracking, a estimated isoform expression values in FPKM Tracking Format
	# Finally, it produces a genes.fpkm_tracking, a etimated gene expression values in FPKM Tracking Format
	cufflinks -o $CUFFLINKS/${f} -p $CORES --library-type fr-firststrand -g $REFERENCE --3-overhang-tolerance 25 $INDIR/$f.3.bam
	# Second, it runs cuffcompare to compare the new assembly to the reference assembly.
	# This produces a .stats file that describes the sensitivity and specificity for detecting nucleotide, exons, introns, transcripts genes
	# It also produces a .combined.gtf which is mostly pointless here.
	# It produces a .tracking file which keeps track of the transcripts, and labels transfrags as novel or not based on the reference annotation
	# It produces a .refmap file which also labels the transcripts as novel or not
	# It produces a .tmap file which lists retains the novel/partial/full transcript match (to the refrence)
	#    and lists the closest matching reference transcript, the coverage, expression (FPKM), confidence intervals for the FPKM, and length of the resulting transcript
	#cuffcompare -r $REFERENCE -o $CUFFCOMPARE/${f} -e 10 -T $CUFFLINKS/${f}/transcripts.gtf

done

find $CUFFLINKS/ -name 'transcripts.gtf' > assemblies.txt
cuffmerge -o $CUFFLINKS -g $REFERENCE -p $CORES -s $REFFASTA assemblies.txt
gtfToGenePred -genePredExt $CUFFLINKS/merged.gtf $CUFFLINKS/tmp.refflat
awk 'BEGIN{FS="\t"};{print $12"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' $CUFFLINKS/tmp.refflat > $CUFFLINKS/merged.refflat
REFFLAT=$CUFFLINKS/merged.refflat
ANNOTATION=/home/mrals/ETP/$CUFFLINKS/merged.gtf
export ANNOTATION

for file in $FILES
do
    ############ OPTIONAL #############
	# FIX ME: CollectRnaSeqMetrics: Program that collects information about alignment of reads to coding, intronic, UTR, intergenic, ribosomal, etc.
	java -jar $PICARD/CollectRnaSeqMetrics.jar INPUT=$file REF_FLAT=$REFFLAT TMP_DIR=$TMPD STRAND_SPECIFICITY=FIRST_READ_TRANSCRIPTION_STRAND MINIMUM_LENGTH=50 CHART_OUTPUT=$LOGDIR/${file##*/*}.RNAseq_metrics.pdf OUTPUT=$LOGDIR/${file##*/*}.RNAseq_metrics.log REFERENCE_SEQUENCE=$REFFASTA ASSUME_SORTED=true
done







## EOF-------------------------------------------
