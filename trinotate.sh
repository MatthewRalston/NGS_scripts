#!/bin/bash
#PBS -N trinotate
#PBS -r n
#PBS -V
#PBS -l nodes=1:ppn=6
#PBS -l walltime=400:00:00
#PBS -d /home/mrals/Final/Trinotate
#------------------------------------------------
# Title: initialqc.sh
#
# Copyright 2014 Matt Ralston
# 
# This script performs a quality check on fastq
# format files in a directory, reporting the 
# results to a separate directory.
# 
#------------------------------------------------

#------------------------------------------------
# Functions
#------------------------------------------------



#------------------------------------------------
# Parameters
#------------------------------------------------
# CORES: the number of processors available
CORES=6
# OUTDIR: where the output should be written
OUTDIR=Trinotate
# FASTA: assembly to annotate
FASTA=Trinity-GG.fasta
# REFPROTEOME: reference proteome
REFPROTEOME=../reference/CAC_proteins.fasta
# MINSIZE: minimum protein size
MINSIZE=80
# PFAM: location of pfam HMM database
PFAM=databases/Pfam-A.hmm
# UNIPROT: location of UNIPROT database
UNIPROT=databases/uniprot_sprot.fasta
# UNIREF: location of UNIREF90 database
UNIREF=databases/uniref90.fasta
# GENE-TRANS MAP: Tab delmited file of 'gene\ttranscript\n'
GENETRANSMAP=genetransmap.txt


#------------------------------------------------
# ORF finding
#------------------------------------------------
#TransDecoder -t $FASTA --train $REFPROTEOME -m $MINSIZE -S --search_pfam $PFAM --CPU $CORES


#------------------------------------------------
# Annotation
#------------------------------------------------
# RBB
#blastx -query $FASTA -db $UNIPROT -num_threads $CORES -max_target_seqs 1 -outfmt 6 > blastx.outfmt6
#blastp -query $FASTA.transdecoder.pep -db $UNIPROT -num_threads $CORES -max_target_seqs 1 -outfmt 6 > blastp.outfmt6
blastx -query $FASTA -db $UNIREF -max_target_seqs 1 -outfmt 6 > uniref90.blastx.outfmt6
blastp -query $FASTA.transdecoder.pep -db $UNIREF -max_target_seqs 1 -outfmt 6 > uniref90.blastp.outfmt6
# HMMer
#hmmscan --cpu $CORES --domtblout TrinotatePFAM.out $PFAM $FASTA.transdecoder.pep > pfam.log
# Signal-P
#signalp -f short -n signalp.out $FASTA.transdecoder.pep
# TM-HMM
#tmhmm --short < $FASTA.transdecoder.pep > tmhmm.out
# Trinotate
Trinotate Trinotate.sqlite init --gene_trans_map $GENETRANSMAP --transcript_fasta $FASTA --transdecoder_pep $FASTA.transdecoder.pep LOAD_swissprot_blastp blastp.outfmt6 LOAD_swissprot_blastx blastx.outfmt6 LOAD_trembl_blastp uniref90.blastp.outfmt6 LOAD_trembl_blastx uniref90.blastx.outfmt6 LOAD_pfam TrinotatePFAM.out LOAD_tmhmm tmhmm.out LOAD_signalp signalp.out
Trinotate Trinotate.sqlite report > trinotate_annotation_report.xls


## EOF-------------------------------------------
