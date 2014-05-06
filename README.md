#NGS_scripts-paired

##Description
This project contains a number of BASH scripts that are used to process data
from next generation sequencing experiments. These scripts perform general
QC, trimming/clipping, mapping, and processing operations on a set of fastq files.



##Contents
* initialqc.sh
* trimfiltercheck.sh
* bowtie.sh
* alignmentsummary.rb
* postprocess.sh
* viewing.sh
* cufflinks.sh
* expression.sh
* countsummary.rb
* bedscript.rb



##Workflow
1. initialqc.sh
   * This script produces an initial quality check with fastqc
2. trimfiltercheck.sh
   * This script trims reads with trimmomatic and then performs a repeated quality check
3. bowtie.sh
   * This script aligns reads to the genome after in-silico rRNA removal
4. alignmentsummary.rb
   * This script produces summary tables for the rRNA removal and alignment
5. postprocess.sh
   * This script processes and checks the integrity of the alignment files
6. viewing.sh
   * This script produces coverage vectors, circos plots, and separates the alignments by strand
7. cufflinks.sh
   * This script performs a reference based transcriptome assembly
8. expression.sh
   * This script calculates fpkm/geometric read counts and raw read count gene expression measurements



