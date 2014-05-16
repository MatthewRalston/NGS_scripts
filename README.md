#NGS_scripts-paired

##Description
This project contains a number of BASH scripts that are used to process data
from next generation sequencing experiments. These scripts perform general
QC, trimming/clipping, mapping, and processing operations on a set of fastq files.



##Contents
* alignmentsummary.rb
* bedscript.rb
* bowtie.sh
* circos/circos.conf
* countsummary.rb
* cufflinks.sh
* expression.sh
* initialqc.sh
* mapqsummary.rb
* postprocess.sh
* summary.r
* trimfiltercheck.sh
* viewing.sh



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
7. mapqsummary.rb
   * This script summarizes the (filtered) mapq qualities for each condition
8. cufflinks.sh
   * This script performs a reference based transcriptome assembly
9. expression.sh
   * This script calculates fpkm/geometric read counts and raw read count gene expression measurements



