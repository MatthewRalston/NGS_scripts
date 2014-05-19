#NGS_scripts-paired

##Description
This project contains a number of BASH scripts that are used to process data
from next generation sequencing experiments. These scripts perform general
QC, trimming/clipping, mapping, and processing operations on a set of fastq files.



##Contents
* alignmentsummary.rb
* assemblyfilter.rb
* bedscript.rb
* bowtie.sh
* circos/circos.conf
* countsummary.rb
* coverage_calc.rb
* cufflinks.sh
* expression.sh
* initialqc.sh
* mapqsummary.rb
* postprocess.sh
* summary/summary.r
* trimfiltercheck.sh
* viewing.sh



##Workflow
1. initialqc.sh
   * Produces an initial quality check with fastqc
2. trimfiltercheck.sh
   * Trims/clips reads and additional qc.
3. bowtie.sh
   * Aligns reads to the genome after in-silico rRNA removal
4. alignmentsummary.rb
   * Produces summary tables for the rRNA removal and alignment
5. postprocess.sh
   * Processes and verifies alignment file integrity for downstream applications
6. viewing.sh
   * Coverage vectors, circos plots, and separates the alignments by strand
7. bedscript.rb
   * Processes coverage vectors to see 'fragment' coverage (vs. read)
8. mapqsummary.rb
   * Summarizes the (filtered) mapq qualities
9. cufflinks.sh
   * Reference based transcriptome assembly
10. expression.sh
   * Calculates fpkm/geometric read counts and raw read count gene expression measurements
11. coverage_calc.rb
   * Custom script: per-gene coverage vectors over the percentage of transcript length
12. summary.r
   * Generate summary and diagnostic plots from the sequencing experiments.


