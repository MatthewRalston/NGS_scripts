#NGS_scripts: paired-end

##Description
This project contains a number of BASH and ruby scripts that are used to process data
from RNA sequencing experiments. These scripts include [Torque](http://www.adaptivecomputing.com/products/open-source/torque) headers for use on a cluster.
These scripts perform general QC, trimming/clipping, mapping, processing, visualization, quantification, coverage-summary
and additional summary operations on a set of fastq/BAM files. This group of files is designed
to be run in the sequence described below. These scripts assume that the software is installed
and accessible in the PATH environment variable.  

##Dependencies
* [Ruby](https://www.ruby-lang.org/en/)
* [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
* [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* [SAMtools](http://samtools.sourceforge.net/)
* [Picard](http://picard.sourceforge.net/)
* [Circos](http://circos.ca/)
* [Bedtools](http://bedtools.readthedocs.org/en/latest/)
* [Cufflinks suite](http://cufflinks.cbcb.umd.edu/)
* [HT-Seq](http://www-huber.embl.de/users/anders/HTSeq/doc/count.html)
* [R](http://www.r-project.org/)
* [Cairo Graphics](http://cran.r-project.org/web/packages/cairoDevice/)




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


