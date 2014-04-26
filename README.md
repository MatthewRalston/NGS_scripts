NGS_scripts-paired
===========

This project contains a number of BASH scripts that are used to process data
from next generation sequencing experiments. These scripts perform general
QC, trimming/clipping, mapping, and processing operations on a set of fastq files.


1. Paired end fastq are checked with initialqc.sh
2. These files are then processed by trimfiltercheck.sh, which clips adapters, trims by quality, and filters the reads.
3. Next, reads are aligned to rRNA sequences (in silico rRNA removal) before aligning the difference to the genome in bowtie.sh.
4. The alignment files are processed (sorting, marking duplicates, Picard metrics, etc.) with postprocess.sh
5. A consensus transcriptome assembly is created from the files using cufflinks.sh.
6. The assembly must be manually curated with the assistance of a genome browser with cufflinks.sh.
7. The expression levels are tabulated through several means, most of which are dependent on a good transcriptome assembly, in expression.sh.



