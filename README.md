NGS_scripts-paired
===========

This project contains a number of BASH scripts that are used to process data
from next generation sequencing experiments. These scripts perform general
QC, trimming/clipping, mapping, and processing operations on a set of fastq files.


First, raw reads are processed through each step, assuming the all the appropriately formated files are available (gtf, etc.)
Next, after the cufflinks.sh script produces the merged.gtf file in Cufflinks_assemblies directory, the columns must be rearranged
so that the cufflinks-produced XLOC ids are replaced with the traditional gene ids. This can be done in excel or other software that 
accepts and produces tab-delimited files. Afterwards, this modified merged.gtf file should be stored in the reference directory.
Next, the gene counts and FPKM counts can be produced from the modified gtf file.

