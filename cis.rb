=begin
------------------------------------------------------------------------------
--                                                                          --
--                                MATTHEW RALSTON                           --
--                                                                          --
--                                   C I S . R B                            --
--                                                                          --
------------------------------------------------------------------------------
-- TITLE                                                                    --
--                                                                          --
--  Summer 2013                                                             --
--                                                                          --
------------------------------------------------------------------------------
-- This file is designed to acquire reads belonging to potential antisense  --
-- transcribed regions.
-- First a gtf file is loaded into memory. Then, a reversed gtf file with   --
-- the strand switched for each gene is produced.
-- Second, the reverse gtf is used with bedtools intersect to find all BAM  --
-- records that overlap the opposite strand on any gene.
-- Third, this BAM file is then filtered for the set difference with the    --
-- traditional CDS annotation. This would produce a BAM file that surely    --
-- contains the antisense transcripts. 
-- Fourth, this BAM file can then be used in a set difference with the      --
-- original combined.gtf file. Any BAM records that are found in this way   --
-- could be 'assembled' and quantified. These BAM records would be on the   --
-- antisense strand, not belong to a CDS or to a gene. If no records are    --
-- found in this way, finding antisense transcription may be difficult.     --
-- However, we see evidence for such reads by glancing at the alignment     --
-- results. It is clear that such fragments/reads were not assembled, and   --
-- represent an unquantified and source of potentially novel information.   --
-- 
 
--                                                                          --
------------------------------------------------------------------------------
=end

################################################
#
#               R E Q U I R E
#
################################################





################################################
#
#               U S E R    V A R I A B L E S
#
################################################





################################################
#
#               S U B - R O U T I N E S
#
################################################

################################################

################################################





#*****************************************************************************#
################################################
#
#-----------------------------------------------
#
#                   M A I N
#-----------------------------------------------
#
################################################







##########################  E O F   ############################################
