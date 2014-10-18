=begin
------------------------------------------------------------------------------
--                                                                          --
--                                MATTHEW RALSTON                           --
--                                                                          --
--                               S K E L E T O N . R B                      --
--                                                                          --
------------------------------------------------------------------------------
-- TITLE                                                                    --
--                                                                          --
--  Summer 2013                                                             --
--                                                                          --
------------------------------------------------------------------------------
-- This file is designed to identify discrepancies between the old CDS      --
-- annotation and the new transcript assembly.                              --
-- records are first divided into 4 categories: transcript, CDS, rrna, trna.--
-- Each transcript is tested to see if it eclipses/includes at least one CDS
-- rrna or trna. This status is included in a description. 
-- Next, partial overlaps are checked for, and these transcripts are marked --
-- for followup.
-- Then, each CDS is given a status: partial overlap, full inclusion, none  --
-- and it's parent transcript is noted.
-- Finally, each CDS and with a parent and its associated parent are output --
-- into the gene<tab>transcript style file.
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
