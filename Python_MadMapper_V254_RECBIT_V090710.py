#!/usr/bin/python

#################################################################
#                                                               #
#                          MAD MAPPER                           #
#                                                               #
#                      MAD MAPPING PROGRAM                      #
#                                                               #
#                     PART 1  ( CLUSTERING )                    #
#                                                               #
#                              and                              #
#                                                               #
#                  PART 2  ( MAP CONSTRUCTION )                 #
#                                                               #
#                      COPYRIGHT  2004 2005                     #
#                        Alexander Kozik                        #
#                                                               #
#                      http://www.atgc.org/                     #
#                                                               #
#             UCD Genome Center. R.Michelmore group             #
#                                                               #
#################################################################

#################################################################
#                                   +-------+                   #
#                                   |  BIT  |                   #
#                  SCORING SYSTEM:  |       |                   #
#                                   |  REC  |                   #
#                                   +-------+                   #
#                                                               #
#    .      +-------+-------+-------+-------+-------+-------+   #
#      .    |       |       |       |       |       |       |   #
#        .  |   A   |   B   |   C   |   D   |   H   |   -   |   #
#          .|       |       |       |       |       |       |   #
#   +-------*-------+-------+-------+-------+-------+-------+   #
#   |       | . 6   |  -6   |  -4   |   4   |  -2   |   0   |   #
#   |   A   |       |       |       |       |       |       |   #
#   |       |   0  .|   1   |   1   |   0   |  0.5  |   0   |   #
#   +-------+-------*-------+-------+-------+-------+-------+   #
#   |       |  -6   | . 6   |   4   |  -4   |  -2   |   0   |   #
#   |   B   |       |       |       |       |       |       |   #
#   |       |   1   |   0  .|   0   |   1   |  0.5  |   0   |   #
#   +-------+-------+-------*-------+-------+-------+-------+   #
#   |       |  -4   |   4   | . 4   |  -4   |   0   |   0   |   #
#   |   C   |       |       |       |       |       |       |   #
#   |       |   1   |   0   |   0  .|   1   |   0   |   0   |   #
#   +-------+-------+-------+-------*-------+-------+-------+   #
#   |       |   4   |  -4   |  -4   | . 4   |   0   |   0   |   #
#   |   D   |       |       |       |       |       |       |   #
#   |       |   0   |   1   |   1   |   0  .|   0   |   0   |   #
#   +-------+-------+-------+-------+-------*-------+-------+   #
#   |       |  -2   |  -2   |   0   |   0   | . 2   |   0   |   #
#   |   H   |       |       |       |       |       |       |   #
#   |       |  0.5  |  0.5  |   0   |   0   |   0  .|   0   |   #
#   +-------+-------+-------+-------+-------+-------*-------+   #
#   |       |   0   |   0   |   0   |   0   |   0   | . 0   |   #
#   |   -   |       |       |       |       |       |       |   #
#   |       |   0   |   0   |   0   |   0   |   0   |   0  .|   #
#   +-------+-------+-------+-------+-------+-------+-------*.  #
#                                                               #
#                                                               #
#   NOTES:                                                      #
#      C - NOT A ( H or B )                                     #
#      D - NOT B ( H or A )                                     #
#      H - A and B                                              #
#                                                               #
#################################################################

#################################################################
#                                                               #
#                        EXAMPLES OF SCORING:                   #
#                                                               #
#                                                               #
#  POSITIVE LINKAGE:                                            #
#                                                               #
#  AAAAAAAAAAAAAAAAAAAA     BIT SCORE = 6*20 = 120              #
#  AAAAAAAAAAAAAAAAAAAA     REC SCORE = 0 (0.0)                 #
#                    ..                                         #
#  AAAAAAAAAAAAAAAAAAAA     BIT SCORE = 6*18 - 6*2 = 96         #
#  AAAAAAAAAAAAAAAAAABB     REC SCORE = 2 (2/20 = 0.1)          #
#                                                               #
#  AAAAAAAAAABBBBBBBBBB     BIT SCORE = 6*10 + 6*10 = 120       #
#  AAAAAAAAAABBBBBBBBBB     REC SCORE = 0 (0.0)                 #
#           ..                                                  #
#  AAAAAAAAABABBBBBBBBB     BIT SCORE = 6*18 - 6*2 = 96         #
#  AAAAAAAAAABBBBBBBBBB     REC SCORE = 2 (2/20 = 0.1)          #
#                                                               #
#                                                               #
#  NO LINKAGE:                                                  #
#            ..........                                         #
#  AAAAAAAAAAAAAAAAAAAA     BIT SCORE = 6*10 - 6*10 = 0         #
#  AAAAAAAAAABBBBBBBBBB     REC SCORE = 10 (10/20 = 0.5)        #
#   . . . . . . . . . .                                         #
#  BBBAABBAAAAAAABAABBB     BIT SCORE = 6*10 - 6*10 = 0         #
#  BABBAABBABABABBBAABA     REC SCORE = 10 (10/20 = 0.5)        #
#                                                               #
#                                                               #
#  NEGATIVE LINKAGE:                                            #
#    ..................                                         #
#  AAAAAAAAAAAAAAAAAAAA     BIT SCORE = 6*2 - 6*18 = -96        #
#  AABBBBBBBBBBBBBBBBBB     REC SCORE = 18 (18/20 = 0.9)        #
#    ..................                                         #
#  ABABABABABABABABABAB     BIT SCORE = 6*2 - 6*18 = -96        #
#  ABBABABABABABABABABA     REC SCORE = 18 (18/20 = 0.9)        #
#                                                               #
#                                                               #
#################################################################

##################################################################################
#1   5   10   15   20   25   30   35   40   45   50   55   60   65   70   75   80#
#-=3+5-=3+5-=3+5-=3+5-=3+5-=3+5-=3+5-=3+5-=3+5-=3+5-=3+5-=3+5-=3+5-=3+5-=3+5-=3+5#
##################################################################################

##################################################################################
#                           --+============+--                                   #
#                              README_FIRST:                                     #
#                            PYTHON MAD MAPPER                                   #
#                        HOW IT WORKS AND WHAT IT DOES:                          #
#                  RULES AND ALGORITHMS OF MADMAPPER APPROACH                    #
#             ---++==========================================++---               #
#                                                                                #
#   --------------------------------------------------------------------------   #
#  NOTE:                                                                         #
#         LABELS LIKE ### README_00X ### ARE INCORPORATED INTO THE SOURCE CODE   #
#         TO HELP TO FIND FUNCTIONS AND VARIABLES DESCRIBED IN THIS DOCUMENT     #
#   --------------------------------------------------------------------------   #
#                                                                                #
#  = ASSUMPTIONS AND RESTRICTIONS =                                              #
#  0. MadMapper is designed to work with RILs of 6-th generation or higher where #
#     the number of heterozygotes in the set of scored individuals is low.       #
#     MadMapper is designed to work with high-density genetic maps where         #
#     recombination values between adjacent markers is 0.2 or better (lower).    #
#                                                                                #
#  = WHAT MAD MAPPER DOES =                                                      #
#                            MadMapper analyses raw marker scores (loc file).    #
#                            Based on this analysis it DOES:                     #
#                                                                                # 
#  A. MadMapper generates pairwise distance scores for all markers in a dataset. #
#  B. MadMapper does grouping/clustering analysis based on pairwise distances.   #
#  C. MadMapper assists to find genetic bins (sets of tightly linked markers).   #
#  D. MadMapper assigns new markers in a dataset to known linkage groups.        #
#  E. MadMapper validates all markers in a dataset and sorts them into           #
#     'good' and 'bad' groups based on several criteria and cutoff values.       #
#  F. MadMapper does so called 'triplet' or 'trio' analysis finding sets of      #
#     tightly linked triplets (3 markers in a row) and their relative order.     #
#  G. MadMapper analyzes markers with negative linkage.                          #
#                                                                                #
#  H. MadMapper DOES NOT construct genetic maps yet. However, triplet analysis   #
#     is close enough to this final goal.                                        #
#                                                                                #
#     MadMapper utilizes BIT scoring system instead of LOD scores to manipulate  #
#     genetic data.                                                              #
#  ____________________________________________________________________________  #
#                                                                                #
#  = INPUT FILES =                                                               #
#  1. Locus file with raw marker scores is required.                             #
#     ### README_001 ###                                                         #
#     Frame work marker list (map) is optional.                                  #
#                                                                                #
#  = BASIC DATA STRUCTURE =                                                      #
#     ### README_002 ###                                                         #
#  2. MadMapper reads locus file into memory and creates two-dimensional array:  #
#        sequence_array[id,q] = data_point                                       #
#             where:    id - marker_id                                           #
#                        q - individual RIL                                      #
#               data_point - raw marker score                                    #
#     sequence_array[id,q] contains marker scores 'as is' in the input file:     #
#     'A' 'B' 'C' 'D' 'H' and '-' (heterozygotes are allowed and considered)     #
#                                                                                #
#     ### README_003 ###                                                         #
#  3. At the same MadMapper creates a second two-dimensional array:              #
#        sequence_array_bin[id,q] = score_point                                  #
#     where all heterozygotes are removed or replaced with homozygotes.          #
#     This binary array has 'A' and 'B' values only.                             #
#     'H' is considered as "no data"                                             #
#     'C' -> 'B'    'D' -> 'A'                                                   #
#     There is an assumption that 'C' and 'D' scores should be homozygous for    #
#     RIL set of the 6-th generation or higher.                                  #
#     This 'sequence_array_bin' will be used to build non-redundant data set     #
#     based on marker scores as well as to filter 'bad' scores (see #4 below).   #
#                                                                                #
#     NOTE: If the initial dataset has 'A' and 'B' values only then              #
#     'array sequence_array_bin' is identical to 'sequence_array'                #
#  ____________________________________________________________________________  #
#                                                                                #
#  = NON-REDUNDANT MARKER SCORES =                                               #
#     ### README_004 ###                                                         #
#  4. '*.z_nr_scores.loc' will be generated with non-redundant binary scores.    #
#     If two or more markers have identical scores based on the analysis of the  #    
#     sequence_array_bin[id,q] then only one of them will be written into        #
#     '*.z_nr_scores.loc file'.                                                  #
#     '*.z_scores_dupl' have info about duplicated markers and which marker      #
#     is used as a master marker.                                                #
#     NOTE: THIS NON-REDUNDANT SET IS FILTERED FOR "BAD" MARKERS (SEE BELOW #5). #
#     IT MEANS THAT THIS SET DOES NOT CONTAIN MARKERS WITH ALLELE DISTORTION     #
#     AND MARKERS WITH THE NUMBER OF MISSING DATA ABOVE DEFINED CUTOFF VALUES.   #
#                                                                                #
#  = FILTERING OF BAD MARKERS =                                                  #
#     ### README_005 ###                                                         #
#  5. MadMapper has two variables:                                               #
#     'allele_dist'  = 0.33 ( allele distortion )                                #
#     'abs_loss' = 50       ( absolute data loss )                               #
#     to filter those markers that do not pass criteria:                         #
#     allele ratio (distortion) should be less than 1:3 (0.33)                   #
#     number of missing scores should be less than 50                            #
#     ### README_006 ###                                                         #
#     these cutoff values can be changed using MadMapper arguments/options       #
#     at the time of script execution.                                           #
#     ### README_007 ###                                                         #
#     variable 'nr_good_list' stores marker IDs which passed criteria above.     #
#     Markers IDs stored in 'nr_good_list' will be used for TRIO analysis.       #
#                                                                                #
#     NOTE: FILTERING DOES NOT AFFECT CALCULATION OF PAIRWISE DISTANCES AND      #
#     GLOBAL MATRIX WITH PAIRWISE DISTANCES AND FURTHER CLUSTERING ANALYSIS.     #
#     NON-REDUNDANT CLEAN (FILTERED) SET CAN BE USED ON A SECOND ROUND OF        #
#     THE ANALYSIS IF NECESSARY. HOWEVER, FILTERED NON-REDUNDANT SET IS USED     #
#     FOR TRIPLET ANALYSIS TO FIND BEST TRIO SETS (SEE BELOW).                   #
#  ____________________________________________________________________________  #
#                                                                                #
#  = CALCULATION OF PAIRWISE DISTANCES =                                         #
#     ### README_008 ###                                                         #
#  6. BIT scoring matrix values and REC (recombination) scoring matrix values    #
#     are defined under function 'Define_Bit_Scores'. For each pair of markers   #
#     haplotype pairwise distances are calculated according to this BIT and REC  #
#     matrix scoring system.                                                     #
#                                                                                #
#  = GLOBAL MATRIX WITH PAIRWISE DISTANCES =                                     #
#     ### README_009 ###                                                         #
#  7. Creation of Global Matrix Data for further marker clustering:              #
#     There are three cutoff values to create Matrix with pairwise distances     #
#     for further clustering:                                                    #
#     'rec_cut' (recombination or distance cutoff) 0.2 - default                 #
#     'bit_cut' (BIT score cutoff) 100 - default                                 #
#     'dat_cut' (datapoints cutoff) 25 - default                                 #
#                                                                                #
#             Explanation of datapoints:                                         #
#                 0    1    1    2    2    3    3    4    4    5                 #
#             ====5====0====5====0====5====0====5====0====5====0                 #
#          M1 ----AAABBBAAABBBAAABBBAAABBBAA--------------------                 #
#          M2 ---------BAAABBBAAABBBAAABBBAAAABBB---------------                 #
#          M3 --------------BBAAABBBAAABBBAAAABBBBBAAA----------                 #
#          M4 -------------------BBBAAABBBAAAABBBBBAAAAABBB-----                 #
#          M5 ------------------------ABBBAAAABBBBBAAAAABBBBBAAA                 #
#                                                                                #
#     In this example the number of individuals (RILs) is 50.                    #
#     If 'data_cut' value for this small set is 12 then                          #
#     at least 12 scores between two markers must be compared                    #   
#     to get pairwise matrix values. In this example pairwise scores for pairs   #
#     of markers: M1-M4 M1-M5 M2-M5 are not assigned because they have           #
#     overlapping scores below data_cut cutoff (12).                             #
#                                                                                #
#  = MATRIX OUTPUT FILES [ PAIRWISE SCORES ] =                                   #
#     ### README_010 ###                                                         #
#     '*.pairs_all' global matrix file will contain pairwise scores for all      #
#     markers if data overlap is greater than 'data_cut' value (25 by default).  #
#     BIT score cutoff and REC score cutoff do not affect this file and will be  #
#     used for further clustering (grouping) analysis only.                      #
#     ### README_011 ###                                                         #
#     '*.pairs_positive' will contain pairwise data for all pairs of markers     #
#     with recombination 0.4 or less (positive linkage)                          #
#     ### README_012 ###                                                         #
#     '*.pairs_negative' will contain pairwise data for all pairs of markers     #
#     with recombination 0.6 or greater (negative linkage)                       #
#                                                                                #
#     STRUCTURE OF '*.pairs_all' '*.pairs_positive' '*.pairs_negative' FILES:    #
#    * first column:   marker ID "A"                                             #
#    * second column:  marker ID "B"                                             #
#    * third column:   recombination value between markers "A" and "B"           #
#    * fourth column:  BIT score for pair "A" and "B"                            #
#    * fifth column:   datapoints value                                          #
#                      (fraction of datapoints for pair "A" and "B")             #
#    * sixth column:   "***" - visual mark                                       #
#    * seventh column: total number of recombination events                      #
#    * eighth column:  total number of datapoints                                #
#    * ninth column:   total number of data loss                                 #
#    * tenth column:   total number of possible comparisons (number of RILs)     #
#                                                                                #
#     PAIRWISE SCORES FOR THE EXAMPLE ABOVE:                                     #
#   M1      M1      0.0     156     0.52    ***     0       26      24      50   #
#   M1      M2      0.0     126     0.42    ***     0       21      29      50   #
#   M1      M3      0.0     96      0.32    ***     0       16      34      50   #
#   M2      M1      0.0     126     0.42    ***     0       21      29      50   #
#   M2      M2      0.0     156     0.52    ***     0       26      24      50   #
#   M2      M3      0.0     126     0.42    ***     0       21      29      50   #
#   M2      M4      0.0     96      0.32    ***     0       16      34      50   #
#   M3      M1      0.0     96      0.32    ***     0       16      34      50   #
#   M3      M2      0.0     126     0.42    ***     0       21      29      50   #
#   M3      M3      0.0     156     0.52    ***     0       26      24      50   #
#   M3      M4      0.0     126     0.42    ***     0       21      29      50   #
#   M3      M5      0.0     96      0.32    ***     0       16      34      50   #
#   M4      M2      0.0     96      0.32    ***     0       16      34      50   #
#   M4      M3      0.0     126     0.42    ***     0       21      29      50   #
#   M4      M4      0.0     156     0.52    ***     0       26      24      50   #
#   M4      M5      0.0     126     0.42    ***     0       21      29      50   #
#   M5      M3      0.0     96      0.32    ***     0       16      34      50   #
#   M5      M4      0.0     126     0.42    ***     0       21      29      50   #
#   M5      M5      0.0     156     0.52    ***     0       26      24      50   #
#                                                                                #
#  ____________________________________________________________________________  #
#                                                                                #
#  = CLUSTERING AND GROUPING =                                                   #
#     ### README_014 ###                                                         #
#  8. Regardless of filtering described above (see paragraphs 4 and 5)           #
#     MadMapper will do clustering with all markers from input dataset.          #
#     Clustering with filtered (selected) markers can be done on a second run    #
#     of MadMapper if necessary.                                                 #
#                                                                                #
#     ### README_015 ###                                                         #
#     DFS (Depth First Search) procedure is repeated 16 times with different     #
#     recombination (haplotype distance) cutoff values starting with 0.2 and     #
#     ending with 0.0. Markers are grouped together based on transitive linkage. #
#     If marker 'M1' is linked to marker 'M2' and marker 'M2' is linked to 'M3'  #
#     then all three markers 'M1' 'M2' and 'M3' belong to the same linkage group #
#     even if 'M1' is not linked directly to marker 'M3'.                        #
#                                                                                #
#  = CLUSTERING/GROUPING OUTPUT FILES =                                          #
#     Information about grouping is stored in three files per iteration:         #
#     '*.matrix'     - pairwise distances for a given group                      #
#     '*.adj_list'   - adjacency list                                            #
#     '*.group_info' - group info                                                #
#                                                                                #
#     STRUCTURE OF '*.matrix' FILE:                                              #
#    * first column:   marker ID "A"                                             #
#    * second column:  marker ID "B"                                             #
#    * third column:   recombination value between markers "A" and "B"           #
#    * fourth column:  BIT score for markers "A" and "B"                         #
#    * fifth column:   datapoints value                                          #
#                     (fraction of datapoints for pair "A" and "B")              #
#                                                                                #
#     STRUCTURE OF '*.group_info' FILE:                                          #
#    * first column:   marker ID                                                 #
#    * second column:  length of an adjacency list for given marker or how many  #
#                      other markers are linked to given marker directly         #
#    * third column:   size of the given group (how many markers in this         #
#                      particular linked group)                                  #
#    * fourth column:  arbitrary group number                                    #
#    * fifth column:   visual mark ("*****" separates different group)           #
#    * sixth column:   information about framework markers                       #
#    * seventh column: type of graph [SINGLETON/LINKED/COMPLETE].                #
#      If node is a singleton (is not connected to any other node) then it is    #
#      labeled as 'SINGLE____NODE'                                               #
#      If nodes form complete graph (all nodes linked to each other directly)    #
#      then such group is labeled as 'COMPLETE_GRAPH'                            #
#      If group is not a complete graph (some nodes do not have direct links     #
#      or connections to other nodes then such group is labeled as               #
#      'LINKED___GROUP'.                                                         #
#    * eighth column:  type of node (SATURATED or DILUTED). If a node has all    #
#      possible connections to all other nodes in a group [in other words: node  #
#      is connected directly to all other nodes in a group] then such node is    #
#      labeled as 'SATURATED_NODE'. 'DILUTED___NODE' is an indication that a     #
#      group does not form complete graph. Group can be considered as a bin if   #
#      all nodes in a given group are SATURATED and graph is COMPLETE.           #
#  ____________________________________________________________________________  #
#                                                                                #
#  = CLUSTERING/GROUPING SUMMARY FOR ALL 16 ITERATIONS =                         #
#     [ DENDRO-CLUSTERING ]                                                      #
#     STRUCTURE OF *.x_tree_clust FILE:                                          #
#    * 1-st column: group ID for clustering with cutoff 0.20                     #
#    * 2-nd column: group ID for clustering with cutoff 0.18                     #
#    * 3-d  column: group ID for clustering with cutoff 0.16                     #
#    * 4-th column: group ID for clustering with cutoff 0.14                     #
#    * 5-th column: group ID for clustering with cutoff 0.12                     #
#    * 6-th column: group ID for clustering with cutoff 0.10                     #
#    * 7-th column: group ID for clustering with cutoff 0.09                     #
#    * 8-th column: group ID for clustering with cutoff 0.08                     #
#    * 9-th column: group ID for clustering with cutoff 0.07                     #
#    *10-th column: group ID for clustering with cutoff 0.06                     #
#    *11-th column: group ID for clustering with cutoff 0.05                     #
#    *12-th column: group ID for clustering with cutoff 0.04                     #
#    *13-th column: group ID for clustering with cutoff 0.03                     #
#    *14-th column: group ID for clustering with cutoff 0.02                     #
#    *15-th column: group ID for clustering with cutoff 0.01                     #
#    *16-th column: group ID for clustering with cutoff 0.00                     #
#    *17-th column: type of graph (COMPLETE or LINKED) for the last iteration    #
#    *18-th column: type of node (SATURATED or DILUTED) for the last iteration   #
#    *19-th column: '***' - visual mark                                          #
#    *20-th column: linkage group from frame work marker list                    #
#    *21-st column: position on the map of frame marker (if available)           #
#    *22-nd column: '***' - visual mark                                          #
#    *23-d  column: "ABC" (alphabetical) order of markers prior clustering       #
#    *24-th column: '***' - visual mark                                          #
#    *25-th column: "LG" - reserved field for manipulation in excel like editor  #
#    *26-th column: marker ID                                                    #
#    *27-th column: '*' - visual mark                                            #
#    *28-th column: sum of 'A' scores (A+D)                                      #
#    *29-th column: sum of 'B' scores (B+C)                                      #
#    *30-th column: total number of possible scores (number of RILs)             #
#    *31-st column: '*' - visual mark                                            #
#    *32-nd column: order of markers according to dendro-clustering              #
#                                                                                #
#  ____________________________________________________________________________  #
#                                                                                #
#  = TRIO/TRIPLET ANALYSIS =                                                     #
#                                                                                #
#  9. MadMapper performs TRIPLET or TRIO analysis. It finds tightly linked       #
#     triplets (3 markers in a row) and their relative order.                    #
#                                                                                #
#     Each marker from 'GOOD NON-REDUNDANT SET' (see #5 above) is checked for    #
#     all possible combinations with other markers. In this case 'test' marker   #
#     takes a middle position in a trio and recombination scores between middle  #
#     and flanking markers are analyzed.                                         #
#                                                                                #
#    DATA STRUCTURE FOR TRIO ANALYSIS:                                           ####################
#                                                                                                   #
#    ALL TRIOS (madmapper_scores.loc.out.z_trio_all):                                               #
#                                                                                                   #
# MV  0.1     240   1.0    MP   0.06    264   1.0    MI  ***  0.16    204   1.0   ***  0            #
# MV  0.1     240   1.0    MP   0.0455  240   0.88   MK  ***  0.1591  180   0.88  ***  0            #
# MV  0.1     240   1.0    MP   0.04    276   1.0    MR  ***  0.06    264   1.0   ***  2            #
# MV  0.1     240   1.0    MP   0.08    252   1.0    MS  ***  0.02    288   1.0   ***  4            #
# ..............................                                                                    #
# MG  0.125   216   0.96   MR   0.14    216   1.0    MT  ***  0.2708  132   0.96  ***  0            #
# MG  0.125   216   0.96   MR   0.06    264   1.0    MV  ***  0.1875  180   0.96  ***  0            #
# MH  0.1111  210   0.9    MR   0.1458  204   0.96   MF  ***  0.0444  246   0.9   ***  5            #
# MH  0.1111  210   0.9    MR   0.125   216   0.96   MG  ***  0.0222  258   0.9   ***  5            #
#                                                                                                   #
#                                                                                                   #
#    GOOD TRIOS (madmapper_scores.loc.out.z_trio_good):                                             #
#                                                                                                   #
# MK  0.0909  216   0.88   MR   0.06    264   1.0    MV  ***  0.1591  180   0.88  ***  0  ---  0.15 #
# MP  0.04    276   1.0    MR   0.04    276   1.0    MS  ***  0.08    252   1.0   ***  0  +++  0.0  #
# MP  0.04    276   1.0    MR   0.14    216   1.0    MT  ***  0.18    192   1.0   ***  0  ---  0.2  #
# MP  0.04    276   1.0    MR   0.06    264   1.0    MV  ***  0.1     240   1.0   ***  0  ---  0.04 #
#                                                                                                   #
#                                                                                                   #
#    BEST TRIOS (madmapper_scores.loc.out.z_trio_best)                                              #
#                                                                                                   #
# MH  0.0     234   0.78   MK   0.0     264   0.88   MI  ***  0.0     270   0.9   ***  0  +++  0.0  #
# MI  0.0     264   0.88   MK   0.0     234   0.78   MH  ***  0.0     270   0.9   ***  0  +++  0.0  #
# MK  0.0455  240   0.88   MP   0.04    276   1.0    MR  ***  0.0909  216   0.88  ***  0  +++  0.0  #
# MR  0.04    276   1.0    MP   0.0455  240   0.88   MK  ***  0.0909  216   0.88  ***  0  +++  0.0  #
# MP  0.04    276   1.0    MR   0.04    276   1.0    MS  ***  0.08    252   1.0   ***  0  +++  0.0  #
# MS  0.04    276   1.0    MR   0.04    276   1.0    MP  ***  0.08    252   1.0   ***  0  +++  0.0  #
#                                                                                                   #
#                                                                                                   #
#    BAD TRIOS (madmapper_scores.loc.out.z_trio_bad)                                                #
#                                                                                                   #
# MS  0.1     240   1.0    MT   0.12    228   1.0    MV  ***  0.02    288   1.0   ***  5  +++  0.0  #
# MV  0.12    228   1.0    MT   0.1     240   1.0    MS  ***  0.02    288   1.0   ***  5  +++  0.0  #
#   |       |     |      |    |       |     |      |   |     |      |     |     |    |   |    |     #
# 1 |  2    |  3  |  4   | 5  |  6    |  7  |  8   | 9 | 10  | 11   |  12 |  13 | 14 | 15  16 | 17  #
# -    -       -     -     -     -       -     -     -   --    --      --    --   --   --  --   --  #
#                                                                                                   #
#                                                                                ####################
#   * 1-st column - 'upper' marker ID (above 'target' marker)                    #
#                                                                                #
#   * 2-nd column - distance  between 'upper' marker and 'target' [DIST_UP]      #
#                                                                                #
#   * 3-d  column - BIT score between 'upper' marker and 'target'                #
#                                                                                #
#   * 4-th column - fraction of datapoints between 'upper' marker and 'target'   #
#                                                                                #
#   * 5-th column - 'target' marker ID (middle position in a triplet)            #
#                                                                                #
#   * 6-th column - distance  between 'lower' marker and 'target' [DIST_DN]      #
#                                                                                #
#   * 7-th column - BIT score between 'lower' marker and 'target'                #
#                                                                                #
#   * 8-th column - fraction of datapoints between 'lower' marker and 'target'   #
#                                                                                #
#   * 9-th column - 'lower' marker ID (below 'target' marker)                    #
#                                                                                #
#   *10-th column - '***' - visual mark                                          #
#                                                                                #
#   *11-th column - distance  between 'upper' marker and 'lower' marker          #
#                  [distance  between flanking markers]           [DIST_FL]      #
#                                                                                #
#   *12-th column -  BIT score between flanking markers                          #
#                                                                                #
#   *13-th column - fraction of datapoints between flanking markers              #
#                                                                                #
#   *14-th column - '***' - visual mark                                          #
#                                                                                #
#   *15-th column - number of double cross-overs in a triplet                    #
#                                                                                #
#   *16-th column - '+++' stands for 'best' case ('---' for other than 'best')   #
#                                                                                #
#   *17-th column - trio distance summary to calculate the 'best' case           #
#                   it is a sum of three values:                                 #
#                   [DIST_UP] + [DIST_DN] + [DIST_FL] minus [BEST_VALUE]         #
#        ( it must be 0.0 for 'best' case and greater than 0.0 for other cases ) #
#        [BEST_VALUE] is defined as a minimum (lowest value) of sum of           #
#                     [DIST_UP] + [DIST_DN] + [DIST_FL] for 'best' case          #
#                                                                                #
#                                                                                #
#  = BRIEF DESCRIPTION OF TRIO APPROACH =                                        #
#                                                                                #
#       ### README_016 ###                                                       #
#   [A] - All possible trios are analyzed and TRIO MATRIX is created             #
#                                                                                #
#       ### README_017 ###                                                       #
#   [B] - best_recomb_trio[item_c,dk] array contains trios that have a number    #
#         of double crossovers less than cutoff value (3 by default)             #
#         ['double_cross' is a last argument/option of MadMapper program]        #
#                                                                                #
#       ### README_018 ###                                                       #
#   [C] - Best recombination value is calculated for each 'target' marker        #
#            [BEST_VALUE] = MIN OF SUM: [DIST_UP] + [DIST_DN] + [DIST_FL]        #
#                                                                                #
#       ### README_019 ###                                                       #
#   [D] - Extraction of BEST TRIOS into '*.z_trio_best' file                     #
#                                                                                #
#  ____________________________________________________________________________  #
#                                                                                #
#                                                                   NUMBER OF    #
#  = BRIEF DESCRIPTION OF ALL OUTPUT FILES =                           FILES:    #
#                                                                                #
#   PAIRWISE DISTANCES AND CLUSTERING:                                           #
#    *.adj_list_(01-16) - adjacency lists (positive clustering)            16    #
#    *.adj_list_N(1-3)  - adjacency lists (negative clustering)             3    #
#    *.group_info_(01-16) - group info (positive clustering)               16    #
#    *.group_info_N(1-3)  - group info (negative clustering)                3    #
#    *.group_info_Summary -                                                 1    #
#               -  summary for all 16 iterations of  positive clustering         #
#    *.matrix_(01-16)   - pairwise distance matrix (positive clustering)   16    #
#    *.matrix_N(1-3)    - pairwise distance matrix (negative clustering)    3    #
#    *.pairs_all    - pairwise distance matrix with all available scores    1    #
#    *.pairs_negative - pairwise distance matrix with all negative scores   1    #
#    *.pairs_positive - pairwise distance matrix with all positive scores   1    #
#                                                                                #
#   ERROR CHECKING:                                                              #
#    *.set_dupl       - duplicated marker IDs  | DO NOT CONFUSE WITH |      1    #
#    *.set_uniq       - unique marker IDs      |  DUPLICATED SCORES  |      1    #
#                                                                                #
#   LOG FILE:                                                                    #
#    *.x_log_file     - log file with run parameters recorded               1    #
#                                                                                #
#   MARKER SCORES INFO:                                                          #
#    *.x_scores_stat  - detailed information about scores/linkage stat      1    #
#                       (see *.x_scores_stat structure below README_013)         #
#                                                                                #
#   GROUPING OF MARKERS ( DENDRO-CLUSTERING ):                                   #
#    *.x_tree_clust   - grouping of markers based on the analysis           1    #
#                       of all 16 *.group_info_(01-16) files                     #
#                                                                                #
#   NON-REDUNDANT SCORES SET AND FILTERING OF 'BAD' MARKERS:                     #
#    *.z_marker_sum                                                         1    #
#    *.z_nr_scores.loc                                                      1    #
#    *.z_scores_dupl                                                        1    #
#                                                                                #
#   TRIO (TRIPLET) ANALYSIS:                                                     #
#    *.z_trio_all                                                           1    #
#    *.z_trio_bad                                                           1    #
#    *.z_trio_best                                                          1    #
#    *.z_trio_good                                                          1    #
#    *.z_trio_graph                                                         1    #
#    *.z_trio_map                                                           1    #
#                                                                                #
#   BIN CONSENSUS SET (EXPERIMENTAL):                                            #
#    *.z_xconsensus.conv_adjc                                               1    #
#    *.z_xconsensus.conv_real                                               1    #
#    *.z_xconsensus.debug                                                   1    #
#    *.z_xconsensus.frame                                                   1    #
#    *.z_xconsensus.loc_all                                                 1    #
#    *.z_xconsensus.loc_nr                                                  1    #
#    *.z_xconsensus.matrix                                                  1    #
#                                                                                #
#                                          TOTAL NUMBER OF OUTPUT FILES:   82    #
#  ____________________________________________________________________________  #
#                                                                                #
#   MARKER SCORES INFO OUTPUT FILE FORMAT:                                       #
#     ### README_013 ###                                                         #
#    * 1-st column (MARKER_ID): marker ID                                        #
#    * 2-nd column (00_01): number of markers linked to a given marker within    #
#                   0.0 - 0.1 recombination frequency (strong positive linkage)  #
#    * 3-d column  (01_02): number of markers linked to a given marker within    #
#                   0.1 - 0.2 recombination frequency (strong positive linkage)  #
#    * 4-th column (02_03): number of markers linked to a given marker within    #
#                   0.2 - 0.3 recombination frequency (positive linkage)         #
#    * 5-th column (03_04): number of markers linked to a given marker within    #
#                   0.3 - 0.4 recombination frequency                            #
#    * 6-th column (04_05): number of markers linked to a given marker within    #
#                   0.4 - 0.5 recombination frequency                            #
#    * 7-th column (05_06): number of markers linked to a given marker within    #
#                   0.5 - 0.6 recombination frequency                            #
#    * 8-th column (06_07): number of markers linked to a given marker within    #
#                   0.6 - 0.7 recombination frequency                            #
#    * 9-th column (07_08): number of markers linked to a given marker within    #
#                   0.7 - 0.8 recombination frequency (negative linkage)         #
#    *10-th column (08_09): number of markers linked to a given marker within    #
#                   0.8 - 0.9 recombination frequency (strong negative linkage)  #
#    *11-th column (09_10): number of markers linked to a given marker within    #
#                   0.9 - 1.0 recombination frequency (strong negative linkage)  #
#    *12-th column (REC:P-N): difference between strong positive and strong      #
#                   negative events [ (00_01+01_02) minus (08_09+09_10) ]        #
#    *13-th column (REC_ABS): absolute recombination value for a given marker    #
#                  (accumulative value for all possible recombination values)    #
#                  [ Summary of all (0.5 - recombination frequency) ]            #
#    *14-th column (***): "***" - visual mark                                    #
#    *15-th column (BIT_POS): number of markers having BIT score 100 or higher   #
#                             to a given marker (positive BIT scores)            #
#    *16-th column (BIT_MED): number of markers having BIT score within          #
#                             100 and -100 range (low or medium BIT scores)      #
#    *17-th column (BIT_NEG): number of markers having BIT score -100 or lower   #
#                             to a given marker (negative BIT scores)            #
#    *18-th column (BIT:P-N): difference between positive and negative BIT       #
#                             scores [ BIT_POS minus BIT_NEG ]                   #
#    *19-th column (BIT_ABS): absolute BIT value (score) for a given marker      #
#                             (accumulative value for all possible BIT scores)   #
#                             [ Summary of all BIT scores ]                      #
#    *20-th column (***): "***" - visual mark                                    #
#     ### README_13LG ###                                                        #
#    *21-st column (LG_SUMMARY): attempt to classify markers by analysis of      #
#                  2-nd through 11-th columns [can be ignored currently]         #
#    *22-nd column (NP_SUMMARY): brief summary about positive/negative dominance #
#    *23-d  column (***): "***" - visual mark                                    #
#    *24-th to 33-d columns: stat summary for 'AD-BC-H' scores ('X' - no data)   #
#                                                                                #
##################################################################################

##################################################################################
#                                                                                #
#  = EXAMPLES OF DATA PROCESSING BY MAD MAPPER =                                 #
#                                                                                #
#   *** INPUT RAW MARKER SCORES:                                                 #
#                                                                                #
#    M1 ----AAABBBAAABBBAAABBBAAABBBAA--------------------                       #
#    M2 ---------BAAABBBAAABBBAAABBBAAAABBB---------------                       #
#    M3 --------------BBAAABBBAAABBBAAAABBBBBAAA----------                       #
#    M4 -------------------BBBAAABBBAAAABBBBBAAAAABBB-----                       #
#    M5 ------------------------ABBBAAAABBBBBAAAAABBBBBAAA                       #
#    M6 -------------------BBBAAABBBAAAABBBBBAAAAABBBBBAAB                       #
#    M7 ------------------ABBBAAABBBAAAABBBBBAAAAABBBBBABB                       #
#    M8 ----------------AAABBBAAABBBAAAAABBBBAAAAABBBBBBBB                       #
#    M9 ----------------AAABBBAAABBBAAAAAABBBAAAAABBBBBBBB                       #
#    MF -AABBB-AABBBBBABBBBAAAAAAAAAAAAAAAABBBBBBBBBBBBBBB                       #
#    MG -AABBB-AABBBBAABBBBAAAAAAAAAAAAAAAABBBBBBBBBBBBBBB                       #
#    MH -AABBB-AABBBAAA-BBBAAAAAAAA-AAAAAAABBBBBBBB-BBBBBB                       #
#    MI AAABBBAAABBBAAABBBBAAAAAAAAAAAAAAAABBBBBBBBBBBBBBB                       #
#    MJ AAABBBAAABBBAAABBBBAAAAAAAAAAAAAAAABBBBBBBBBBBBBBB                       #
#    MK AA-BBBA-ABBBAAABBBB---AAAAAAAAAAAAABBBBBBBBBBBBBB-                       #
#    ML DD-CCCD-DCCCDDDCCCC---DDDDDDDDDDDDDCCCCCCCCCCCCCC-                       #
#    MN BBBAAABBBAAABBBAAAABBBBBBBBBBBBBBBBAAAAAAAAAAAAAAB                       #
#    MP AAABBBAAABBBAAABBBBAAAAAAAAAAAAAAAABBBBBBBBBBBBAAA                       #
#    MR AAABBBAAABBBAAABBBBAAAAAAAAAAAAAAAABBBBBBBBBBAAAAA                       #
#    MS AAABBBAAABBBAAABBBBAAAAAAAAAAAAAAAABBBBBBBBAAAAAAA                       #
#    MT ABABBBAAABABAAABBBBAABAAAABAAAAAAAABBBBABBBAAAAAAA                       #
#    MV AAABBBAAABBBAAABBBBAAAAAAAAAAAAAAAABBBBBBBAAAAAAAA                       #
#    MW -AAAAAAAAAAAAAAAAAAAAAAAABAAAAAAAAAAAAAAAAAABBBBBB                       #
#    MX ------------A------------B------------A----------B                       #
#    MY -----------BA-----------AB-----------AA---------BB                       #
#    MZ ------------AA----------ABB-----------AAA--------B                       #
#                                                                                #
#                                                                                #
#   *** PROGRAM EXECUTION:                                                       #
#    $python Python_MadMapper_V248_RECBIT.py madmapper_scores.loc                #
#                                            madmapper_scores.loc.out            #
#                                            0.2 100 12 X 0.33 25 TRIO 3         #
#                                                                                #
#                                                                                #
#   *** X_LOG (madmapper_scores.loc.out.x_log_file) OUTPUT FILE:                 #
#    =============================================                               #
#        RUN PARAMETERS:                                                         #
#     1. INPUT  FILE:  madmapper_scores.loc                                      #
#     2. OUTPUT FILE:  madmapper_scores.loc.out                                  #
#     3. RECM CUTOFF:  0.2                                                       #
#     4. BITS CUTOFF:  100                                                       #
#     5. DATA CUTOFF:  12                                                        #
#     6. ALLELE DIST:  0.33                                                      #
#     7. MISSING DAT:  25                                                        #
#     8. FRAME  LIST:  X                                                         #
#     9. TRIO ANALYS:  TRUE                                                      #
#    10. DOUBLE LIMT:  3                                                         #
#    =============================================                               #
#    =============================================                               #
#    26 UNIQ IDs IN THE SET FOUND                                                #
#    0 IDs ARE DUPLICATED                                                        #
#    =============================================                               #
#    CONTINUE ANALYSIS WITH 26 SEQUENCES OUT OF 26                               #
#    =============================================                               #
#    7  GROUPS WERE FOUND (STEP 01)                                              #
#    7  GROUPS WERE FOUND (STEP 02)                                              #
#    ........                                                                    #
#                                                                                #
#                                                                                #
#   *** NON-REDUNDANT SCORES (madmapper_scores.loc.out.z_nr_scores.loc):         #
#                                                                                #
#    M1 ----AAABBBAAABBBAAABBBAAABBBAA--------------------                       #
#    M2 ---------BAAABBBAAABBBAAABBBAAAABBB---------------                       #
#    M3 --------------BBAAABBBAAABBBAAAABBBBBAAA----------                       #
#    M4 -------------------BBBAAABBBAAAABBBBBAAAAABBB-----                       #
#    M5 ------------------------ABBBAAAABBBBBAAAAABBBBBAAA                       #
#    M6 -------------------BBBAAABBBAAAABBBBBAAAAABBBBBAAB                       #
#    M7 ------------------ABBBAAABBBAAAABBBBBAAAAABBBBBABB                       #
#    M8 ----------------AAABBBAAABBBAAAAABBBBAAAAABBBBBBBB                       #
#    M9 ----------------AAABBBAAABBBAAAAAABBBAAAAABBBBBBBB                       #
#    MF -AABBB-AABBBBBABBBBAAAAAAAAAAAAAAAABBBBBBBBBBBBBBB                       #
#    MG -AABBB-AABBBBAABBBBAAAAAAAAAAAAAAAABBBBBBBBBBBBBBB                       #
#    MH -AABBB-AABBBAAA-BBBAAAAAAAA-AAAAAAABBBBBBBB-BBBBBB                       #
#    MI AAABBBAAABBBAAABBBBAAAAAAAAAAAAAAAABBBBBBBBBBBBBBB                       #
#    MK AA-BBBA-ABBBAAABBBB---AAAAAAAAAAAAABBBBBBBBBBBBBB-                       #
#    MN BBBAAABBBAAABBBAAAABBBBBBBBBBBBBBBBAAAAAAAAAAAAAAB                       #
#    MP AAABBBAAABBBAAABBBBAAAAAAAAAAAAAAAABBBBBBBBBBBBAAA                       #
#    MR AAABBBAAABBBAAABBBBAAAAAAAAAAAAAAAABBBBBBBBBBAAAAA                       #
#    MS AAABBBAAABBBAAABBBBAAAAAAAAAAAAAAAABBBBBBBBAAAAAAA                       #
#    MT ABABBBAAABABAAABBBBAABAAAABAAAAAAAABBBBABBBAAAAAAA                       #
#    MV AAABBBAAABBBAAABBBBAAAAAAAAAAAAAAAABBBBBBBAAAAAAAA                       #
#  NOTE THAT MARKERS:    MJ and ML were removed because of scores duplication    #
#                        MW was removed because of allele distortion             #
#                        MX MY MZ were removed because of missing data           #
#                                                                                #
#                                                                                #
#   *** INFO ABOUT DUPLICATED SCORES (madmapper_scores.loc.out.z_scores_dupl):   #
#    MI      ==      MJ      ***     MI                                          #
#    MJ      ==      MI      ***     MI                                          #
#    MK      ==      ML      ***     MK                                          #
#    ML      ==      MK      ***     MK                                          #
#     (last column - 'master' marker)                                            #
#                                                                                #
#                                                                                #
#   *** GROUP INFO FILE (iteration #16 madmapper_scores.loc.out.group_info_16):  #
#                                                                                #
# M1   1    7     1     *****   _NONE_  LINKED___GROUP_00001    DILUTED___NODE   #
# M2   3    7     1     -----   _NONE_  LINKED___GROUP_00001    DILUTED___NODE   #
# M3   4    7     1     -----   _NONE_  LINKED___GROUP_00001    DILUTED___NODE   #
# M4   4    7     1     -----   _NONE_  LINKED___GROUP_00001    DILUTED___NODE   #
# M5   1    7     1     -----   _NONE_  LINKED___GROUP_00001    DILUTED___NODE   #
# M6   2    7     1     -----   _NONE_  LINKED___GROUP_00001    DILUTED___NODE   #
# M7   3    7     1     -----   _NONE_  LINKED___GROUP_00001    DILUTED___NODE   #
# M8   0    1     2     *****   _NONE_  SINGLE____NODE_00002    SATURATED_NODE   #
# M9   0    1     3     *****   _NONE_  SINGLE____NODE_00003    SATURATED_NODE   #
# MF   0    1     4     *****   _NONE_  SINGLE____NODE_00004    SATURATED_NODE   #
# MG   0    1     5     *****   _NONE_  SINGLE____NODE_00005    SATURATED_NODE   #
# MH   4    5     6     *****   _NONE_  COMPLETE_GRAPH_00006    SATURATED_NODE   #
# MI   4    5     6     -----   _NONE_  COMPLETE_GRAPH_00006    SATURATED_NODE   #
# MJ   4    5     6     -----   _NONE_  COMPLETE_GRAPH_00006    SATURATED_NODE   #
# MK   4    5     6     -----   _NONE_  COMPLETE_GRAPH_00006    SATURATED_NODE   #
# ML   4    5     6     -----   _NONE_  COMPLETE_GRAPH_00006    SATURATED_NODE   #
# MN   0    1     7     *****   _NONE_  SINGLE____NODE_00007    SATURATED_NODE   #
# MP   0    1     8     *****   _NONE_  SINGLE____NODE_00008    SATURATED_NODE   #
# MR   0    1     9     *****   _NONE_  SINGLE____NODE_00009    SATURATED_NODE   #
# MS   0    1    10     *****   _NONE_  SINGLE____NODE_00010    SATURATED_NODE   #
# MT   0    1    11     *****   _NONE_  SINGLE____NODE_00011    SATURATED_NODE   #
# MV   0    1    12     *****   _NONE_  SINGLE____NODE_00012    SATURATED_NODE   #
# MW   0    1    13     *****   _NONE_  SINGLE____NODE_00013    SATURATED_NODE   #
# MX   0    1    14     *****   _NONE_  SINGLE____NODE_00014    SATURATED_NODE   #
# MY   0    1    15     *****   _NONE_  SINGLE____NODE_00015    SATURATED_NODE   #
# MZ   0    1    16     *****   _NONE_  SINGLE____NODE_00016    SATURATED_NODE   #
#                                                                                #
#   NOTE THE MAJOR DIFFERENCE BETWEEN GROUP 1 AND GROUP 6:                       #
#       GROUP 1 is a 'LINKED GROUP'                                              #
#       GROUP 6 is a 'COMPLETE GRAPH'                                            #
#                                                                                #
#                                                                                #
#   *** CONSENSUS (madmapper_scores.loc.out.z_xconsensus.loc_nr):                #
#                                                                                #
#    GC_00001 --------------------------------------------------                 #
#    GC_00008 ----------------AAABBBAAABBBAAAAABBBBAAAAABBBBBBBB                 #
#    GC_00009 ----------------AAABBBAAABBBAAAAAABBBAAAAABBBBBBBB                 #
#    GC_00010 -AABBB-AABBBBBABBBBAAAAAAAAAAAAAAAABBBBBBBBBBBBBBB                 #
#    GC_00011 -AABBB-AABBBBAABBBBAAAAAAAAAAAAAAAABBBBBBBBBBBBBBB                 #
#    GC_00012 AAABBBAAABBBAAABBBBAAAAAAAAAAAAAAAABBBBBBBBBBBBBBB                 #
#    GC_00015 BBBAAABBBAAABBBAAAABBBBBBBBBBBBBBBBAAAAAAAAAAAAAAB                 #
#    GC_00016 AAABBBAAABBBAAABBBBAAAAAAAAAAAAAAAABBBBBBBBBBBBAAA                 #
#    GC_00017 AAABBBAAABBBAAABBBBAAAAAAAAAAAAAAAABBBBBBBBBBAAAAA                 #
#    GC_00018 AAABBBAAABBBAAABBBBAAAAAAAAAAAAAAAABBBBBBBBAAAAAAA                 #
#    GC_00019 ABABBBAAABABAAABBBBAABAAAABAAAAAAAABBBBABBBAAAAAAA                 #
#    GC_00020 AAABBBAAABBBAAABBBBAAAAAAAAAAAAAAAABBBBBBBAAAAAAAA                 #
#                                                                                #
#                                                                                #
#   *** CONSENSUS MATRIX (madmapper_scores.loc.out.z_xconsensus.matrix):         #
#                                                                                #
#          GC_00008        M8      1.0                                           #
#          GC_00009        M9      1.0                                           #
#          GC_00010        MF      1.0                                           #
#          GC_00011        MG      1.0                                           #
#          GC_00012        MH      1.0                                           #
#          GC_00012        MJ      1.0                                           #
#          GC_00012        MK      1.0                                           #
#          GC_00012        ML      1.0                                           #
#          GC_00012        MI      1.0                                           #
#          GC_00013        MI      1.0                                           #
#          GC_00013        MK      1.0                                           #
#          GC_00013        MH      1.0                                           #
#          GC_00013        MJ      1.0                                           #
#          GC_00013        ML      1.0                                           #
#          GC_00014        MK      1.0                                           #
#          GC_00014        MI      1.0                                           #
#          GC_00014        MH      1.0                                           #
#          GC_00014        ML      1.0                                           #
#          GC_00014        MJ      1.0                                           #
#          GC_00015        MN      1.0                                           #
#          GC_00016        MP      1.0                                           #
#          GC_00017        MR      1.0                                           #
#          GC_00018        MS      1.0                                           #
#          GC_00019        MT      1.0                                           #
#          GC_00020        MV      1.0                                           #
#                                                                                #
#    NOTE THAT GC_00012, GC_00013 and GC_00014 are identical (redundant).        #
#    This is a reason why only GC_00012 was written into non-redundant loc file  #
#                                                                                #
#    Markers M1 - M7 were unable to generate a consensus because they did not    #
#    form 'COMPLETE GRAPH'                                                       #
#                                                                                #
#                                                                                #
#   *** CONSENSUS EXAMPLE (madmapper_scores.loc.out.z_xconsensus.debug):         #
#    ........                                                                    #
#              ===================================================               #
#    MH         -AABBB-AABBBAAA-BBBAAAAAAAA-AAAAAAABBBBBBBB-BBBBBB               #
#    MJ         AAABBBAAABBBAAABBBBAAAAAAAAAAAAAAAABBBBBBBBBBBBBBB               #
#    MK         AA-BBBA-ABBBAAABBBB---AAAAAAAAAAAAABBBBBBBBBBBBBB-               #
#    ML         AA-BBBA-ABBBAAABBBB---AAAAAAAAAAAAABBBBBBBBBBBBBB-               #
#    MI         AAABBBAAABBBAAABBBBAAAAAAAAAAAAAAAABBBBBBBBBBBBBBB               #
#    GC_00012   AAABBBAAABBBAAABBBBAAAAAAAAAAAAAAAABBBBBBBBBBBBBBB               #
#              ===================================================               #
#                                                                                #
#                                                                                #
##################################################################################

def Define_Bit_Scores():

	###################

	global bit_score_AA
	global bit_score_AB
	global bit_score_AC
	global bit_score_AD
	global bit_score_AH

	global bit_score_BA
	global bit_score_BB
	global bit_score_BC
	global bit_score_BD
	global bit_score_BH

	global bit_score_CA
	global bit_score_CB
	global bit_score_CC
	global bit_score_CD
	global bit_score_CH

	global bit_score_DA
	global bit_score_DB
	global bit_score_DC
	global bit_score_DD
	global bit_score_DH

	global bit_score_HA
	global bit_score_HB
	global bit_score_HC
	global bit_score_HD
	global bit_score_HH

	###################

	global rec_score_AA
	global rec_score_AB
	global rec_score_AC
	global rec_score_AD
	global rec_score_AH

	global rec_score_BA
	global rec_score_BB
	global rec_score_BC
	global rec_score_BD
	global rec_score_BH

	global rec_score_CA
	global rec_score_CB
	global rec_score_CC
	global rec_score_CD
	global rec_score_CH

	global rec_score_DA
	global rec_score_DB
	global rec_score_DC
	global rec_score_DD
	global rec_score_DH

	global rec_score_HA
	global rec_score_HB
	global rec_score_HC
	global rec_score_HD
	global rec_score_HH

	###################

	bit_score_AA =  6
	bit_score_AB = -6
	bit_score_AC = -4
	bit_score_AD =  4
	bit_score_AH = -2

	bit_score_BA = -6
	bit_score_BB =  6
	bit_score_BC =  4
	bit_score_BD = -4
	bit_score_BH = -2

	bit_score_CA = -4
	bit_score_CB =  4
	bit_score_CC =  4
	bit_score_CD = -4
	bit_score_CH =  0

	bit_score_DA =  4
	bit_score_DB = -4
	bit_score_DC = -4
	bit_score_DD =  4
	bit_score_DH =  0

	bit_score_HA = -2
	bit_score_HB = -2
	bit_score_HC =  0
	bit_score_HD =  0
	bit_score_HH =  2

	#################

	rec_score_AA =  0
	rec_score_AB =  1
	rec_score_AC =  1
	rec_score_AD =  0
	rec_score_AH =  0.5

	rec_score_BA =  1
	rec_score_BB =  0
	rec_score_BC =  0
	rec_score_BD =  1
	rec_score_BH =  0.5

	rec_score_CA =  1
	rec_score_CB =  0
	rec_score_CC =  0
	rec_score_CD =  1
	rec_score_CH =  0

	rec_score_DA =  0
	rec_score_DB =  1
	rec_score_DC =  1
	rec_score_DD =  0
	rec_score_DH =  0

	rec_score_HA =  0.5
	rec_score_HB =  0.5
	rec_score_HC =  0
	rec_score_HD =  0
	rec_score_HH =  0

	#################

def Read_Data_File(in_name, out_name, rec_cut, bit_cut, dat_cut, frame_file_name):

	###################

	global bit_score_AA
	global bit_score_AB
	global bit_score_AC
	global bit_score_AD
	global bit_score_AH

	global bit_score_BA
	global bit_score_BB
	global bit_score_BC
	global bit_score_BD
	global bit_score_BH

	global bit_score_CA
	global bit_score_CB
	global bit_score_CC
	global bit_score_CD
	global bit_score_CH

	global bit_score_DA
	global bit_score_DB
	global bit_score_DC
	global bit_score_DD
	global bit_score_DH

	global bit_score_HA
	global bit_score_HB
	global bit_score_HC
	global bit_score_HD
	global bit_score_HH

	###################

	global rec_score_AA
	global rec_score_AB
	global rec_score_AC
	global rec_score_AD
	global rec_score_AH

	global rec_score_BA
	global rec_score_BB
	global rec_score_BC
	global rec_score_BD
	global rec_score_BH

	global rec_score_CA
	global rec_score_CB
	global rec_score_CC
	global rec_score_CD
	global rec_score_CH

	global rec_score_DA
	global rec_score_DB
	global rec_score_DC
	global rec_score_DD
	global rec_score_DH

	global rec_score_HA
	global rec_score_HB
	global rec_score_HC
	global rec_score_HD
	global rec_score_HH

	###################

	global deep_clustering
	global map_construction
	global double_limit
	global allele_dist
	global abs_loss
	global round_scale
	global print_all_pairs

	###################

	print "============================================="
	print "    RUN PARAMETERS:  "
	print " 1. INPUT  FILE:  " + in_name
	print " 2. OUTPUT FILE:  " + out_name
	print " 3. RECM CUTOFF:  " + str(rec_cut)
	print " 4. BITS CUTOFF:  " + str(bit_cut)
	print " 5. DATA CUTOFF:  " + str(dat_cut)
	print " 6. ALLELE DIST:  " + str(allele_dist)
	print " 7. MISSING DAT:  " + str(abs_loss)
	print " 8. FRAME  LIST:  " + frame_file_name
	print " 9. TRIO ANALYS:  " + map_construction
	print "10. DOUBLE LIMT:  " + str(double_limit)

	in_file  = open(in_name,  "rb")
	if map_construction != "NR_SET":
		out_file0 = open(out_name + '.pairs_positive', "wb")
		out_file1 = open(out_name + '.pairs_all', "wb")
		out_file2 = open(out_name + '.pairs_negative', "wb")
	############# GOOD AND BAD SET ###############
	out_file4 = open(out_name + '.set_uniq', "wb")
	out_file5 = open(out_name + '.set_dupl', "wb")
	##############################################
	out_file6 = open(out_name + '.x_log_file', "wb")
	##############################################
	if map_construction != "NR_SET":
		out_file7 = open(out_name + '.x_scores_stat', "wb")
	##############################################
	global out_file8
	if map_construction != "NR_SET":
		out_file8 = open(out_name + '.x_tree_clust', "wb")
	##############################################
	# if map_construction == "TRUE":
	out_file3a = open(out_name + '.z_nr_scores.loc', "wb")
	out_file3b = open(out_name + '.z_scores_dupl', "wb")
	if map_construction != "NR_SET":
		out_file3i = open(out_name + '.z_marker_sum', "wb")
	if deep_clustering  == "TRUE" and map_construction != "NR_SET":
		global out_file3k
		global out_file3la
		global out_file3lb
		global out_file3m
		global out_file3p
		global out_file3r
		global out_file3s
		# global out_file3t
		out_file3k = open(out_name + '.z_xconsensus.loc_all', "wb")
		out_file3la = open(out_name + '.z_xconsensus.conv_adjc', "wb")
		out_file3lb = open(out_name + '.z_xconsensus.conv_real', "wb")
		out_file3m = open(out_name + '.z_xconsensus.debug', "wb")
		out_file3p = open(out_name + '.z_xconsensus.matrix', "wb")
		out_file3r = open(out_name + '.z_xconsensus.frame', "wb")
		out_file3s = open(out_name + '.z_xconsensus.loc_nr', "wb")
		# out_file3t = open(out_name + '.z_xconsensus.dupl', "wb")
	if map_construction == "TRUE":
		# global out_file3k
		# global out_file3l
		# global out_file3m
		# global out_file3p
		# global out_file3r
		out_file3c = open(out_name + '.z_trio_all', "wb")
		out_file3d = open(out_name + '.z_trio_good', "wb")
		out_file3e1 = open(out_name + '.z_trio_graph1', "wb")
		out_file3e2 = open(out_name + '.z_trio_graph2', "wb")
		out_file3f = open(out_name + '.z_trio_map', "wb")
		out_file3g = open(out_name + '.z_trio_best', "wb")
		out_file3h = open(out_name + '.z_trio_bad', "wb")
		# out_file3i = open(out_name + '.z_marker_sum', "wb")
		# out_file3k = open(out_name + '.z_consensus.loc', "wb")
		# out_file3l = open(out_name + '.z_consensus.conv', "wb")
		# out_file3m = open(out_name + '.z_consensus.debug', "wb")
		# out_file3p = open(out_name + '.z_consensus.matrix', "wb")
		# out_file3r = open(out_name + '.z_consensus.frame', "wb")
	##############################################

	out_file6.write("=============================================" + '\n')
	out_file6.write("    RUN PARAMETERS: " + '\n')
	out_file6.write(" 1. INPUT  FILE:  " + in_name + '\n')
	out_file6.write(" 2. OUTPUT FILE:  " + out_name + '\n')
	out_file6.write(" 3. RECM CUTOFF:  " + str(rec_cut) + '\n')
	out_file6.write(" 4. BITS CUTOFF:  " + str(bit_cut) + '\n')
	out_file6.write(" 5. DATA CUTOFF:  " + str(dat_cut) + '\n')
	out_file6.write(" 6. ALLELE DIST:  " + str(allele_dist) + '\n')
	out_file6.write(" 7. MISSING DAT:  " + str(abs_loss) + '\n')
	out_file6.write(" 8. FRAME  LIST:  " + frame_file_name + '\n')
	out_file6.write(" 9. TRIO ANALYS:  " + map_construction + '\n')
	out_file6.write("10. DOUBLE LIMT:  " + str(double_limit) + '\n')

	time.sleep(2)

	global id_list
	global id_array
	global pairs_array
	global pairs_array_N
	global frame_array
	global tree_clust_array
	global graph_depth
	global marker_depth
	global sequence_array_bin
	global init_len
	global nr_good_list

	id_list = []		; # LIST OF ALL NON-REDUNDANT IDs
	id_array = {}		; # TO CHECK AND CREATE NON-REDUNDANT LIST
	matrix_array = {}	; # PAIRWISE DATA FOR GOOD PAIRS ( <= 0.4 )
	matrix_array_X = {}	; # PAIRWISE DATA FOR GOOD PAIRS ( <= 0.25 )
	matrix_array_1 = {}	; # 0.20
	matrix_array_2 = {}	; # 0.18
	matrix_array_3 = {}	; # 0.16
	matrix_array_4 = {}	; # 0.14
	matrix_array_5 = {}	; # 0.12
	matrix_array_6 = {}	; # 0.10
	matrix_array_7 = {}	; # 0.09
	matrix_array_8 = {}	; # 0.08
	matrix_array_9 = {}	; # 0.07
	matrix_array_A = {}	; # 0.06
	matrix_array_B = {}	; # 0.05
	matrix_array_C = {}	; # 0.04
	matrix_array_D = {}	; # 0.03
	matrix_array_E = {}	; # 0.02
	matrix_array_F = {}	; # 0.01
	matrix_array_G = {}	; # 0.00
	matrix_array_N = {}	; # PAIRWISE DATA FOR NEGATIVE LINKAGE ( >= 0.6 )
	matrix_array_N1 = {}	; # -0.7
	matrix_array_N2 = {}	; # -0.8
	matrix_array_N3 = {}	; # -0.9
	pairs_array  = {}	; # ARRAY OF ALL POSITIVE PAIRS IN THE MATRIX
	pair_counter = 0	; # COUNTER OF POSITIVE PAIRS IN THE MATRIX
	pairs_array_N = {}	; # ARRAY OF ALL NEGATIVE PAIRS IN THE MATRIX
	pair_counter_N = 0	; # COUNTER OF NEGATIVE PAIRS IN THE MATRIX
	sequence_array = {}	; # ARRAY OF ALL SEQUENCES
	sequence_array_bin = {}	; # ARRAY OF ALL SEQUENCES IN BINARY MODE: A or B and -
	frame_array = {}	; # ARRAY OF FRAME MARKERS
	frame_position = {}	; # ARRAY WITH FRAME MARKERS COORDINATES
	bit_sum_array = {}	; # SUMMARY OF BIT SCORES
	rec_sum_array = {}	; # SUMMARY OF REC SCORES
	bit_abs_array = {}	; # SUMMARY OF BIT SCORES
	rec_abs_array = {}	; # SUMMARY OF REC SCORES
	tree_clust_array = {}	; # SUMMARY OF CLUSTERING ARRAY

	nr_scores_array = {}	; # ARRAY NON-REDUNDANT MARKER SCORES I
	nr_markers_array = {}	; # ARRAY NON-REDUNDANT MARKER SCORES II
	all_scores_array = {}	; # ARRAY WITH ALL MARKER SCORES

	graph_depth  = {}	; # TYPE OF GRAPH: COMPLETE or CONNECTED
	marker_depth = {}	; # MARKER DEPTH: IS IT CONNECTED DIRECTLY TO ALL OTHER MARKERS IN A GROUP?

	count_A = {}
	count_B = {}
	count_C = {}
	count_D = {}
	count_H = {}
	count_X = {}

	count_nr_A = {}
	count_nr_B = {}
	count_nr_X = {}

	rec_list = ["00_01", "01_02", "02_03", "03_04", "04_05", "05_06", "06_07", "07_08", "08_09", "09_10"]
	bit_list = ["BIT_POS", "BIT_MED", "BIT_NEG"]

	if map_construction != "NR_SET":

		### HEADER IN SUMMARY FILE ###
		out_file7.write("MARKER_ID" + '\t')
		for item in rec_list:
			out_file7.write(item + '\t')
		out_file7.write("REC:P-N" + '\t')
		out_file7.write("REC_ABS" + '\t')
		out_file7.write("***" + '\t')
		for item in bit_list:
			out_file7.write(item + '\t')
		out_file7.write("BIT:P-N" + '\t')
		out_file7.write("BIT_ABS" + '\t')
		out_file7.write("***" + '\t')
		out_file7.write("LG_SUMMARY" + '\t')
		# out_file7.write("NP_SUMMARY" + '\n')

		out_file7.write("NP_SUMMARY" + '\t')
		out_file7.write("***" + '\t')
		out_file7.write("A" + '\t')
		out_file7.write("D" + '\t')
		out_file7.write("B" + '\t')
		out_file7.write("C" + '\t')
		out_file7.write("H" + '\t')
		out_file7.write("X" + '\t')
		out_file7.write("*" + '\t')
		out_file7.write("AD" + '\t')
		out_file7.write("BC" + '\t')
		out_file7.write("ALL" + '\n')

	global print_frame
	print_frame = "FALSE"
	### TRY TO READ FRAME MARKERS LIST ###
	### README_001 ###
	try:
		frame_file = open(frame_file_name, "rb")
		print "USING FRAME MARKERS LIST"
		while 1:
			u = frame_file.readline()
			if u == '':
				break
			if '\n' in u:
				u = u[:-1]
			if '\r' in u:
				u = u[:-1]
			u = u.split('\t')
			fm = u[1]
			fl = u[0]
			try:
				fp = u[2]
			except:
				fp = "-"
			frame_array[fm] = fl
			frame_position[fm] = fp
			print_frame = "TRUE"
	except:
		print "DID NOT FIND FRAME MARKERS FILE:  " + frame_file_name
		print "WORKING WITHOUT FRAME MARKERS LIST"
		# continue

	### READ DATA FILE ###
	print "============================================="
	time.sleep(2)
	out_file6.write("=============================================" + '\n')
	n = 0
	d = 0
	l = 1
	init_len = 10000
	duplic_id = []
	while 1:
		### README_001 ###
		t = in_file.readline()
		if t == '':
			break
		if '\n' in t:
			t = t[:-1]
		if '\r' in t:
			t = t[:-1]
		tl = t.split('\t')
		curr_len = len(tl)
		if tl[0][0] == ";":
			print "============================================="
			print tl
			print "JUNK LINE"
			print "============================================="
			out_file4.write(t + '\n')
			out_file3a.write(t + '\n')
			time.sleep(2)
		if l == 1 and curr_len >= 12 and tl[0][0] != ";":
			init_len = curr_len
		if curr_len == init_len and tl[0][0] != ";":
			####
			id = tl[0]
			scores_string = []
			# if id not in id_list:
			try:
				id_test = id_array[id]
				duplic_id.append(id)
				out_file5.write(t + '\n')
				print '\n'
				print id + "    IS DUPLICATED -=- CHECK DATA FOR DUPLICATION"
				out_file6.write( id + "    IS DUPLICATED -=- CHECK DATA FOR DUPLICATION" + '\n')
				d = d + 1
				time.sleep(2)
			except:
				id_array[id] = 1
				id_list.append(id)
				out_file4.write(t + '\n')
				count_A[id] = 0
				count_D[id] = 0
				count_B[id] = 0
				count_C[id] = 0
				count_H[id] = 0
				count_X[id] = 0
				q = 1

				## COLLECTING DATAPOINTS
				while q < init_len:
					data_point = tl[q]
					if data_point == "A" or data_point == "B" or data_point == "C" \
						or data_point == "D" or data_point == "H" or data_point == "-":

						### README_002 ###
						sequence_array[id,q] = data_point

						if data_point == "A":
							count_A[id] = count_A[id] + 1
						if data_point == "D":
							count_D[id] = count_D[id] + 1
						if data_point == "B":
							count_B[id] = count_B[id] + 1
						if data_point == "C":
							count_C[id] = count_C[id] + 1
						if data_point == "H":
							count_H[id] = count_H[id] + 1
						if data_point == "-":
							count_X[id] = count_X[id] + 1

						### IT WORKS ONLY FOR INBRED LINES ! ###
						### SCORES STRING
						score_point = data_point
						if score_point == "C":
							score_point = "B"
						if score_point == "D":
							score_point = "A"
						if score_point == "H":
							score_point = "-"

						scores_string.append(score_point)

						### README_003 ###
						sequence_array_bin[id,q] = score_point

					else:
						print '\n'
						print "WRONG DATA FORMAT"
						print "CHECK LINE:  " + `l` + '\n'
						print t
						print "============================================="
						print "VALUE:   " + data_point
						print "============================================="
						print "TERMINATED.........."
						print "============================================="
						out_file6.write("TERMINATED")
						out_file6.write(".........." + '\n')
						sys.exit()
					q = q + 1
				sys.stdout.write(".")
				all_scores_array[id] = scores_string
				n = n + 1
			l = l + 1
		if curr_len != init_len and l > 1:
			print '\n'
			print "WRONG NUMBER OF DATA POINTS"
			print "CHECK LINE:  " + `l` + '\n'
			print t
			print "============================================="
			print "TERMINATED.........."
			print "============================================="
			out_file6.write("TERMINATED")
			out_file6.write(".........." + '\n')
			sys.exit()

	print '\n'
	print "============================================="
	print `n` + " UNIQ IDs IN THE SET FOUND"
	print `d` + " IDs ARE DUPLICATED"
	out_file6.write("=============================================" + '\n')
	out_file6.write(`n` + " UNIQ IDs IN THE SET FOUND" + '\n')
	out_file6.write(`d` + " IDs ARE DUPLICATED" + '\n')

	duplic_id.sort()
	print duplic_id
	out_file6.write("=============================================" + '\n')

	print "CONTINUE ANALYSIS WITH " + `n` + " SEQUENCES OUT OF " + `n+d`
	out_file6.write("CONTINUE ANALYSIS WITH " + `n` + " SEQUENCES OUT OF " + `n+d` + '\n')
	print "============================================="
	print "SUCCESS!!!"
	print "============================================="
	out_file6.write("=============================================" + '\n')

	if n >= 5001:
		print_all_pairs = "FALSE"
		print "TOO MANY MARKERS"
		print "ALL PAIRS FILE WILL NOT BE GENERATED"

	time.sleep(2)

	######### PAIRWISE COMPARISON #########

	## SORTED ID LIST ##
	id_list.sort()
	## IS IT BETTER NOT TO SORT ? ##

	### SET SUMMARY VALUES FOR MARKERS ####

	nr_scores_list = []
	nr_marker_list = []
	nr_good_list   = []
	really_bad     = []	; # ALL BAD MARKERS
	really_bad_allele = []	; # BAD WITH ALLELE DISTORTION
	really_bad_data   = []	; # BAD WITH LOSS OF DATA
	r_count = 0
	dupl_markers_list = {}
	good_string = {}
	ratio_AB = {}
	bias_AB = {}

	for item in id_list:
		tree_clust_array[item] = []
		rec_abs_array[item] = 0
		bit_abs_array[item] = 0
		for rec_range in rec_list:
			rec_sum_array[item,rec_range] = 0
		for bit_range in bit_list:
			bit_sum_array[item,bit_range] = 0

		### MARKER SCORE DUPLICATION ###
		score_string = "\t".join(all_scores_array[item])
		text_string  = "".join(all_scores_array[item])
		count_nr_A[item] = text_string.count("A")
		count_nr_B[item] = text_string.count("B")
		count_nr_X[item] = text_string.count("-")
		good_string[item] = "BAD"
		ratio_AB[item] = -1
		bias_AB[item] = "XXX"
		if count_nr_A[item] == 0:
			good_string[item] = "BAD"
		if count_nr_B[item] == 0:
			good_string[item] = "BAD"
		if count_nr_A[item] > 0 and count_nr_B[item] > 0:
			if count_nr_A[item] == count_nr_B[item]:
				ratio_AB[item] = count_nr_A[item]*1.0/count_nr_B[item]
				bias_AB[item] = "A=B"
			if count_nr_A[item] > count_nr_B[item]:
				ratio_AB[item] = count_nr_B[item]*1.0/count_nr_A[item]
				bias_AB[item] = "A>B"
			if count_nr_A[item] < count_nr_B[item]:
				ratio_AB[item] = count_nr_A[item]*1.0/count_nr_B[item]
				bias_AB[item] = "A<B"
		### README_005 ###
		if ratio_AB[item] >= allele_dist and count_nr_X[item] <= abs_loss:
			good_string[item] = "GOOD"
		if good_string[item] == "BAD":
			really_bad.append(item)
			if ratio_AB[item] < allele_dist:
				really_bad_allele.append(item)
			if count_nr_X[item] > abs_loss:
				really_bad_data.append(item)
		if score_string in nr_scores_list:
			print "DATA DUPLICATION FOR MARKER:  " + item
			# time.sleep(1)
		### README_004 ###
		if score_string not in nr_scores_list:
			nr_scores_list.append(score_string)
			nr_marker_list.append(item)
			if good_string[item] == "GOOD":
				### README_007 ###
				nr_good_list.append(item)
				out_file3a.write(item + '\t' + score_string + '\n')
			nr_scores_array[item] = all_scores_array[item]
			nr_markers_array[score_string] = item

	for item_x in id_list:
		# dupl_markers_list[item_x] = []
		for item_y in id_list:
			if all_scores_array[item_x] == all_scores_array[item_y] and item_x != item_y:
				# print "IDENTICAL DATA FOR PAIR: " + item_x + "  " + item_y
				working_id = "++++++"
				# if item_x in nr_marker_list:
				#	working_id = item_x
				# if item_y in nr_marker_list:
				#	working_id = item_y
				score_string = "\t".join(all_scores_array[item_x])
				try:
					working_id = nr_markers_array[score_string]
					# dupl_markers_list[item_y] = working_id
					dupl_markers_list[item_x] = working_id
				except:
					working_id = "------"
				### README_004 ###
				print "IDENTICAL DATA FOR PAIR: " + item_x + "  " + item_y + "    WORKING ID: " + working_id
				out_file3b.write(item_x + '\t' + "==" + '\t' + item_y + '\t' + "***" + '\t' + working_id + '\n')
				# dupl_markers_list[item_y] = working_id
				# dupl_markers_list[item_x] = working_id
				


	if map_construction == "NR_SET":
		print "                                 "
		print "  EXIT HERE                      "
		print "  RUN MadMapper WITH NR SET NOW  "
		print "                                 "
		sys.exit()

	### START ###

	dummy_cat = 1
	dummy_len = len(id_list)
	dummy_value = dummy_len*dummy_len
	print `dummy_value` + "  PAIRWISE COMPARISONS HAVE TO BE DONE"

	### README_008 ###
	for item_a in id_list:
		for item_b in id_list:
			p = 1
			bit_score = 0
			recomb = 0
			data_loss = 0
			data_points = 0

			line_tick_update = math.fmod(dummy_cat, 1000)
			if line_tick_update == 0:
				print `dummy_cat` + "  COMPARISONS DONE OUT OF  " + `dummy_value`

			while p < init_len:

				score_a = sequence_array[item_a,p]
				score_b = sequence_array[item_b,p]
				####  SCORES FOR PAIRWISE COMPARISON ####
				####    NO DATA    ####
				if score_a == "-" or score_b == "-":
					data_loss = data_loss + 1
				#### PERFECT MATCH ####
				if score_a == "A" and score_b == "A":
					recomb = recomb + rec_score_AA
					bit_score  = bit_score + bit_score_AA
					data_points = data_points + 1
				if score_a == "B" and score_b == "B":
					recomb = recomb + rec_score_BB
					bit_score  = bit_score + bit_score_BB
					data_points = data_points + 1
				####   REPULSION   ####
				if score_a == "A" and score_b == "B":
					recomb = recomb + rec_score_AB
					bit_score  = bit_score + bit_score_AB
					data_points = data_points + 1
				if score_a == "B" and score_b == "A":
					recomb = recomb + rec_score_BA
					bit_score  = bit_score + bit_score_BA
					data_points = data_points + 1
				#### MONO - HETERO ####
				if score_a == "A" and score_b == "H":
					recomb = recomb + rec_score_AH
					bit_score  = bit_score + bit_score_AH
					data_points = data_points + 1
				if score_a == "B" and score_b == "H":
					recomb = recomb + rec_score_BH
					bit_score  = bit_score + bit_score_BH
					data_points = data_points + 1
				#### HETERO - MONO ####
				if score_a == "H" and score_b == "A":
					recomb = recomb + rec_score_HA
					bit_score  = bit_score + bit_score_HA
					data_points = data_points + 1
				if score_a == "H" and score_b == "B":
					recomb = recomb + rec_score_HB
					bit_score  = bit_score + bit_score_HB
					data_points = data_points + 1
				### HETERO - HETERO ###
				if score_a == "H" and score_b == "H":
					recomb = recomb + rec_score_HH
					bit_score  = bit_score + bit_score_HH
					data_points = data_points + 1
				####  C CASE ALL  ####
				if score_a == "C" and score_b == "C":
					recomb = recomb + rec_score_CC
					bit_score  = bit_score + bit_score_CC
					data_points = data_points + 1
				#### C - CASE (L)  ####
				if score_a == "C" and score_b == "A":
					recomb = recomb + rec_score_CA
					bit_score  = bit_score + bit_score_CA
					data_points = data_points + 1
				if score_a == "C" and score_b == "B":
					recomb = recomb + rec_score_CB
					bit_score  = bit_score + bit_score_CB
					data_points = data_points + 1
				if score_a == "C" and score_b == "H":
					recomb = recomb + rec_score_CH
					bit_score  = bit_score + bit_score_CH
					data_points = data_points + 1
				#### C - CASE (R)  ####
				if score_a == "A" and score_b == "C":
					recomb = recomb + rec_score_AC
					bit_score  = bit_score + bit_score_AC
					data_points = data_points + 1
				if score_a == "B" and score_b == "C":
					recomb = recomb + rec_score_BC
					bit_score  = bit_score + bit_score_BC
					data_points = data_points + 1
				if score_a == "H" and score_b == "C":
					recomb = recomb + rec_score_HC
					bit_score  = bit_score + bit_score_HC
					data_points = data_points + 1
				####  D CASE ALL  ####
				if score_a == "D" and score_b == "D":
					recomb = recomb + rec_score_DD
					bit_score  = bit_score + bit_score_DD
					data_points = data_points + 1
				#### D - CASE (L)  ####
				if score_a == "D" and score_b == "A":
					recomb = recomb + rec_score_DA
					bit_score  = bit_score + bit_score_DA
					data_points = data_points + 1
				if score_a == "D" and score_b == "B":
					recomb = recomb + rec_score_DB
					bit_score  = bit_score + bit_score_DB
					data_points = data_points + 1
				if score_a == "D" and score_b == "H":
					recomb = recomb + rec_score_DH
					bit_score  = bit_score + bit_score_DH
					data_points = data_points + 1
				#### D - CASE (R)  ####
				if score_a == "A" and score_b == "D":
					recomb = recomb + rec_score_AD
					bit_score  = bit_score + bit_score_AD
					data_points = data_points + 1
				if score_a == "B" and score_b == "D":
					recomb = recomb + rec_score_BD
					bit_score  = bit_score + bit_score_BD
					data_points = data_points + 1
				if score_a == "H" and score_b == "D":
					recomb = recomb + rec_score_HD
					bit_score  = bit_score + bit_score_HD
					data_points = data_points + 1
				#####################################
				####    D -=- C    ####
				if score_a == "D" and score_b == "C":
					recomb = recomb + rec_score_DC
					bit_score  = bit_score + bit_score_DC
					data_points = data_points + 1
				####    C -=- D    ####
				if score_a == "C" and score_b == "D":
					recomb = recomb + rec_score_CD
					bit_score  = bit_score + bit_score_CD
					data_points = data_points + 1
				#####################################

				p = p + 1

			dummy_cat = dummy_cat + 1
			### DATA POINT TEST ###
			if data_points + data_loss != init_len - 1:
				print "...SOMETHING IS WRONG..."
				sys.exit()

			### README_009 ###
			### ANALYZE ONLY VALUABLE CRAP ###
			data_fr = data_points*1.00/(data_points + data_loss)
			data_fr = round(data_fr,round_scale)
			data_fr_str = str(data_fr)
			# if data_fr >= dat_cut:
			if data_points >= dat_cut:

				rec_dec = recomb*1.00/data_points
				# rec_dec = round(rec_dec,2)
				rec_dec = round(rec_dec,round_scale)
				rec_str = str(rec_dec)
				# data_fr = round(data_fr,2)
				# data_fr = round(data_fr,round_scale)
				# data_fr_str = str(data_fr)

				### SUMMARY DATA ###

				rec_abs_array[item_a] = rec_abs_array[item_a] + (0.5 - rec_dec)
				bit_abs_array[item_a] = bit_abs_array[item_a] + bit_score

				rec_data = rec_dec

				### README_013 ###
				if rec_data <= 0.1 and rec_data >= 0:
					rec_sum_array[item_a,"00_01"] = rec_sum_array[item_a,"00_01"] + 1
				if rec_data <= 0.2 and rec_data > 0.1:
					rec_sum_array[item_a,"01_02"] = rec_sum_array[item_a,"01_02"] + 1
				if rec_data <= 0.3 and rec_data > 0.2:
					rec_sum_array[item_a,"02_03"] = rec_sum_array[item_a,"02_03"] + 1
				if rec_data <= 0.4 and rec_data > 0.3:
					rec_sum_array[item_a,"03_04"] = rec_sum_array[item_a,"03_04"] + 1
				if rec_data <= 0.5 and rec_data > 0.4:
					rec_sum_array[item_a,"04_05"] = rec_sum_array[item_a,"04_05"] + 1
				if rec_data <= 0.6 and rec_data > 0.5:
					rec_sum_array[item_a,"05_06"] = rec_sum_array[item_a,"05_06"] + 1
				if rec_data <= 0.7 and rec_data > 0.6:
					rec_sum_array[item_a,"06_07"] = rec_sum_array[item_a,"06_07"] + 1
				if rec_data <= 0.8 and rec_data > 0.7:
					rec_sum_array[item_a,"07_08"] = rec_sum_array[item_a,"07_08"] + 1
				if rec_data <= 0.9 and rec_data > 0.8:
					rec_sum_array[item_a,"08_09"] = rec_sum_array[item_a,"08_09"] + 1
				if rec_data <= 1.0 and rec_data > 0.9:
					rec_sum_array[item_a,"09_10"] = rec_sum_array[item_a,"09_10"] + 1

				if bit_score >= 100:
					bit_sum_array[item_a,"BIT_POS"] = bit_sum_array[item_a,"BIT_POS"] + 1
				if bit_score <= -100:
					bit_sum_array[item_a,"BIT_NEG"] = bit_sum_array[item_a,"BIT_NEG"] + 1
				if bit_score > -100 and bit_score < 100:
					bit_sum_array[item_a,"BIT_MED"] = bit_sum_array[item_a,"BIT_MED"] + 1

				### MATRIX DATA  ###
				### README_012 ###
				if rec_dec >= 0.6:

					pair_counter_N = pair_counter_N + 1
					current_pair = [item_a, item_b]
					pairs_array_N[pair_counter_N] = current_pair
					matrix_values = [rec_dec, bit_score, data_fr, data_points]
					matrix_array_N[item_a,item_b] = matrix_values
					matrix_array_N1[item_a,item_b] = matrix_values
					matrix_array_N2[item_a,item_b] = matrix_values
					matrix_array_N3[item_a,item_b] = matrix_values

					out_file2.write(item_a + '\t' + item_b + '\t' + rec_str + '\t' + \
						`bit_score` + '\t' + data_fr_str + '\t' + "***" + '\t' + `recomb` + \
						'\t' + `data_points` + '\t' + `data_loss` + '\t' + `data_points + data_loss` + '\n')

				### README_011 ###
				if rec_dec <= 0.4:

					pair_counter = pair_counter + 1
					current_pair = [item_a, item_b]
					pairs_array[pair_counter] = current_pair
					matrix_values = [rec_dec, bit_score, data_fr, data_points]
					matrix_array[item_a,item_b] = matrix_values
					### if rec_dec <= 0.25:
					if rec_dec <= 0.25 and map_construction == "TRUE":
						matrix_array_X[item_a,item_b] = matrix_values	; # 0.25
					if rec_dec <= rec_cut:
						matrix_array_1[item_a,item_b] = matrix_values	; # 0.20
					if rec_dec <= rec_cut - 0.02:
						matrix_array_2[item_a,item_b] = matrix_values	; # 0.18
					if rec_dec <= rec_cut - 0.04:
						matrix_array_3[item_a,item_b] = matrix_values	; # 0.16
					if rec_dec <= rec_cut - 0.06:
						matrix_array_4[item_a,item_b] = matrix_values	; # 0.14
					if rec_dec <= rec_cut - 0.08:
						matrix_array_5[item_a,item_b] = matrix_values	; # 0.12
					if rec_dec <= rec_cut - 0.10:
						matrix_array_6[item_a,item_b] = matrix_values	; # 0.10
					if rec_dec <= rec_cut - 0.11:
						matrix_array_7[item_a,item_b] = matrix_values	; # 0.09
					if rec_dec <= rec_cut - 0.12:
						matrix_array_8[item_a,item_b] = matrix_values	; # 0.08
					if rec_dec <= rec_cut - 0.13:
						matrix_array_9[item_a,item_b] = matrix_values	; # 0.07
					if rec_dec <= rec_cut - 0.14:
						matrix_array_A[item_a,item_b] = matrix_values	; # 0.06
					if rec_dec <= rec_cut - 0.15:
						matrix_array_B[item_a,item_b] = matrix_values	; # 0.05
					if rec_dec <= rec_cut - 0.16:
						matrix_array_C[item_a,item_b] = matrix_values	; # 0.04
					if rec_dec <= rec_cut - 0.17:
						matrix_array_D[item_a,item_b] = matrix_values	; # 0.03
					if rec_dec <= rec_cut - 0.18:
						matrix_array_E[item_a,item_b] = matrix_values	; # 0.02
					if rec_dec <= rec_cut - 0.19:
						matrix_array_F[item_a,item_b] = matrix_values	; # 0.01
					if rec_dec <= rec_cut - 0.20:
						matrix_array_G[item_a,item_b] = matrix_values	; # 0.00

					out_file0.write(item_a + '\t' + item_b + '\t' + rec_str + '\t' + \
						`bit_score` + '\t' + data_fr_str + '\t' + "***" + '\t' + `recomb` + \
						'\t' + `data_points` + '\t' + `data_loss` + '\t' + \
						`data_points + data_loss` + '\n')

				### README_010 ###
				if print_all_pairs == "TRUE":
					out_file1.write(        item_a + '\t' + item_b + '\t' + rec_str + '\t' + \
							`bit_score` + '\t' + data_fr_str + '\t' + "***" + '\t' + `recomb` + \
							'\t' + `data_points` + '\t' + `data_loss` + '\t' + \
							`data_points + data_loss` + '\n')
				########################################
			# if data_fr < dat_cut:
			if data_points < dat_cut:
				print "DATA POINTS ARE LESS THAN CUTOFF FOR PAIR: " + item_a + " " + item_b + " " + \
										`data_points` + "   " + data_fr_str
			########################################
		### SUMMARY  FILE ###
		### README_013 ###
		rec_sum_diff1 = rec_sum_array[item_a,"00_01"] + rec_sum_array[item_a,"01_02"] \
			- rec_sum_array[item_a,"08_09"] - rec_sum_array[item_a,"09_10"]
		bit_sum_diff1 = bit_sum_array[item_a,"BIT_POS"] - bit_sum_array[item_a,"BIT_NEG"]
		rec_sum_diff2 = rec_abs_array[item_a]
		bit_sum_diff2 = bit_abs_array[item_a]
		sum_text = "---"
		lnkg_sum = "CLASS_X"
		### NP SUMMARY ###
		if rec_sum_diff2 > 0 and bit_sum_diff2 >  0:
			sum_text = "POSITIVE"
		if rec_sum_diff2 ==0 and bit_sum_diff2 == 0:
			sum_text = "NEUTRAL"
		if rec_sum_diff2 < 0 and bit_sum_diff2 >= 0:
			sum_text = "NEGATIVE REC"
		if bit_sum_diff2 < 0 and rec_sum_diff2 >= 0:
			sum_text = "NEGATIVE BIT"
		if bit_sum_diff2 < 0 and rec_sum_diff2 <  0:
			sum_text = "NEGATIVE ALL"
		### LINKAGE SUMMARY ###
		### README_13LG ###
		min_best = 6
		min_good = 3
		min_step = 3
		### CLASS E
		if rec_sum_array[item_a,"00_01"] >= 0 and rec_sum_array[item_a,"00_01"] < min_good:
			if rec_sum_array[item_a,"01_02"] >= rec_sum_array[item_a,"00_01"] and rec_sum_array[item_a,"02_03"] >= rec_sum_array[item_a,"01_02"]:
				lnkg_sum = "CLASS_E"
		### CLASS C AND D
		if rec_sum_array[item_a,"00_01"] >= min_good and rec_sum_array[item_a,"00_01"] < min_best:
			if rec_sum_array[item_a,"01_02"] >= rec_sum_array[item_a,"00_01"] and rec_sum_array[item_a,"02_03"] >= rec_sum_array[item_a,"01_02"]:
				lnkg_sum = "CLASS_D"
			if rec_sum_array[item_a,"01_02"] >= (rec_sum_array[item_a,"00_01"]+min_step) and rec_sum_array[item_a,"02_03"] >= (rec_sum_array[item_a,"01_02"]+min_step):
				lnkg_sum = "CLASS_C"
		### CLASS A AND B
		if rec_sum_array[item_a,"00_01"] >= min_best:
			if rec_sum_array[item_a,"01_02"] >= rec_sum_array[item_a,"00_01"] and rec_sum_array[item_a,"02_03"] >= rec_sum_array[item_a,"01_02"]:
				lnkg_sum = "CLASS_B"
			if rec_sum_array[item_a,"01_02"] >= (rec_sum_array[item_a,"00_01"]+min_step) and rec_sum_array[item_a,"02_03"] >= (rec_sum_array[item_a,"01_02"]+min_step):
				lnkg_sum = "CLASS_A"
		### CLASS F
		if rec_sum_array[item_a,"00_01"] == 0 and rec_sum_array[item_a,"01_02"] == 0:
			lnkg_sum = "CLASS_F"
		#######################
		### README_013 ###
		count_ABCDHX = count_A[item_a]+count_D[item_a]+count_B[item_a]+count_C[item_a]+count_H[item_a]+count_X[item_a]
		out_file7.write(                item_a + '\t' + str(rec_sum_array[item_a,"00_01"]) + '\t' + \
						str(rec_sum_array[item_a,"01_02"]) + '\t' + \
						str(rec_sum_array[item_a,"02_03"]) + '\t' + \
						str(rec_sum_array[item_a,"03_04"]) + '\t' + \
						str(rec_sum_array[item_a,"04_05"]) + '\t' + \
						str(rec_sum_array[item_a,"05_06"]) + '\t' + \
						str(rec_sum_array[item_a,"06_07"]) + '\t' + \
						str(rec_sum_array[item_a,"07_08"]) + '\t' + \
						str(rec_sum_array[item_a,"08_09"]) + '\t' + \
						str(rec_sum_array[item_a,"09_10"]) + '\t' + \
						str(rec_sum_diff1) + '\t' + \
						str(round(rec_sum_diff2,round_scale)) + '\t' + \
						"***" + '\t' + \
						str(bit_sum_array[item_a,"BIT_POS"]) + '\t' + \
						str(bit_sum_array[item_a,"BIT_MED"]) + '\t' + \
						str(bit_sum_array[item_a,"BIT_NEG"]) + '\t' + \
						str(bit_sum_diff1) + '\t' + \
						str(bit_sum_diff2) + '\t' + \
						"***" + '\t' + \
						lnkg_sum + '\t' + \
						sum_text + '\t' + "***" + '\t' + `count_A[item_a]` + '\t' + \
						`count_D[item_a]` + '\t' + `count_B[item_a]` + '\t' + `count_C[item_a]` + '\t' + \
						`count_H[item_a]` + '\t' + `count_X[item_a]` + '\t' "*" + '\t' + \
						`count_A[item_a]+count_D[item_a]` + '\t' + `count_B[item_a]+count_C[item_a]` + '\t' + \
						`count_ABCDHX` + '\t' + '\n')
						### END OF MARKER SUMMARY FILE ###
		#####################
		sys.stdout.write(".")
		#####################
	print ""
	print "============================================="
	print "GLOBAL MATRIX CREATED"
	print "============================================="

	#### POSITIVE CLUSTERING ####

	### README_014 ###

	print ""
	print "STARTING POSITIVE CLUSTERING"
	print "============================================="

	global out_file14

	out_file14 = open(out_name + '.group_info' + "_Summary", "wb")

	time.sleep(2)

	Seqs_Clustering(rec_cut, bit_cut, dat_cut, "01", matrix_array_1)	; # 0.20
	print "ROUND 1 DONE"
	print `group_count` + "  GROUPS WERE FOUND (STEP 01)"
	out_file6.write(`group_count` + "  GROUPS WERE FOUND (STEP 01)" + '\n')

	time.sleep(2)

	Seqs_Clustering(rec_cut, bit_cut, dat_cut, "02", matrix_array_2)	; # 0.18
	print "ROUND 2 DONE"
	print `group_count` + "  GROUPS WERE FOUND (STEP 02)"
	out_file6.write(`group_count` + "  GROUPS WERE FOUND (STEP 02)" + '\n')

	time.sleep(2)

	Seqs_Clustering(rec_cut, bit_cut, dat_cut, "03", matrix_array_3)	; # 0.16
	print "ROUND 3 DONE"
	print `group_count` + "  GROUPS WERE FOUND (STEP 03)"
	out_file6.write(`group_count` + "  GROUPS WERE FOUND (STEP 03)" + '\n')

	time.sleep(2)

	Seqs_Clustering(rec_cut, bit_cut, dat_cut, "04", matrix_array_4)	; # 0.14
	print "ROUND 4 DONE"
	print `group_count` + "  GROUPS WERE FOUND (STEP 04)"
	out_file6.write(`group_count` + "  GROUPS WERE FOUND (STEP 04)" + '\n')

	time.sleep(2)

	Seqs_Clustering(rec_cut, bit_cut, dat_cut, "05", matrix_array_5)	; # 0.12
	print "ROUND 5 DONE"
	print `group_count` + "  GROUPS WERE FOUND (STEP 05)"
	out_file6.write(`group_count` + "  GROUPS WERE FOUND (STEP 05)" + '\n')

	time.sleep(2)

	Seqs_Clustering(rec_cut, bit_cut, dat_cut, "06", matrix_array_6)	; # 0.10
	print "ROUND 6 DONE"
	print `group_count` + "  GROUPS WERE FOUND (STEP 06)"
	out_file6.write(`group_count` + "  GROUPS WERE FOUND (STEP 06)" + '\n')

	time.sleep(2)

	if deep_clustering == "TRUE":

		Seqs_Clustering(rec_cut, bit_cut, dat_cut, "07", matrix_array_7)	; # 0.09
		print "ROUND 7 DONE"
		print `group_count` + "  GROUPS WERE FOUND (STEP 07)"
		out_file6.write(`group_count` + "  GROUPS WERE FOUND (STEP 07)" + '\n')

		time.sleep(2)

		Seqs_Clustering(rec_cut, bit_cut, dat_cut, "08", matrix_array_8)	; # 0.08
		print "ROUND 8 DONE"
		print `group_count` + "  GROUPS WERE FOUND (STEP 08)"
		out_file6.write(`group_count` + "  GROUPS WERE FOUND (STEP 08)" + '\n')

		time.sleep(2)

		Seqs_Clustering(rec_cut, bit_cut, dat_cut, "09", matrix_array_9)	; # 0.07
		print "ROUND 9 DONE"
		print `group_count` + "  GROUPS WERE FOUND (STEP 09)"
		out_file6.write(`group_count` + "  GROUPS WERE FOUND (STEP 09)" + '\n')

		time.sleep(2)

		Seqs_Clustering(rec_cut, bit_cut, dat_cut, "10", matrix_array_A)	; # 0.06
		print "ROUND 10 DONE"
		print `group_count` + "  GROUPS WERE FOUND (STEP 10)"
		out_file6.write(`group_count` + "  GROUPS WERE FOUND (STEP 10)" + '\n')

		time.sleep(2)

		Seqs_Clustering(rec_cut, bit_cut, dat_cut, "11", matrix_array_B)	; # 0.05
		print "ROUND 11 DONE"
		print `group_count` + "  GROUPS WERE FOUND (STEP 11)"
		out_file6.write(`group_count` + "  GROUPS WERE FOUND (STEP 11)" + '\n')

		time.sleep(2)

		Seqs_Clustering(rec_cut, bit_cut, dat_cut, "12", matrix_array_C)	; # 0.04
		print "ROUND 12 DONE"
		print `group_count` + "  GROUPS WERE FOUND (STEP 12)"
		out_file6.write(`group_count` + "  GROUPS WERE FOUND (STEP 12)" + '\n')

		time.sleep(2)

		Seqs_Clustering(rec_cut, bit_cut, dat_cut, "13", matrix_array_D)	; # 0.03
		print "ROUND 13 DONE"
		print `group_count` + "  GROUPS WERE FOUND (STEP 13)"
		out_file6.write(`group_count` + "  GROUPS WERE FOUND (STEP 13)" + '\n')

		time.sleep(2)

		Seqs_Clustering(rec_cut, bit_cut, dat_cut, "14", matrix_array_E)	; # 0.02
		print "ROUND 14 DONE"
		print `group_count` + "  GROUPS WERE FOUND (STEP 14)"
		out_file6.write(`group_count` + "  GROUPS WERE FOUND (STEP 14)" + '\n')

		time.sleep(2)

		Seqs_Clustering(rec_cut, bit_cut, dat_cut, "15", matrix_array_F)	; # 0.01
		print "ROUND 15 DONE"
		print `group_count` + "  GROUPS WERE FOUND (STEP 15)"
		out_file6.write(`group_count` + "  GROUPS WERE FOUND (STEP 15)" + '\n')

		time.sleep(2)

		Seqs_Clustering(rec_cut, bit_cut, dat_cut, "16", matrix_array_G)	; # 0.00
		print "ROUND 16 DONE"
		print `group_count` + "  GROUPS WERE FOUND (STEP 16)"
		out_file6.write(`group_count` + "  GROUPS WERE FOUND (STEP 16)" + '\n')

		time.sleep(2)

	#### NEGATIVE CLUSTERING ####

	print ""
	print "============================================="
	print "STARTING NEGATIVE CLUSTERING"
	print "============================================="

	time.sleep(2)

	Seqs_Neg_Clustering(rec_cut, bit_cut, dat_cut, "N1", matrix_array_N1)	; # -0.7
	print "ROUND N1 DONE"
	print `group_count` + "  GROUPS WERE FOUND (STEP N1)"
	out_file6.write(`group_count` + "  GROUPS WERE FOUND (STEP N1)" + '\n')

	time.sleep(2)

	Seqs_Neg_Clustering(rec_cut, bit_cut, dat_cut, "N2", matrix_array_N2)	; # -0.8
	print "ROUND N2 DONE"
	print `group_count` + "  GROUPS WERE FOUND (STEP N2)"
	out_file6.write(`group_count` + "  GROUPS WERE FOUND (STEP N2)" + '\n')

	time.sleep(2)

	Seqs_Neg_Clustering(rec_cut, bit_cut, dat_cut, "N3", matrix_array_N3)	; # -0.9
	print "ROUND N3 DONE"
	print `group_count` + "  GROUPS WERE FOUND (STEP N3)"
	out_file6.write(`group_count` + "  GROUPS WERE FOUND (STEP N3)" + '\n')

	time.sleep(2)

	out_file6.write("=============================================" + '\n')

	### GROUP LIST ###
	f = 0
	for item in id_list:
		### TRY FRAME MARKER ID ###
		if print_frame == "TRUE":
			try:
				frame_group = frame_array[item]
			except:
				frame_group = "-"
			try:
				frame_coords = frame_position[item]
			except:
				frame_coords = "-"
		if print_frame == "FALSE":
			frame_group  = "-"
			frame_coords = "-"
		tree_clust_array[item].append("***")
		tree_clust_array[item].append(frame_group)
		tree_clust_array[item].append(frame_coords)
		tree_clust_array[item].append("***")
		tree_clust_array[item].append(str(f))
		tree_clust_array[item].append("***")
		tree_clust_array[item].append("LG")
		tree_clust_array[item].append(item)
		tree_clust_array[item].append("*")
		tree_clust_array[item].append(str(count_A[item]+count_D[item]))
		tree_clust_array[item].append(str(count_B[item]+count_C[item]))
		tree_clust_array[item].append(str(count_A[item]+count_D[item]+count_B[item]+count_C[item]+count_H[item]+count_X[item]))
		tree_clust_array[item].append("*")
		f = f + 1
		print tree_clust_array[item]

	print ""
	print "============================================="	
	print "            FINAL STEP                       "
	print "       PROCESSING MEGALIST                   "
	print "         DENDRO SORTING                      "
	print "============================================="
	print ""

	time.sleep(2)

	### MEGA LIST ###
	mega_list = []
	for item in id_list:
		mega_list.append(tree_clust_array[item])
	# print mega_list
	mega_list.sort()
	f = 0
	for item in mega_list:
		current_list = []
		item.append(str(f))
		for subitem in item:
			subitem = str(subitem)
			current_list.append(subitem)
		current_list = "\t".join(current_list)
		print current_list
		out_file8.write(current_list + '\n')
		f = f + 1

	#### MARKER SUMMARY ####
	acceptable_markers = []
	not_in_list_markers  = []
	very_bad_markers   = []
	core_list = []

	if map_construction == "TRUE":
		#### DOUBLE RECOMBINATION FOR TRIPLETS ###

		double_recomb_trio = {}	; # THREE DIM ARRAW WITH DOUBLE CROSS DATA
		best_recomb_trio = {}
		min_double_recomb = {}
		min_flank_recomb = {}
		dk_array = {}

		acceptable_markers = []
		not_in_list_markers  = []
		very_bad_markers   = []
		core_list = []

		nr_marker_list.sort()

		dummy_cat = 1
		# dummy_len = len(nr_marker_list)
		dummy_len = len(nr_good_list)
		dummy_value = dummy_len*dummy_len*dummy_len
		print "========================================================="
		print "GOOD MARKERS:  "
		print len(nr_good_list)
		print nr_good_list
		time.sleep(2)
		print "========================================================="
		print `dummy_value` + "  TRIPLET COMPARISONS HAVE TO BE DONE"
		print "========================================================="
		time.sleep(2)
		print "BAD MARKERS:  "
		print len(really_bad)
		print really_bad
		print "========================================================="
		time.sleep(2)

		# for item_c in nr_marker_list:
		### README_016 ###
		for item_c in nr_good_list:
			dk_array[item_c] = 0
			min_double_recomb[item_c] = 999
			# for item_u in nr_marker_list:
			for item_u in nr_good_list:
				# for item_d in nr_marker_list:
				for item_d in nr_good_list:
					mooba = "GOOD"
					booba = "GOOD"
					try:
						test_up = matrix_array_X[item_c,item_u]
						test_dn = matrix_array_X[item_c,item_d]
						# test_fl = matrix_array_X[item_u,item_d]
						test_fl = matrix_array[item_u,item_d]
						rec_cu = test_up[0]
						bit_cu = test_up[1]
						dat_cu_fr = test_up[2]
						dat_cu = test_up[3]
						rec_dn = test_dn[0]
						bit_dn = test_dn[1]
						dat_dn_fr = test_dn[2]
						dat_dn = test_dn[3]
						rec_fl = test_fl[0]
						bit_fl = test_fl[1]
						dat_fl_fr = test_fl[2]
						dat_fl = test_fl[3]
						if rec_cu > rec_cut:
							mooba = "BAD"
						if rec_dn > rec_cut:
							mooba = "BAD"
						if bit_cu < bit_cut:
							mooba = "BAD"
						if bit_dn < bit_cut:
							mooba = "BAD"
						if dat_cu < dat_cut:
							mooba = "BAD"
						if dat_dn < dat_cut:
							mooba = "BAD"
							
					except:
						mooba = "BAD"
						booba = "BAD"

					if mooba == "GOOD" and item_c != item_u and item_c != item_d and item_u != item_d:
						p = 1
						double_recomb = 0

						line_tick_update = math.fmod(dummy_cat, 1000)
						if line_tick_update == 0:
							print `dummy_cat` + "  COMPARISONS DONE OUT OF  " + `dummy_value`

						while p < init_len:

							score_u = sequence_array_bin[item_u,p]
							score_c = sequence_array_bin[item_c,p]
							score_d = sequence_array_bin[item_d,p]

							if score_u == "A" and score_c == "B" and score_d == "A":
								double_recomb = double_recomb + 1
							if score_u == "B" and score_c == "A" and score_d == "B":
								double_recomb = double_recomb + 1

							p = p + 1

						### README_017 ###
						if double_recomb <= min_double_recomb[item_c]:
							dk_array[item_c] = dk_array[item_c] + 1
							dk = dk_array[item_c]
							min_double_recomb[item_c] = double_recomb
							best_recomb_trio[item_c,dk] = item_u + '\t' + str(rec_cu) + '\t' + \
								str(bit_cu) + '\t' + str(dat_cu_fr) + '\t' + item_c + '\t' + \
								str(rec_dn) + '\t' + str(bit_dn) + '\t' + str(dat_dn_fr) + '\t' + \
								item_d + '\t' + "***" + '\t' + str(rec_fl) + '\t' + str(bit_fl) + '\t' + \
								str(dat_fl_fr) + '\t' + "***" + '\t' + str(double_recomb)

						double_recomb_trio[item_u,item_c,item_d] = double_recomb
						out_file3c.write(item_u + '\t' + str(rec_cu) + '\t' + str(bit_cu) + '\t' + \
							str(dat_cu_fr) + '\t' + item_c + '\t' + str(rec_dn) + '\t' + str(bit_dn) + '\t' + \
							str(dat_dn_fr) + '\t' + item_d + '\t' + "***" + '\t' + str(rec_fl) + '\t' + \
							str(bit_fl) + '\t' + str(dat_fl_fr) + '\t' + "***" + '\t' + \
							str(double_recomb) + '\n')

					### NUMBER OF COMPARISONS
					dummy_cat = dummy_cat + 1

				#####################
				sys.stdout.write(".")
				#####################

		### README_018 ###
		### CALCULATION OF BEST CASE / MIN RECOMBINATION VALUE FOR TRIO ###
		# for item_x in nr_marker_list:
		for item_x in nr_good_list:
			dk = 0
			min_flank_recomb[item_x] = 9.0
			while dk <= dk_array[item_x]:
				try:
					trio_string = best_recomb_trio[item_x,dk]
					trio_list = trio_string.split('\t')
					trio_value = int(trio_list[14])
					upper_value  = float(trio_list[1])
					# upper_value = round(upper_value,2)
					upper_value = round(upper_value,round_scale)
					### 'middle_value' IS A RECOMBINATION VALUE BETWEEN FLANKING MARKERS ###
					middle_value = float(trio_list[10])
					# middle_value = round(middle_value,2)
					middle_value = round(middle_value,round_scale)
					lower_value  = float(trio_list[5])
					# lower_value = round(lower_value,2)
					lower_value = round(lower_value,round_scale)
					### 'flank_value' IS A SUM OF ALL THREE VALUES IN TRIO ###
					flank_value = upper_value + middle_value + lower_value
					# flank_value = round(flank_value,2)
					flank_value = round(flank_value,round_scale)
					### FIND THE BEST (MINIMUM) ###
					if flank_value <= min_flank_recomb[item_x] and trio_value <= min_double_recomb[item_x]:
						min_flank_recomb[item_x] = flank_value
					dk = dk + 1
				except:
					mooba = booba
					dk = dk + 1

		# for item_x in nr_marker_list:

		### ARRAYS FOR PHYLOGRAPHER ###
		trio_best_counter = {}
		map_array_list = []
		map_array_matrix = {}
		###  BEST TRIOS EXTRACTION  ###
		### README_019 ###
		for item_x in nr_good_list:
			trio_best_counter[item_x] = 0
			dk = 0
			while dk <= dk_array[item_x]:
				try:
					flank_mark = "---"
					trio_string = best_recomb_trio[item_x,dk]
					trio_list = trio_string.split('\t')
					trio_value = int(trio_list[14])
					upper_marker = trio_list[0]
					upper_value  = float(trio_list[1])
					# upper_value = round(upper_value,2)
					upper_value = round(upper_value,round_scale)
					middle_value = float(trio_list[10])
					# middle_value = round(middle_value,2)
					middle_value = round(middle_value,round_scale)
					lower_marker = trio_list[8]
					lower_value  = float(trio_list[5])
					# lower_value = round(lower_value,2)
					lower_value = round(lower_value,round_scale)
					flank_value = upper_value + middle_value + lower_value
					# flank_value = round(flank_value,2)
					flank_value = round(flank_value,round_scale)

					flank_diff = flank_value - min_flank_recomb[item_x]
					if flank_diff <= 0.0001 and flank_diff >= -0.0001:
					# if (flank_value - min_flank_recomb[item_x]) == 0:
						flank_mark = "+++"
						trio_best_counter[item_x] = trio_best_counter[item_x] + 1
						# if trio_value <= double_limit:
						#	out_file3g.write(best_recomb_trio[item_x,dk] + '\t' + flank_mark + '\t' + str(flank_diff) + '\t' + '\n')
						# if trio_value > double_limit:
						#	out_file3h.write(best_recomb_trio[item_x,dk] + '\t' + flank_mark + '\t' + str(flank_diff) + '\t' + '\n')
					if trio_value <= min_double_recomb[item_x]:
						out_file3d.write(best_recomb_trio[item_x,dk] + '\t' + flank_mark + '\t' + \
										str(flank_diff) + '\t' + '\n')
						if trio_value <= double_limit and flank_mark == "+++":
							out_file3g.write(best_recomb_trio[item_x,dk] + '\t' + flank_mark + '\t' + \
											str(flank_diff) + '\t' + '\n')
							if item_x not in acceptable_markers:
								acceptable_markers.append(item_x)
								map_array_list.append(item_x)
							if upper_marker not in acceptable_markers:
								acceptable_markers.append(upper_marker)
								map_array_list.append(upper_marker)
							if lower_marker not in acceptable_markers:
								acceptable_markers.append(lower_marker)
								map_array_list.append(lower_marker)

							if item_x not in core_list:
                                                                core_list.append(item_x)

							### MATRIX DATA FOR PHYLOGRAPHER ###
							try:
								matrix_try = map_array_matrix[item_x,upper_marker]
							except:
								map_array_matrix[item_x,upper_marker] = upper_value

							try:
								matrix_try = map_array_matrix[upper_marker,item_x]
							except:
								map_array_matrix[upper_marker,item_x] = upper_value

							try:
								matrix_try = map_array_matrix[item_x,lower_marker]
							except:
								map_array_matrix[item_x,lower_marker] = lower_value

							try:
								matrix_try = map_array_matrix[lower_marker,item_x]
							except:
								map_array_matrix[lower_marker,item_x] = lower_value

							#####################################

						if trio_value > double_limit and flank_mark == "+++":
							out_file3h.write(best_recomb_trio[item_x,dk] + '\t' + flank_mark + '\t' + \
										str(flank_diff) + '\t' + '\n')
							if item_x not in very_bad_markers:
								very_bad_markers.append(item_x)
							# if upper_marker not in very_bad_markers:
							#	very_bad_markers.append(upper_marker)
							# if lower_marker not in very_bad_markers:
							#	very_bad_markers.append(lower_marker)

					dk = dk + 1
				except:
					mooba = booba
					dk = dk + 1

		### PHYLOGRAPHER MAP AND MATRIX FILES ###
		already_done_p1p2 = {}
		for item_p1 in map_array_list:
			trio_quality = "UNDEFINED"
			if trio_best_counter[item_p1] == 2:
				trio_quality = "UNIQ_TRIO_" + str(trio_best_counter[item_p1])
			if trio_best_counter[item_p1] >  2:
				trio_quality = "MULT_TRIO_" + str(trio_best_counter[item_p1])
			out_file3f.write(item_p1 + '\t' + trio_quality + '\n')
			for item_p2 in map_array_list:
				try:
					p1p2 = map_array_matrix[item_p1,item_p2]
					# already_done_p1p2[item_p1,item_p2] = "DONE"
					# already_done_p1p2[item_p2,item_p1] = "DONE"
				except:
					p1p2 = "NOTHING"

				try:
					query1 = already_done_p1p2[item_p1,item_p2]
				except:
					query1 = "NOT_YET"
				try:
					query2 = already_done_p1p2[item_p2,item_p1]
				except:
					query2 = "NOT_YET"

				if p1p2 != "NOTHING" and query1 == "NOT_YET" and query2 == "NOT_YET":
					dist_matrix_value = 1.0 - p1p2
					out_file3e1.write(item_p1 + '\t' + item_p2 + '\t' + str(round(dist_matrix_value,4)) + '\n')
					if trio_best_counter[item_p1] == 2 and trio_best_counter[item_p2] == 2:
						out_file3e2.write(item_p1 + '\t' + item_p2 + '\t' + str(round(dist_matrix_value,4)) + '\n')
					already_done_p1p2[item_p1,item_p2] = "DONE"
					already_done_p1p2[item_p2,item_p1] = "DONE"

	### SUMMARY OF MARKER SCORES / QUALITY ###
	for item_x in id_list:
		core_status   = "xxxxxxxxxxx"
		marker_status = "DEFAULT__STATUS"
		# marker_status = "DISTANT__MARKER"
		if map_construction == "TRUE":
			marker_status = "DISTANT__MARKER"
		try:
			master_id = dupl_markers_list[item_x]
			if master_id != item_x:
				dupl_status = " -- dupl -- "
			if master_id == item_x:
				dupl_status = "-=+master+=-"
		except:
			master_id = item_x
			dupl_status = "-=+unique+=-"
		try:
			drecomb = min_double_recomb[master_id]
		except:
			drecomb = 999
		# marker_status =       "DISTANT__MARKER"
		if master_id in very_bad_markers:
			marker_status = "BETTER_TO_AVOID"
		if master_id in acceptable_markers:
			marker_status = "GOOD_____MARKER"
		if master_id in core_list:
			core_status   = "CORE_MARKER"
		if master_id in really_bad_allele and master_id in really_bad_data:
			marker_status = "REALLY_ALL__BAD"
		if master_id in really_bad_allele and master_id not in really_bad_data:
			marker_status = "ALLELE_DIST_BAD"
		if master_id not in really_bad_allele and master_id in really_bad_data:
			marker_status = "DATA_LOSS___BAD"
		out_file3i.write(item_x + '\t' + master_id + '\t' + dupl_status + '\t' + marker_status + '\t' \
			+ `drecomb` + '\t' + core_status + '\t' + "***" + '\t' + `count_nr_A[master_id]` + '\t' + \
			`count_nr_B[master_id]` + '\t' + str(round(ratio_AB[master_id],round_scale)) + '\t' + bias_AB[master_id] + '\n')

	##########################################

	print ""
	print "============================================="	
	print "            ANALYSIS DONE                    "
	print "            ENJOY MAPPING                    "
	print "============================================="
	print ""

	in_file.close()
	out_file0.close()
	out_file1.close()
	out_file2.close()
	out_file4.close()
	out_file5.close()
	out_file6.close()
	out_file7.close()
	out_file8.close()
	out_file14.close()
	# if map_construction == "TRUE":
	out_file3a.close()
	out_file3b.close()
	if deep_clustering  == "TRUE":
		out_file3k.close()
		out_file3la.close()
		out_file3lb.close()
		out_file3m.close()
		out_file3p.close()
		out_file3r.close()
		out_file3s.close()
	if map_construction == "TRUE":
		out_file3c.close()
		out_file3d.close()
		out_file3e1.close()
		out_file3e2.close()
		out_file3f.close()
		out_file3g.close()
		out_file3h.close()
		out_file3i.close()
		# out_file3k.close()
		# out_file3l.close()
		# out_file3m.close()
		# out_file3p.close()
		# out_file3r.close()

#### POSITIVE CLUSTERING ####

def Seqs_Clustering(rec_cut, bit_cut, dat_cut, x_counter, matrix_array_current):

	global id_list
	global id_array
	global adj_array
	global already_done
	global group_count
	global dfs_counter
	global working_group
	global group_depth
	global node_count
	global node_array
	global sequence_array
	global sequence_array_bin
	global tree_clust_array
	global nr_good_list

	global pairs_array
	global frame_array
	global print_frame
	global round_scale

	global graph_depth
	global marker_depth

	global out_file14
	global out_file8

	global out_file3k
	global out_file3la
	global out_file3lb
	global out_file3m
	global out_file3p
	global out_file3r
	global out_file3s
	global init_len

	node_array = {}
	adj_array    = {}
	already_done = []
	group_count = 0
	pair_matrix_count = 0
	pair_matrix_array = {}

	print ""
	print "ALREADY DONE"
	print already_done
	print ""
	print ""
	print "MATRIX ARRAY"
	print len(matrix_array_current)
	print ""
	time.sleep(2)

	out_file10 = open(out_name + '.matrix' + "_" + x_counter, "wb")
	out_file12 = open(out_name + '.adj_list' + "_" + x_counter, "wb")
	out_file13 = open(out_name + '.group_info' + "_" + x_counter, "wb")

	if x_counter == "01":
		rec_cut = rec_cut
		out_file14.write("### CLUSTERING RESULTS ITERATION 01 REC CUTOFF: " + str(rec_cut) + " ###" + '\n')
	if x_counter == "02":
		rec_cut = rec_cut - 0.02
		out_file14.write("### CLUSTERING RESULTS ITERATION 02 REC CUTOFF: " + str(rec_cut) + " ###" + '\n')
	if x_counter == "03":
		rec_cut = rec_cut - 0.04
		out_file14.write("### CLUSTERING RESULTS ITERATION 03 REC CUTOFF: " + str(rec_cut) + " ###" + '\n')
	if x_counter == "04":
		rec_cut = rec_cut - 0.06
		out_file14.write("### CLUSTERING RESULTS ITERATION 04 REC CUTOFF: " + str(rec_cut) + " ###" + '\n')
	if x_counter == "05":
		rec_cut = rec_cut - 0.08
		out_file14.write("### CLUSTERING RESULTS ITERATION 05 REC CUTOFF: " + str(rec_cut) + " ###" + '\n')
	if x_counter == "06":
		rec_cut = rec_cut - 0.10
		out_file14.write("### CLUSTERING RESULTS ITERATION 06 REC CUTOFF: " + str(rec_cut) + " ###" + '\n')
	if x_counter == "07":
		rec_cut = rec_cut - 0.11
		out_file14.write("### CLUSTERING RESULTS ITERATION 07 REC CUTOFF: " + str(rec_cut) + " ###" + '\n')
	if x_counter == "08":
		rec_cut = rec_cut - 0.12
		out_file14.write("### CLUSTERING RESULTS ITERATION 08 REC CUTOFF: " + str(rec_cut) + " ###" + '\n')
	if x_counter == "09":
		rec_cut = rec_cut - 0.13
		out_file14.write("### CLUSTERING RESULTS ITERATION 09 REC CUTOFF: " + str(rec_cut) + " ###" + '\n')
	if x_counter == "10":
		rec_cut = rec_cut - 0.14
		out_file14.write("### CLUSTERING RESULTS ITERATION 10 REC CUTOFF: " + str(rec_cut) + " ###" + '\n')
	if x_counter == "11":
		rec_cut = rec_cut - 0.15
		out_file14.write("### CLUSTERING RESULTS ITERATION 11 REC CUTOFF: " + str(rec_cut) + " ###" + '\n')
	if x_counter == "12":
		rec_cut = rec_cut - 0.16
		out_file14.write("### CLUSTERING RESULTS ITERATION 12 REC CUTOFF: " + str(rec_cut) + " ###" + '\n')
	if x_counter == "13":
		rec_cut = rec_cut - 0.17
		out_file14.write("### CLUSTERING RESULTS ITERATION 13 REC CUTOFF: " + str(rec_cut) + " ###" + '\n')
	if x_counter == "14":
		rec_cut = rec_cut - 0.18
		out_file14.write("### CLUSTERING RESULTS ITERATION 14 REC CUTOFF: " + str(rec_cut) + " ###" + '\n')
	if x_counter == "15":
		rec_cut = rec_cut - 0.19
		out_file14.write("### CLUSTERING RESULTS ITERATION 15 REC CUTOFF: " + str(rec_cut) + " ###" + '\n')
	if x_counter == "16":
		rec_cut = rec_cut - 0.20
		out_file14.write("### CLUSTERING RESULTS ITERATION 16 REC CUTOFF: " + str(rec_cut) + " ###" + '\n')


	########################################
	###          CLUSTERING              ###
	########################################

	# print id_list
	print "------------------------------"
	print "FOUND " + `len(id_list)` + " UNIQ IDs"
	print "------------------------------"

	r = 0

	### CREATE NON-REDUNDANT MATRIX ###

	for key in pairs_array:
		id_a = pairs_array[key][0]
		id_b = pairs_array[key][1]
		try:
			cur_rec1 = float(matrix_array_current[id_a,id_b][0])
			cur_bit1 = int(matrix_array_current[id_a,id_b][1])
			cur_dat1_fr = float(matrix_array_current[id_a,id_b][2])
			cur_dat1 = float(matrix_array_current[id_a,id_b][3])
			query1 = 1
			print key
		except:
			# print "ALREADY PROCESSED"
			sys.stdout.write(".")
			query1 = 0

		try:
			cur_rec2 = float(matrix_array_current[id_b,id_a][0])
			cur_bit2 = int(matrix_array_current[id_b,id_a][1])
			cur_dat2_fr = float(matrix_array_current[id_b,id_a][2])
			cur_dat2 = float(matrix_array_current[id_b,id_a][3])
			r = r + 1
			# print `r` + " REVERSE PAIR FOUND"
			query2 = 1
		except:
			# print "NO REVERSE PAIR FOUND"
			query2 = 0

		if query1 == 1 and query2 == 0 and cur_rec1 <= rec_cut and cur_bit1 >= bit_cut and \
					cur_dat1 >= dat_cut and id_a != id_b:
			### STRING CONVERSION ###
			# cur_rec1_str = str(round(cur_rec1,2))
			# cur_dat1_str = str(round(cur_dat1,2))
			cur_rec1_str = str(round(cur_rec1,round_scale))
			cur_dat1_fr_str = str(round(cur_dat1_fr,round_scale))
			#########################
			out_file10.write(id_a + '\t' + id_b + '\t' + cur_rec1_str + '\t' + \
				`cur_bit1` + '\t' + cur_dat1_fr_str + '\n')
			pair_matrix_count = pair_matrix_count + 1
			current_matrix_pair = [id_a, id_b]
			pair_matrix_array[id_a,id_b] = current_matrix_pair
			### NEED TO UNSET ARRAY
			try:
				del matrix_array_current[id_a,id_b]
			except:
				# print "ALREADY REMOVED"
				sys.stdout.write(".")

		if query1 == 1 and query2 == 1 and id_a != id_b:
			if cur_rec1 <= cur_rec2:
				# print "CASE 1"
				if cur_dat1 >= dat_cut and cur_bit1 >= bit_cut and cur_rec1 <= rec_cut:
					### STRING CONVERSION ###
					# cur_rec1_str = str(round(cur_rec1,2))
					# cur_dat1_str = str(round(cur_dat1,2))
					cur_rec1_str = str(round(cur_rec1,round_scale))
					cur_dat1_fr_str = str(round(cur_dat1_fr,round_scale))
					#########################
					out_file10.write(id_a + '\t' + id_b + '\t' + cur_rec1_str + '\t' + \
						`cur_bit1` + '\t' + cur_dat1_fr_str + '\n')
					pair_matrix_count = pair_matrix_count + 1
					current_matrix_pair = [id_a, id_b]
					pair_matrix_array[id_a,id_b] = current_matrix_pair
					### NEED TO UNSET ARRAY
					try:
						del matrix_array_current[id_a,id_b]
						del matrix_array_current[id_b,id_a]
					except:
						# print "ALREADY REMOVED"
						sys.stdout.write(".")
			if cur_rec1 > cur_rec2:
				# print "CASE 2"
				if cur_dat2 >= dat_cut and cur_bit2 >= bit_cut and cur_rec2 <= rec_cut:
					### STRING CONVERSION ###
					# cur_rec2_str = str(round(cur_rec2,2))
					# cur_dat2_str = str(round(cur_dat2,2))
					cur_rec2_str = str(round(cur_rec2,round_scale))
					cur_dat2_fr_str = str(round(cur_dat2_fr,round_scale))
					#########################
					out_file10.write(id_b + '\t' + id_a + '\t' + cur_rec2_str + '\t' + \
					`cur_bit2` + '\t' + cur_dat2_fr_str + '\n')
					pair_matrix_count = pair_matrix_count + 1
					current_matrix_pair = [id_b, id_a]
					pair_matrix_array[id_b,id_a] = current_matrix_pair
					### NEED TO UNSET ARRAY
					try:
						del matrix_array_current[id_a,id_b]
						del matrix_array_current[id_b,id_a]
					except:
						# print "ALREADY REMOVED"
						sys.stdout.write(".")

	print "-------------------------------------"
	print `pair_matrix_count` + " PAIRS IN REDUNDANT MATRIX"
	print "-------------------------------------"
	print "BEGIN CLUSTERING"

	time.sleep(2)

	item_count = 0
	id_list.sort()
	### CREATE ADJACENCY LIST ###
	for item in id_list:
		item_count = item_count + 1
		print `item_count` + '\t' + item

		item_list = [item]
		for key in pair_matrix_array:
			id_a = pair_matrix_array[key][0]
			id_b = pair_matrix_array[key][1]
			if id_a == item:
				item_list.append(id_b)
			if id_b == item:
				item_list.append(id_a)
		adj_array[item] = item_list
		item_string = " ".join(item_list)
		out_file12.write(item_string + '\n')

	print "GROUP ANALYSIS"
	time.sleep(2)

	### GROUP ANALYSIS ###
	node_count = 0
	for item in id_list:
		# if item not in already_done:
		try:
			node_test = node_array[item]
		except:
			node_array[item] = 1
			group_count = group_count + 1
			already_done.append(item)
			node_count = node_count + 1
			working_group = [item]
			current_adj_list = adj_array[item]
			current_adj_len = len(current_adj_list)

			q = 0
			while q <= (current_adj_len - 1):
				current_adj_item = current_adj_list[q]
				# print `q` + '\t' + current_adj_item
				if current_adj_item in already_done:
					go_to_dfs = 0
				# if current_adj_item not in already_done:
				try:
					node_test = node_array[current_adj_item]
				except:
					node_array[current_adj_item] = 1
					already_done.append(current_adj_item)
					node_count = node_count + 1
					if current_adj_item not in working_group:
						working_group.append(current_adj_item)
					go_to_dfs = 1
					dfs_counter = 0
					# print 'Processing Group:  ' + `group_count`
					### README_015 ###
					DFS_procedure(current_adj_item)
				q = q + 1
			# if item not in already_done:
			#	already_done.append(item)
			working_group.sort()
			# print working_group
			print 'Processing Group:  ' + `group_count`
			print 'Number of processed nodes:  ' + `node_count`
			i = 0
			group_suffix = str(group_count)
			suffix_len = len(group_suffix)
			if suffix_len < 5:
				if suffix_len == 1:
					group_suffix = "0000" + group_suffix
				if suffix_len == 2:
					group_suffix = "000"  + group_suffix
				if suffix_len == 3:
					group_suffix = "00"   + group_suffix
				if suffix_len == 4:
					group_suffix = "0"    + group_suffix
				if suffix_len == 5:
					group_suffix = group_suffix
			## GRAPH TYPE: CONNECTED OR COMPLETE
			## NODE  TYPE: SATURATED OR DILUTED
			graph_type = "COMPLETE_GRAPH_" + group_suffix
			graph_depth[group_count] = graph_type
			for egg in working_group:
				current_adj_list7 = adj_array[egg]
				current_adj_len7 = len(current_adj_list7)
				current_group_len7 = len(working_group)
				# marker_depth[egg] = "SATURATED_NODE"
				marker_depth[egg] = "UNDEFINED_NODE"
				if current_adj_len7 == current_group_len7:
					marker_depth[egg] = "SATURATED_NODE"
				if current_adj_len7 != current_group_len7:
					graph_type = "LINKED___GROUP_" + group_suffix
					graph_depth[group_count] = graph_type
					marker_depth[egg] = "DILUTED___NODE"
				if current_group_len7 == 1:
					graph_type = "SINGLE____NODE_" + group_suffix
					graph_depth[group_count] = graph_type
					# marker_depth[egg] = "SINGLE____NODE"
			for element in working_group:
				current_adj_list1 = adj_array[element]
				current_adj_len1 = len(current_adj_list1)
				current_group_len1 = len(working_group)
				### TRY FRAME MARKER ID ###
				if print_frame == "TRUE":
					try:
						frame_marker = frame_array[element]
						frame_marker = ' ' + frame_marker + '_LinkageGroup '
					except:
						frame_marker = " __LinkageGroup "
				if print_frame == "FALSE":
					frame_marker = "_NONE_"
				###########################
				# tree_clust_array[element].append(str(group_count))
				tree_clust_array[element].append(group_count)
				if x_counter == "16":
					tree_clust_array[element].append(graph_depth[group_count])
					tree_clust_array[element].append(marker_depth[element])
				print tree_clust_array[element]
				###########################
				### WRITE DATA TO OUTPUT GROUP-INFO FILES ###
				if i == 0:
					out_file13.write(element + '\t' + `(current_adj_len1 - 1)` + '\t' \
						+ `current_group_len1` + '\t' + `group_count` + '\t' + '*****' \
						+ '\t' + frame_marker + '\t' + graph_depth[group_count] + '\t' \
						+ marker_depth[element] + '\n')
					out_file14.write(element + '\t' + `(current_adj_len1 - 1)` + '\t' \
						+ `current_group_len1` + '\t' + `group_count` + '\t' + '*****' \
						+ '\t' + frame_marker + '\t' + "*" + x_counter + "*" + '\t' \
						+ graph_depth[group_count] + '\t' + marker_depth[element] + '\n')
				if i != 0:
					out_file13.write(element + '\t' + `(current_adj_len1 - 1)` + '\t' \
						+ `current_group_len1` + '\t' + `group_count` + '\t' + '-----' \
						+ '\t' + frame_marker + '\t' + graph_depth[group_count] + '\t' \
						+ marker_depth[element] + '\n')
					out_file14.write(element + '\t' + `(current_adj_len1 - 1)` + '\t' \
						+ `current_group_len1` + '\t' + `group_count` + '\t' + '-----' \
						+ '\t' + frame_marker + '\t' + "*" + x_counter + "*" + '\t' \
						+ graph_depth[group_count] + '\t' + marker_depth[element] + '\n')
				i = i + 1
				###########################
	print '======================'
	# print already_done
	print len(already_done)
	print `group_count` + '\t' + 'GROUPS FOUND'
	print '======================'

	### BUILDING OF THE DUMMY CONSENSUS ###
	# print "BUILDING BINARY CONSENSUS DATA FILE"
	# time.sleep(2)

	if x_counter == "16" and rec_cut == 0:

		print "BUILDING BINARY CONSENSUS DATA FILE"
		time.sleep(2)

		consensus_array      = {}	; # TWO DIMENSIONAL ARRAY
		consensus_dict       = {}
		consensus_list       = []
		consensus_frame      = {}
		consensus_nr_string  = []	; # LIST OF STRINGS
		consensus_nr         = []	; # LIST OF NON_REDUNDANT IDs
		consensus_all        = {}	; # ARRAY OF SEQUENCES
		b = 1
		# for item in id_list:
		for item in nr_good_list:
			group_prefix = "GC_"
			group_id = "ABCDEFGH"
			len_b = len(str(b))
			### DUMMY GROUP COUNTER ###
			if len_b <= 5:
				if len_b == 1:
					group_id = group_prefix + "0000" + `b`
				if len_b == 2:
					group_id = group_prefix + "000"  + `b`
				if len_b == 3:
					group_id = group_prefix + "00"   + `b`
				if len_b == 4:
					group_id = group_prefix + "0"    + `b`
				if len_b == 5:
					group_id = group_prefix + ""     + `b`
			consensus_dict[group_id] = adj_array[item]
			consensus_list.append(group_id)

			### UNSET ALL VALUES IN ALL POSSIBLE CONSENSUSES ###
			p = 1
			while p < init_len:
				consensus_array[group_id,p] = "-"
				p = p + 1

			### GET REAL CONSENSUS SCORES ###
			current_adj_list = adj_array[item]
			current_adj_len = len(current_adj_list)
			c = 0
			while c <= (current_adj_len - 1):
				current_adj_item = current_adj_list[c]
				if marker_depth[current_adj_item] == "SATURATED_NODE":
					p = 1
					while p < init_len:
						inter_score = sequence_array_bin[current_adj_item,p]
						if inter_score != "-":
							consensus_array[group_id,p] = inter_score
						p = p + 1
				c = c + 1
			### INCR GROUP ID ###
			b = b + 1

		### SUMMARIZE THE CONSENSUS ###
		f = 1
		out_file3k.write(";" + '\t')
		out_file3m.write(";" + '\t')
		out_file3s.write(";" + '\t')
		while f < init_len:
			if f <  init_len-1:
				out_file3k.write(`f` + '\t')
				out_file3m.write(`f` + '\t')
				out_file3s.write(`f` + '\t')
			if f == init_len-1:
				out_file3k.write(`f` + '\n')
				out_file3m.write(`f` + '\n')
				out_file3s.write(`f` + '\n')
			f = f + 1
		
		for group_id in consensus_list:
			h = 1
			out_file3k.write(group_id + '\t')
			# out_file3m.write(group_id + '\t')
			out_file3la.write(group_id + '\t')
			out_file3lb.write(group_id + '\t')
			item_string = " ".join(consensus_dict[group_id])
			out_file3la.write(item_string + '\n')
			current_adj_list = consensus_dict[group_id]
			current_adj_len = len(current_adj_list)
			c = 0
			consensus_frame[group_id] = "-"
			while c <= (current_adj_len - 1):
				current_adj_item = current_adj_list[c]
				if marker_depth[current_adj_item] == "SATURATED_NODE":
					out_file3lb.write(current_adj_item + ' ')
					out_file3p.write(group_id + '\t' + current_adj_item + '\t' + "1.0" + '\n')
					t = 1
					out_file3m.write(current_adj_item + '\t')
					while t < init_len:
						current_score = sequence_array_bin[current_adj_item,t]
						if t <  init_len-1:
							out_file3m.write(current_score + '\t')
						if t == init_len-1:
							out_file3m.write(current_score + '\n')
						t = t + 1
					try:
						lg_fr = frame_array[current_adj_item]
						consensus_frame[group_id] = lg_fr
					except:
						lg_fr = "-"
				c = c + 1
			out_file3lb.write('\n')
			### WRITE CONSENSUS TO FILE ###
			out_file3m.write(group_id + '\t')
			dummy_consensus_l = []
			while h < init_len:
				dummy_consensus_l.append(consensus_array[group_id,h])
				if h <  init_len-1:
					out_file3k.write(consensus_array[group_id,h] + '\t')
					out_file3m.write(consensus_array[group_id,h] + '\t')
				if h == init_len-1:
					out_file3k.write(consensus_array[group_id,h] + '\n')
					out_file3m.write(consensus_array[group_id,h] + '\n')
				h = h + 1
			y = 0
			while y < init_len:
				if y <  init_len-1:
						out_file3m.write("=" + '\t')
				if y == init_len-1:
						out_file3m.write("=" + '\n')
				y = y + 1
			dummy_consensus_s = "\t".join(dummy_consensus_l)
			consensus_all[group_id] = dummy_consensus_s
			if dummy_consensus_s not in consensus_nr_string:
				consensus_nr_string.append(dummy_consensus_s)
				out_file3s.write(group_id + '\t' + dummy_consensus_s + '\n')
			out_file3r.write(consensus_frame[group_id] + '\t' + group_id + '\t' + "ND" + '\n')

	#######################################

	out_file10.close()
	out_file12.close()
	out_file13.close()

#### NEGATIVE CLUSTERING ####

def Seqs_Neg_Clustering(rec_cut, bit_cut, dat_cut, x_counter, matrix_array_current):

	global id_list
	global id_array
	global adj_array
	global already_done
	global group_count
	global dfs_counter
	global working_group
	global group_depth
	global node_count
	global node_array
	global sequence_array

	global pairs_array_N
	global frame_array
	global print_frame
	global round_scale

	node_array = {}
	adj_array    = {}
	already_done = []
	group_count = 0
	pair_matrix_count = 0
	pair_matrix_array = {}

	print ""
	print "ALREADY DONE"
	print already_done
	print ""
	print ""
	print "MATRIX ARRAY"
	print len(matrix_array_current)
	print ""
	time.sleep(2)

	if x_counter == "N1":
		rec_cut = 1.0 - rec_cut - 0.1
		bit_cut = -1*bit_cut
	if x_counter == "N2":
		rec_cut = 1.0 - rec_cut
		bit_cut = -1*bit_cut
	if x_counter == "N3":
		rec_cut = 1.0 - rec_cut + 0.1
		bit_cut = -1*bit_cut


	########################################
	###          CLUSTERING              ###
	########################################

	out_file10 = open(out_name + '.matrix' + "_" + x_counter, "wb")
	out_file12 = open(out_name + '.adj_list' + "_" + x_counter, "wb")
	out_file13 = open(out_name + '.group_info' + "_" + x_counter, "wb")

	# print id_list
	print "------------------------------"
	print "FOUND " + `len(id_list)` + " UNIQ IDs"
	print "------------------------------"

	r = 0

	### CREATE NON-REDUNDANT MATRIX ###

	for key in pairs_array_N:
		id_a = pairs_array_N[key][0]
		id_b = pairs_array_N[key][1]
		try:
			cur_rec1 = float(matrix_array_current[id_a,id_b][0])
			cur_bit1 = int(matrix_array_current[id_a,id_b][1])
			cur_dat1_fr = float(matrix_array_current[id_a,id_b][2])
			cur_dat1 = float(matrix_array_current[id_a,id_b][3])
			query1 = 1
			print key
		except:
			# print "ALREADY PROCESSED"
			sys.stdout.write(".")
			query1 = 0

		try:
			cur_rec2 = float(matrix_array_current[id_b,id_a][0])
			cur_bit2 = int(matrix_array_current[id_b,id_a][1])
			cur_dat2_fr = float(matrix_array_current[id_b,id_a][2])
			cur_dat2 = float(matrix_array_current[id_b,id_a][3])
			r = r + 1
			# print `r` + " REVERSE PAIR FOUND"
			query2 = 1
		except:
			# print "NO REVERSE PAIR FOUND"
			query2 = 0

		if query1 == 1 and query2 == 0 and cur_rec1 >= rec_cut and cur_bit1 <= bit_cut and \
					cur_dat1 >= dat_cut and id_a != id_b:
			### STRING CONVERSION ###
			# cur_rec1_str = str(round(cur_rec1,2))
			# cur_dat1_str = str(round(cur_dat1,2))
			cur_rec1_str = str(round(cur_rec1,round_scale))
			cur_dat1_fr_str = str(round(cur_dat1_fr,round_scale))
			#########################
			out_file10.write(id_a + '\t' + id_b + '\t' + cur_rec1_str + '\t' + \
				`cur_bit1` + '\t' + cur_dat1_fr_str + '\n')
			pair_matrix_count = pair_matrix_count + 1
			current_matrix_pair = [id_a, id_b]
			pair_matrix_array[id_a,id_b] = current_matrix_pair
			### NEED TO UNSET ARRAY
			try:
				del matrix_array_current[id_a,id_b]
			except:
				# print "ALREADY REMOVED"
				sys.stdout.write(".")

		if query1 == 1 and query2 == 1 and id_a != id_b:
			if cur_rec1 >= cur_rec2:
				# print "CASE 1"
				if cur_dat1 >= dat_cut and cur_bit1 <= bit_cut and cur_rec1 >= rec_cut:
					### STRING CONVERSION ###
					# cur_rec1_str = str(round(cur_rec1,2))
					# cur_dat1_str = str(round(cur_dat1,2))
					cur_rec1_str = str(round(cur_rec1,round_scale))
					cur_dat1_fr_str = str(round(cur_dat1_fr,round_scale))
					#########################
					out_file10.write(id_a + '\t' + id_b + '\t' + cur_rec1_str + '\t' + \
						`cur_bit1` + '\t' + cur_dat1_fr_str + '\n')
					pair_matrix_count = pair_matrix_count + 1
					current_matrix_pair = [id_a, id_b]
					pair_matrix_array[id_a,id_b] = current_matrix_pair
					### NEED TO UNSET ARRAY
					try:
						del matrix_array_current[id_a,id_b]
						del matrix_array_current[id_b,id_a]
					except:
						# print "ALREADY REMOVED"
						sys.stdout.write(".")
			if cur_rec1 < cur_rec2:
				# print "CASE 2"
				if cur_dat2 >= dat_cut and cur_bit2 <= bit_cut and cur_rec2 >= rec_cut:
					### STRING CONVERSION ###
					# cur_rec2_str = str(round(cur_rec2,2))
					# cur_dat2_str = str(round(cur_dat2,2))
					cur_rec2_str = str(round(cur_rec2,round_scale))
					cur_dat2_fr_str = str(round(cur_dat2_fr,round_scale))
					#########################
					out_file10.write(id_b + '\t' + id_a + '\t' + cur_rec2_str + '\t' + \
					`cur_bit2` + '\t' + cur_dat2_fr_str + '\n')
					pair_matrix_count = pair_matrix_count + 1
					current_matrix_pair = [id_b, id_a]
					pair_matrix_array[id_b,id_a] = current_matrix_pair
					### NEED TO UNSET ARRAY
					try:
						del matrix_array_current[id_a,id_b]
						del matrix_array_current[id_b,id_a]
					except:
						# print "ALREADY REMOVED"
						sys.stdout.write(".")

	print "-------------------------------------"
	print `pair_matrix_count` + " PAIRS IN REDUNDANT MATRIX"
	print "-------------------------------------"
	print "BEGIN CLUSTERING"

	time.sleep(2)

	item_count = 0
	id_list.sort()
	### CREATE ADJACENCY LIST ###
	for item in id_list:
		item_count = item_count + 1
		print `item_count` + '\t' + item

		item_list = [item]
		for key in pair_matrix_array:
			id_a = pair_matrix_array[key][0]
			id_b = pair_matrix_array[key][1]
			if id_a == item:
				item_list.append(id_b)
			if id_b == item:
				item_list.append(id_a)
		adj_array[item] = item_list
		item_string = " ".join(item_list)
		out_file12.write(item_string + '\n')


	print "GROUP ANALYSIS"

	time.sleep(2)

	### GROUP ANALYSIS ###
	node_count = 0
	for item in id_list:
		# if item not in already_done:
		try:
			node_test = node_array[item]
		except:
			node_array[item] = 1
			group_count = group_count + 1
			already_done.append(item)
			node_count = node_count + 1
			working_group = [item]
			current_adj_list = adj_array[item]
			current_adj_len = len(current_adj_list)

			q = 0
			while q <= (current_adj_len - 1):
				current_adj_item = current_adj_list[q]
				# print `q` + '\t' + current_adj_item
				if current_adj_item in already_done:
					go_to_dfs = 0
				# if current_adj_item not in already_done:
				try:
					node_test = node_array[current_adj_item]
				except:
					node_array[current_adj_item] = 1
					already_done.append(current_adj_item)
					node_count = node_count + 1
					if current_adj_item not in working_group:
						working_group.append(current_adj_item)
					go_to_dfs = 1
					dfs_counter = 0
					# print 'Processing Group:  ' + `group_count`
					DFS_procedure(current_adj_item)
				q = q + 1
			# if item not in already_done:
			#	already_done.append(item)
			working_group.sort()
			# print working_group
			print 'Processing Group:  ' + `group_count`
			print 'Number of processed nodes:  ' + `node_count`
			i = 0
			for element in working_group:
				current_adj_list1 = adj_array[element]
				current_adj_len1 = len(current_adj_list1)
				current_group_len1 = len(working_group)
				### TRY FRAME MARKER ID ###
				if print_frame == "TRUE":
					try:
						frame_marker = frame_array[element]
						frame_marker = ' ' + frame_marker + '_LinkageGroup '
					except:
						frame_marker = " __LinkageGroup "
				if print_frame == "FALSE":
					frame_marker = ""
				###########################
				if i == 0:
					out_file13.write(element + '\t' + `(current_adj_len1 - 1)` + '\t' \
						+ `current_group_len1` + '\t' + `group_count` + '\t' + '*****' \
						+ '\t' + frame_marker + '\n')
				if i != 0:
					out_file13.write(element + '\t' + `(current_adj_len1 - 1)` + '\t' \
						+ `current_group_len1` + '\t' + `group_count` + '\t' + '-----' \
						+ '\t' + frame_marker + '\n')
				i = i + 1
				###########################
	print '======================'
	# print already_done
	print len(already_done)
	print `group_count` + '\t' + 'GROUPS FOUND'
	print '======================'

	out_file10.close()
	out_file12.close()
	out_file13.close()

#############################
### README_015 ###
def DFS_procedure(current_adj_item):

	global adj_array
	global already_done
	global group_count
	global dfs_counter
	global working_group
	global group_depth
	global node_count
	global node_array

	print "DFS" + '\t' + `dfs_counter`

	for key in adj_array:
		current_adj_list = adj_array[key]
		current_adj_len = len(current_adj_list)
		if current_adj_item in current_adj_list:
			q = 0
			while q <= (current_adj_len - 1):
				current_adj_item1 = current_adj_list[q]
				# print `q` + '\t' + current_adj_item1
				if current_adj_item1 in already_done:
					go_to_dfs = 0
				# if current_adj_item1 not in already_done:
				try:
					node_test = node_array[current_adj_item1]
				except:
					node_array[current_adj_item1] = 1
					already_done.append(current_adj_item1)
					node_count = node_count + 1
					if current_adj_item1 not in working_group:
						working_group.append(current_adj_item1)
					go_to_dfs = 1
					DFS_procedure(current_adj_item1)
					### HOW MANY TIMES LOOP INSIDE ITSELF ###
					dfs_counter = dfs_counter + 1
					print "DFS" + '\t' + `dfs_counter`
					group_depth = dfs_counter
				q = q + 1

	######################################################################

import math
import re
import sys
import string
import time

global deep_clustering
global map_construction
global double_limit
global allele_dist
global abs_loss
global round_scale
global print_all_pairs

deep_clustering = "TRUE"
# deep_clustering = "FALSE"

map_construction = "TRUE"
# map_construction = "FALSE"

# print_all_pairs = "TRUE"
print_all_pairs = "FALSE"

double_limit = 3

### OLDREADME_006 ###
### THESE VARIABLES ARE DEFINED BY ARGUMENTS/OPTIONS IN THIS VERSION [SEE BELOW] ###
allele_dist  = 0.33	# 1:3
# allele_dist  = 0.25	# 1:4
abs_loss = 50

##################

# round_scale = 2
round_scale = 4
# round_scale = 6

if __name__ == "__main__":
	### README_006 ###
	if len(sys.argv) <= 10 or len(sys.argv) > 11:
		print "                                                                     "
		print "   PROGRAM  USAGE:                                                   "
		print "   MAD MAPPER TAKES 10 ARGUMENTS/OPTIONS IN THE FOLLOWING ORDER:     "
		print "  (1)input_file[LOC_DATA/MARKER SCORES]          (2)output_file[NAME]"
		print "  (3)rec_cut[0.2]          (4)bit_cut[100]       (5)data_cut[25]     "
		print "  (6)group_file[OPTIONAL]  (7)allele_dist[0.33]  (8)missing_data[50] "
		print "  (9)trio_analysis[TRIO/NOTRIO/NR_SET]           (10)double_cross[3] "
		print "                                                                     "
		print "   if group_file does not exist just enter [X]                       "
		print "                                                                     "
		print "   DEFAULT VALUES:    IN  OUT  0.2  100  25  X  0.33  50  NOTRIO  2  "
		print "                                                                     "
		input = raw_input("   TYPE \"HELP\" FOR HELP [ \"EXIT\" TO EXIT ] : ")
		if input == "HELP" or input == "help":
			print "                                                                             "
			print "              MAD MAPPER ARGUMENTS/OPTIONS - BRIEF EXPLANATION:              "
			print "                                                                             "
			print "       INPUT/OUTPUT FILES:                                                   "
			print "   [1] - Input File Name  (locus file with raw marker scores)                "
			print "   [2] - Output File Name (master name [prefix] for 80 or so output files)   "
			print "                                                                             "
			print "       CLUSTERING PARAMETERS (WILL AFFECT CLUSTERING/GROUPING ONLY):         "
			print "   [3] - Recombination Value (Haplotype Distance) cutoff:  0.20 - 0.25       "
			print "         (NOTE: TRIO analysis (see below) works with 0.2 rec_cut value only) "
			print "   [4] - BIT Score cutoff: 60-1000 [100 is default and highly recommended]   "
			print "         (Check README_MADMAPPER for BIT Scoring Matrix system and values)   "
			print "   [5] - Overlap Data cutoff (data_cut): minimum number of scores between    "
			print "         two markers to be compared to assign pairwise distance              "
			print "                                                                             "
			print "   [6] - Optional Frame Work Marker Map (very useful for clustering analysis "
			print "                              to assign new markers to known linkage groups) "
			print "                                                                             "
			print "       FILTERING PARAMETERS (WILL AFFECT MARKER FILTERING, CREATION OF GOOD  "
			print "                             NON-REDUNDANT SET OF MARKERS AND TRIO ANALYSIS):"
			print "   [7] - Allele Distortion: to filter markers with high allele distortion    "
			print "   [8] - Missing Data: how many missing scores are allowed per marker        "
			print "                                                                             "
			print "       TRIO (TRIPLET) ANALYSIS -                                             "
			print "       FINDING OF TIGHTLY LINKED MARKERS AND THEIR RELATIVE ORDER:           "
			print "   [9] - TRIO/NOTRIO (If TRIO option is chosen then TRIPLET analysis will    "
			print "                      take place. Is not recommended to use for large set    "
			print "                      markers, 1000 or greater)                              "
			print "         if NR_SET chosen there is no CLUSTERING at all,                     "
			print "            Non-Redundant locus file will be created only                    "
			print "   [10] - Number of Double Crossovers cutoff value for TRIPLET analysis:     "
			print "          3 - default for noisy data; 0 is recommended for perfect scores    "
			print "                                                                             "
			print "   CHECK README_MADMAPPER FOR DETAILED DESCRIPTION OF OPTIONS                "
			print "                           AND OUTPUT FILES FORMATS/STRUCTURE                "
			print "                                                                             "
			sys.exit()
		if input != "HELP" and input != "help":
			sys.exit()
		sys.exit()
	if len(sys.argv) == 11:
		in_name  = sys.argv[1]
		out_name = sys.argv[2]
		rec_cut  = sys.argv[3]
		bit_cut  = sys.argv[4]
		dat_cut  = sys.argv[5]
		frame_file_name = sys.argv[6]
		allele_d = sys.argv[7]
		missingd = sys.argv[8]
		trio_analysis = sys.argv[9]
		double_x = sys.argv[10]
		try:
			rec_cut = float(rec_cut)
		except:
			print ""
			print "3-d argument rec_cut must be a number between 0.20 - 0.25"
			print "    default:  0.2 is recommended                         "
			print ""
			sys.exit()
		try:
			bit_cut = int(bit_cut)
		except:
			print ""
			print "4-th argument bit_cut must be a number between 60 - 1000"
			print "    default:   100 is recommended                       "
			print ""
			sys.exit()
		# dat_cut = float(dat_cut)
		try:
			dat_cut = int(dat_cut)
		except:
			print ""
			print "5-th argument data_cut must be a number between 10 - 1000"
			print "    default:    25 is recommended                        "
			print "    default:    50 for large sets (50 is better than 25) "
			print ""
			sys.exit()
		try:
			allele_d = float(allele_d)
		except:
			print ""
			print "7-th argument allele_dist must be a number between 0.1 - 0.9"
			print "    default:  0.33 is recommended                           "
			print "              0.25 can be used with caution                 "
			print ""
			sys.exit()
		try:
			missingd = int(missingd)
		except:
			print ""
			print "8-th argument missing_data must be a number between 0 - 100"
			print "    default:  50 is recommended                            "
			print "    default:  10 for perfect scores                        "
			print ""
			sys.exit()
		try:
			double_x = int(double_x)
		except:
			print ""
			print "10-th argument double_cross must be a number between 0 - 25"
			print "    default:    3 for noisy data (mis-scored markers)      "
			print "    default:    0 for perfect scores                       "
			print ""
			sys.exit()

		if rec_cut > 0.25 or rec_cut < 0.2:
			print "======================================================"
			print "  Recombination cutoff must be between 0.25 and 0.20  "
			print "          default:  0.20 is recommended               "
			print "======================================================"
			sys.exit()

		# if dat_cut > 1.00 or dat_cut < 0.1:
		if dat_cut > 1000 or dat_cut < 10:
			print "======================================================"
			print "   DataPoints cutoff must be between 10 and 1000      "
			print "          default:   25 is recommended                "
			print "          default:   50 for perfect scores            "
			print "======================================================"
			sys.exit()

		if bit_cut > 1000 or bit_cut < 60:
			print "======================================================"
			print "     BIT Score cutoff must be between 1000 and 60     "
			print "          default:  100 is recommended                "
			print "======================================================"
			sys.exit()

		if allele_d < 0.1 or allele_d > 0.9:
			print "======================================================"
			print "     allele distortion must be between 0.1 and 0.9    "
			print "          default:  0.33 is recommended               "
			print "                    0.25 can be used with caution     "
			print "======================================================"
			sys.exit()

		if missingd < 0 or missingd > 100:
			print "======================================================"
			print "     missing data must be between 0 and 100           "
			print "          default:   50 is recommended                "
			print "          default:   10 for perfect scores            "
			print "======================================================"
			sys.exit()

		if double_x < 0 or double_x > 25:
			print "======================================================"
			print "  double rec trio cutoff must be between 0 and 25     "
			print "          default:   3 for noisy data                 "
			print "          default:   0 for perfect scores             "
			print "======================================================"
			sys.exit()

		if trio_analysis != "TRIO" and trio_analysis != "NOTRIO" and trio_analysis != "NR_SET":
			print "======================================================"
			print "      trio_analysis must be TRIO or NOTRIO            "
			print "     for large set of markers (greater than 1000)     "
			print "  TRIO analysis may work really long (several hours)  "
			print "======================================================"
			sys.exit()

		if trio_analysis == "TRIO":
			map_construction = "TRUE"
		if trio_analysis == "NOTRIO":
			map_construction = "FALSE"
		if trio_analysis == "NR_SET":
			map_construction = "NR_SET"

		allele_dist = allele_d
		abs_loss = missingd
		double_limit = double_x

		sys.setrecursionlimit(10000)

		Define_Bit_Scores()

		Read_Data_File(in_name, out_name, rec_cut, bit_cut, dat_cut, frame_file_name)

#### THE END ####
