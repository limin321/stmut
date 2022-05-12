The source data for figure 2A can be found in the figure 1 folder. The only difference between figure 2A and figure 1A are the zoomed insets that we chose to highlight.

The source data for figure 2B is included in two files: Fig2Bheatmap.cdt and MutationSpotReadCounts.csv. The Fig2Bheatmap.cdt file was used to create the tiling plot in figure 2B while the MutationSpotReadCounts.csv file contains additional information pertaining to the distribution of mutant versus non-mutant reads in each spot.

Fig2Bheatmap.cdt		
A 164 columns cdt file stores the presence of point mutations in a specific spot.
The first column "CLID" is for Javatree view visualization.
The second column is a list of point mutations. Note that we exclude mutations that were detected in the DNA-sequencing data but for which there was no coverage in the spatial transcriptomics data.
The rest columns are spots. For example, rep2_spot195 is the spot195 from rep2 of patient 4. See the source data in the figure 1 folder for a barcode legend describing the spots. Note that we did not include spots for which there were no mutant reads detected. 
The value 1 indicates there is at least 1 mutant reads detected in a spot at a specific mutation locus. 
The number 0 means no mutant reads are detected.

MutationSpotReadCounts.csv
This file is more comprehensive than the Fig2Bheatmap.cdt file. The heat map contains 1's and 0's to denote whether a given spot has any mutant reads detected or not. This file provides more detailed information with respect to the number of mutant and reference reads in each spot. This file also contains spots/mutations with reference reads but no mutant reads. It does not, however, list spots or mutations with 0 mutant and 0 reference reads.
A 97 row and 613 columns tables stores the number of a mutant basepairs of a spot.
Column1: the list of 97 point mutations. That is, each row is a point mutation.
the rest columns: the number of reads possessing a specific point mutations in a spot. For example, 
0mut5ref means 0 mutant reads and 5 reference reads are detected in a spot at a specific mutation 
locus.


