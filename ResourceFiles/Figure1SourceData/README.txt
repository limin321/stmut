The source data has two types of files.
1. [Patient]_[replicate]_BarcodeLegend.csv. The legend connects various identifiers used to describe individual spots.
2. [Patient]_[replicate]_SpotMutCount.csv. The SpotMutCount file is formatted so that it can be opened directly in the Loupe browser to visualize the localization of spots with various combinations of mutant versus non-mutant reads, as shown in figure 1.

[Patient]_[replicate]_BarcodeLegend.csv has 6 columns:
Column1: "barcode" - the 10x spatial barcode.                
Column2: "gene.expression.cluster" - the gene expression cluster, assigned during the space ranger workflow, that each spot represented by the spatial barcode belongs to.
Column3: "spot.name" -  a name for each spot/barcode that we developed to simplify references to a given spot.           
Column4: "Total.Reads" - the number of total reads in a spot.         
Column5: "array.row" - the row position of a each barcode in the 10x visium capture area.             
Column6: "array.col" - the column position of a each barcode in the 10x visium capture area.            

[Patient]_[replicate]_SpotMutCount.csv has 2 column:
Column1: "barcode" - the 10x spatial barcode.
Column2: "TumorGp" - the tumor group each barcode is assigned to based on the number of mutant 
reads being detected in the spot. For example, "1Mut" indicates the spot has 1 mutant read; "2+Mut" indicates the spot contains 2 or more mutant reads. Spots with 0 mutant reads are further divided into 6 groups, listed as "0Mut0Ref" (0 mutant reads and 0 Reference reads detected), ..., "0Mut4Ref" (0 mutant reads and 4 reference reads detected); however, "0Mut5+Ref" indicates the spot has 0 mutant reads and 5 or more reference reads.



