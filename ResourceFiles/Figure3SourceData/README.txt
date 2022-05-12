BarcodeLegend.csv has 6 columns:
Column1: "barcode" - the 10x spatial barcode.                
Column2: "New spot name (after grouping)" - the new spot name after 
grouping in the format of spot-name_rep-cluster.
Column3: "Original spot name" -  the spot name that is given to each 
spatial barcode.           
Column4: "Total.Reads" - the number of total reads in a spot.         
Column5: "array.row" - the row position of a each barcode in the 10x 
visium capture area.             
Column6: "array.col" - the column position of a each barcode in the 10x 
visium capture area.    

[Patient]_[replicate]CNVs_Loupe.csv for Loupe visualization has 2 column:
Column1: "barcode" - 10x spatial barcode.
Column2: "Group" - one of the "Likely Tumor", "Possibly Tumor", "Unclear 
Identity", "Not Likely Tumor", "Out" five groups each spatial barcode 
belongs to based on CNVs signals in a spot.

BulkCellDNACopyNumberInference.cdt has 389450 rows and 3 colums.
Column 1: "CLID" - required for Java Tree visualization.
Column 2: "NAME" - the genomic regions and corresponding annotations.
Column 3: "Patient6 log2(tumor/reference)" - the CNVs signals (log2(tumor/reference)) in the bulk data of Patient 6.

SpotCopyNumberInference.cdt has 1506 rows and 1386 columns.
Column 1: GLID required for Java Tree visualization.
Column 2: genomic regions and annotations.
Column 3 and ...: the CNVs signals in each Spot at specific genomic 
regions.
