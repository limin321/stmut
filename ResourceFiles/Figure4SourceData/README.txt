Patient]_[replicate]_SNPMajorAlleleCount.csv has 8 columns
Column1: "SNP" - the germline SNP within the spot in the form of 
"Chromosome:start-ReferenceAllele/AlternativeAllele"
Column2: "Spot" - the spot name of each spatial barcode.
Column3: "Major_Count" - the number of reads mapping to the major allele.
Column4: "Minor_Count" - the number of reads mapping to the minor allele.
Column5: "Total_Read_Count" - the toal reads mapping to that allele, which 
is the sum of # of Major_Count and Minor_Count.
Column6: "Fraction_Major" - the ratio of Major_Count / Total_Read_Count.
Column7: "Cluster" - the gene expression cluster the spot belongs to.
Column8: "Tumor" - either a tumor spot or a non-tumor spot.
Column9: "MajorAllele" - either TumorRef or TumorAlt to specify which allele is the major allele.
Column10: "Reference_Allele" - the reference allele.
Column11: "Tumor_Seq_Allele2" - the alternative allele.