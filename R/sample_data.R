#' Sample Data for sptBClstRds Function
#'
#' @description First input for sptBClstRds
#'
#' @format A vector containing all spotIndex path
#' \describe{
#'   \item{files}{a character}
#' }
"files"


#' Sample Data for sptBClstRds Function
#'
#' @description Second input for sptBClstRds
#'
#' @format A data frame with 36601 rows and 11 variables
#' \describe{ ## document each column/variable; no "" is needed for the variable name, otherwise it will cause warning.
#'   \item{X}{the first colum is the gene name}
#'   \item{AAACAAGTATCTCCCA.1}{spatial barcode of spot1}
#'   \item{AGATACTCAAGATCGA.1}{spatial barcode of spot2}
#'   \item{AGATCTCAGGTGTGAT.1}{spatial barcode of spot3}
#'   \item{AGATGACTCGCCCACG.1}{spatial barcode of spot4}
#'   \item{AGCAACATATCTTATT.1}{spatial barcode of spot5}
#'   \item{CAGAACTTAGCCCTCT.1}{spatial barcode of spot6}
#'   \item{CCGTAGGGTTGTTTAC.1}{spatial barcode of spot7}
#'   \item{CCTCCTGTTGTGTCGT.1}{spatial barcode of spot8}
#'   \item{TGATGTCAATTAAGTG.1}{spatial barcode of spot9}
#'   \item{TTACACGATCTGCGAC.1}{spatial barcode of spot10}
#' }
"data1"

#' Sample Data for sptBClstRds Function
#'
#' @description Third input for sptBClstRds
#'
#' @format A data frame with 10 rows and 2 variables
#' \describe{
#'   \item{Barcode}{the spatial barcode of spots}
#'   \item{Graph.based}{the cluster of 10X Loupe}
#' }
"data2"

#' Sample Data for sptBClstRds Function
#'
#' @description Fourth input for sptBClstRds
#'
#' @format A path to tissue_positions_list.csv
#' \describe{
#'   \item{path}{a character}
#' }
"path"




#' Sample Data for sptMutCt Function
#'
#' @description First input for sptMutCt
#'
#' @format A data frame with 10 rows and 3 variables
#' \describe{
#'   \item{barcode}{the spatial barcode of spots}
#'   \item{spot}{the spot name}
#'   \item{TotalRDs}{the total number reads of each spot}
#' }
"index"


#' Sample Data for sptMutCt Function
#'
#' @description Second input for sptBClstRds
#'
#' @format A vector containing all spots samtools mpileup output path
#' \describe{
#'   \item{files2}{a character}
#' }
"files2"



#' Sample Data for nonZeRdCts Function
#'
#' @description First input for nonZeRdCts
#'
#' @format A data frame with 10 rows and 6 variables
#' \describe{
#'   \item{barcode}{the spatial barcode of spots}
#'   \item{cluster}{the cluster the spot belongs to}
#'   \item{spot}{the spot name}
#'   \item{TotalRDs}{the total number of reads of each spot}
#'   \item{array_row}{spot row position}
#'   \item{array_col}{spot col position}
#' }
"d1"

#' Sample Data for nonZeRdCts Function
#'
#' @description Second input for nonZeRdCts
#'
#' @format A data frame with 5 rows and 6 variables
#' \describe{
#'   \item{spot}{the spot name}
#'   \item{RreadC}{reference read counts}
#'   \item{MreadC}{Mutant reads counts}
#'   \item{Treads}{the sum of RreadC and MreadC}
#'   \item{groups}{if MreadC == 0, the spot will given "G0", if MreadC == 1, the spot will given "G1",
#'   .... if MreadC > 9, the spot will given "G10"}
#'   \item{barcode}{the spatial barcode of spots}
#' }
"d2"





#' Sample Data for spotSummary Function
#'
#' @description First input for spotSummary
#'
#' @format A data frame with 5 rows and 10 variables
#' \describe{
#'   \item{barcode}{the spatial barcode of spots}
#'   \item{spot}{the spot name}
#'   \item{cluster}{the cluster the spot belongs to}
#'   \item{TotalRDs}{the total number of reads of each spot}
#'   \item{array_row}{spot row position}
#'   \item{array_col}{spot col position}
#'   \item{RreadC}{reference read counts}
#'   \item{MreadC}{Mutant reads counts}
#'   \item{Treads}{the sum of RreadC and MreadC}
#'   \item{groups}{if MreadC == 0, the spot will given "G0", if MreadC == 1, the spot will given "G1",
#'   .... if MreadC > 9, the spot will given "G10"}
#' }
"df"

#' Sample Data for spotSummary Function
#'
#' @description Second input for spotSummary
#'
#' @format A path to spot Mpileup file
#' \describe{
#'   \item{path}{a character}
#' }
"path1"



#' Sample Data for wtArmMedianOne Function
#'
#' @description First input for wtArmMedianOne
#'
#' @format A data frame with 1516 rows and 9 variables
#' \describe{
#'   \item{chromosome}{the chromosome information}
#'   \item{start}{the start position}
#'   \item{end}{the end position}
#'   \item{gene}{the gene name}
#'   \item{log2}{log2 value}
#'   \item{depth}{depth}
#'   \item{gc}{gc content}
#'   \item{tx_length}{tx_length}
#'   \item{weight}{weight}
#' }
"cnr"

#' Sample Data for wtArmMedianOne Function
#'
#' @description Second input for wtArmMedianOne
#'
#' @format A data frame with 24 rows and 3 variables
#' \describe{
#'   \item{V1}{Cchromosome column}
#'   \item{V2}{the start position}
#'   \item{V3}{the end position}
#' }
"centm"



#' Sample Data for CtArmGenes and cdt_filt_sort Function
#'
#' @description First input for CtArmGenes and cdt_filt_sort
#'
#' @format A data frame with 1506 rows and 5 variables
#' \describe{
#'   \item{CLID}{CLID}
#'   \item{NAME}{Gene chr:start-end:name}
#'   \item{spot19_rep1_cluster1}{first spot}
#'   \item{spot424_rep1_cluster1}{second spot}
#'   \item{spot325_rep1_cluster1}{third spot}
#' }
"cdt"


#' Sample Data for CtArmGenes Function
#'
#' @description First input for CtArmGenes
#'
#' @format A data frame with 23 rows and 5 variables
#' \describe{
#'   \item{chromosome}{the chromosome information}
#'   \item{p_Genes}{number of genes in p arm}
#'   \item{q_Genes}{number of genes in q arm}
#'   \item{pArmEnds}{position pArm ends}
#'   \item{CM_row_pos}{which row the centromere is from}
#' }
"data3"



#' Sample Data for cdt_filt_sort Function
#'
#' @description First input for cdt_filt_sort
#'
#' @format A data frame with 23 rows and 5 variables
#' \describe{
#'   \item{arms}{arm name}
#'   \item{genes}{number of genes in each arm}
#'   \item{chromosome}{chromosome}
#'   \item{gene_row}{representative gene row}
#' }
"data4"

#' Sample Data for accStartCNR_CNS Function
#'
#' @description First input for accStartCNR_CNS
#'
#' @format A data frame with 389448 rows and 7 variables
#' \describe{
#'   \item{chromosome}{chromosome name}
#'   \item{start}{start position of DNA region}
#'   \item{end}{end position of DNA region}
#'   \item{log2}{log2 ratio}
#'   \item{weight}{gene weight}
#' }
"data5"

#' Sample Data for accStartCNR_CNS Function
#'
#' @description Second input for accStartCNR_CNS
#'
#' @format A data frame with 122 rows and 8 variables
#' \describe{
#'   \item{chromosome}{chromosome name}
#'   \item{start}{start position of DNA region}
#'   \item{end}{end position of DNA region}
#'   \item{log2}{log2 ratio}
#'   \item{weight}{gene weight}
#' }
"data6"

#' Sample Data for accumStartPos Function
#'
#' @description First input for accumStartPos
#'
#' @format A data frame with 13401 rows and 4 variables
#' \describe{
#'   \item{chromosome}{chromosome name}
#'   \item{start}{start position of DNA region}
#'   \item{end}{end position of DNA region}
#'   \item{tumorshift}{tumorshift deviation from 0.5 VAF}
#' }
"data7"

#' Sample Data for groupSpots and newSptBC Function
#'
#' @description First input for groupSpots
#'
#' @format A data frame with 389 rows and 6 variables
#' \describe{
#'   \item{barcode}{spot barcode}
#'   \item{cluster}{the cluster the spot belongs to}
#'   \item{spot}{spot name}
#'   \item{TotalRDs}{the number of reads in the spot}
#'   \item{array_row}{X-axis position}
#'   \item{array_col}{Y-axis position}
#' }
"data8"

#' Sample Data for proposTumLoup Function
#'
#' @description First input for proposTumLoup
#'
#' @format A data frame with 67 rows and 15 variables
#' \describe{
#'   \item{barcode}{spot barcode}
#'   \item{TotalRDs}{TotalRDs}
#'   \item{spot}{spot}
#'   \item{RreadC}{RreadC}
#'   \item{MreadC}{MreadC}
#'   \item{Treads}{Treads}
#'   \item{groups}{groups}
#'   \item{propM2Total}{propM2Total}
#'   \item{propM2SumRM}{propM2SumRM}
#'   \item{MoreThan1Mrd}{MoreThan1Mrd}
#'   \item{Score}{Score}
#'   \item{GenesWithMutRead}{GenesWithMutRead}
#'   \item{array_row}{array_row}
#'   \item{array_col}{array_col}
#'   \item{cluster}{cluster}
#' }
"data9"

#' Sample Data for proposTumLoup Function
#'
#' @description Second input for proposTumLoup
#'
#' @format A data frame with 744 rows and 3 variables
#' \describe{
#'   \item{barcode}{spot barcode}
#'   \item{spot}{spot}
#'   \item{TotalRDs}{TotalRDs}
#' }
"data10"
