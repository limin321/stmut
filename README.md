
<!-- README.md is generated from README.Rmd. Please edit that file -->

# stmut: Somatic Mutation Investigation of Spatial Transcriptomics Data

------------------------------------------------------------------------

<!-- badges: start -->
<!-- badges: end -->

Characterizing gene expression profiles throughout tissue space provides
key insights in investigating biological processes and disease
development, including cancer. Bioinformatic tools exploring and
interpreting spatial transcriptomics data are in great need -
especially, approaches to visualize point mutations, allelic imbalance,
and copy number variations (CNVs). CNVkit is a popular toolkit used to
investigate the copy number alterations in both DNA-seq and RNA-seq
data. Based on
[CNVkit-RNA](https://cnvkit.readthedocs.io/en/stable/rna.html) and
SAMtools, we provide an R package called stmut via this github page. The
stmut package includes a series of functions to visualize copy number
variations (CNVs), point mutations, and allelic imbalance in spatial
transcriptomics data. We also provide the [scripts producing the
figures](https://github.com/limin321/stmut/blob/master/FigTableScripts/FigTables.md)
in the manuscript, which also serves as a user guide for this package.
In addition, this package is also applicable to 10x single cell data
analyses. <br />

The functions in the stmut package are organized into 3 parts: CNVs,
point mutations, and allelic imbalance.

This package was tested using R version 4.1.1, a macOS Monterey, Apple
M1, 16G Memory. Given that spatial transcriptomics data normally have
more than hundreds or thousands spots, we recommend using a high
performance cluster to obtain point mutation and allelic imbalance for
each spot.

## Installation

------------------------------------------------------------------------

You can install the development version of stmut from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("limin321/stmut")
library(stmut)
```

#### Notes

-   Bash scripts displayed in `echo` command are for your reference when
    you run your own data. <br />
-   This package relies on previously sequenced DNAseq data, for
    example, exome data. That is, you need to have your bulk CNVs,
    germline SNPs, and somatic mutations list ready before using this
    package. <br />
-   Prepare the following 5 files from the spaceranger pipeline output:
    <br />

1.  filtered_feature_bc.csv <br />
2.  Graph-Based.csv, this file is exported from 10X Loupe Browser as
    shown below. <br />
    ![](./FigTableScripts/Images/Export_Graph_Based.png)
3.  possorted_genome_bam.bam <br />
4.  spatial/tissue_positions_list.csv <br />
5.  raw_feature_bc_matrix/barcodes.tsv.gz <br />

## I. Point Mutation Detection

------------------------------------------------------------------------

-   spotIndex generation: you can also run splitSpot() to generate an
    individual spot barcode and gene expression file, and each file is
    named numerically. For example, the first spot is spot000.txt, the
    next is spot001.txt and so forth. <br />

``` r
file <- read.csv("./Rep1/Data/SpacerangerOutput/CloupeFilesManualAlignment/filtered_feature_bc.csv")
splitSpot(file = file)
```

output of splitSpot(). spotIndex contains individual spot barcode txt
file; txt directory contains individual spot gene expression profile.
![](./FigTableScripts/Images/splitSpot.png)

-   spotBam generation: the spot bam is generated as suggested by
    [10xGenomics subset-bam](https://github.com/10XGenomics/subset-bam)
    <br />

``` bash
echo "subset-bam_linux --bam possorted_genome_bam.bam --cell-barcodes spot000.txt --out-bam spot000.bam"
echo "samtools index spot000.bam"
#> subset-bam_linux --bam possorted_genome_bam.bam --cell-barcodes spot000.txt --out-bam spot000.bam
#> samtools index spot000.bam
```

-   Count point mutations for each spot: we count the number of ref and
    mut reads using Mpileup_RNA.pl script [found
    here](https://github.com/limin321/stmut/tree/master/FigTableScripts).
    This scripts takes 3 inputs as shown in the following example. The
    first is the somatic mutation list; the second is the spot bam file;
    the third is the reference fasta file, which should be the same used
    either SpaceRanger or CellRanger. Make sure `samtools` is installed
    before running: <br />

``` bash
echo "perl Mpileup_RNA.pl Patient4SomaticSNPs.txt spot000/spot000.bam ./refdata-gex-GRCh38-2020-A/fasta/genome.fa"
#> perl Mpileup_RNA.pl Patient4SomaticSNPs.txt spot000/spot000.bam ./refdata-gex-GRCh38-2020-A/fasta/genome.fa
```

-   [spaPointMutation](https://github.com/limin321/stmut/blob/master/FigTableScripts/FigTables.md#figure-1-figure-2a)
    creates a folder in your working directory including 8 files related
    to spot point mutations exploration. The AllSptTumPropsed.csv file
    contains a list of point mutations for visualization on the 10X
    Loupe Browser. The color scheme can be customized in the 10X Loupe
    Browser. The figures generated should be similar to Figure 1 in our
    manuscript. Make sure the format of your input files matches the
    examples provided by the package to ensure the smooth running of the
    codes. <br />

#### Examples of running individual functions

The following are basic examples on how to run each individual function:
1. sptBClstRds()

``` r
library(stmut)
## basic example code
files <- list.files(path = "./inst/extdata/spotIndex", pattern = ".txt", full.names = TRUE, recursive = FALSE)
data1 <- read.csv("./inst/extdata/filtered_feature_bc.csv", header = TRUE)
data2 <- read.csv("./inst/extdata/Graph-Based.csv", header = TRUE)
path <- "./inst/extdata/tissue_positions_list.csv"

df <- sptBClstRds(files=files,data1=data1, data2=data2, path=path)

head(df[[1]]) # sptBcRds.csv
#>              barcode    spot TotalRDs
#> 1 AAACAAGTATCTCCCA-1 spot000    22655
#> 2 AGATACTCAAGATCGA-1 spot100    10631
#> 3 AGATCTCAGGTGTGAT-1 spot101      797
#> 4 AGATGACTCGCCCACG-1 spot102      260
#> 5 AGCAACATATCTTATT-1 spot103      180
#> 6 CAGAACTTAGCCCTCT-1 spot212     5019
head(df[[2]]) # sptBcSummary.csv
#>              barcode   cluster    spot TotalRDs array_row array_col
#> 1 AAACAAGTATCTCCCA-1 Cluster 3 spot000    22655        50       102
#> 2 AGATACTCAAGATCGA-1 Cluster 4 spot100    10631        46        82
#> 3 AGATCTCAGGTGTGAT-1 Cluster 2 spot101      797        29        75
#> 4 AGATGACTCGCCCACG-1 Cluster 2 spot102      260        21        55
#> 5 AGCAACATATCTTATT-1 Cluster 2 spot103      180        41        93
#> 6 CAGAACTTAGCCCTCT-1 Cluster 4 spot212     5019        53        85
```

2.  sptMutCt()

``` r
index <- read.csv("./inst/extdata/spotBC.csv", header = TRUE)
files2 <- list.files(path = "./inst/extdata/mpileup",pattern = "MpileupOutput_RNA.txt", full.names = TRUE, recursive = TRUE, include.dirs = TRUE)
df1 <- sptMutCt(index = index, files = files2)

head(df1[[1]]) # spotRdCtFinal.csv
#>       spot RreadC MreadC Treads groups            barcode
#> 8  spot266      0      5      5     G5 CCTCCTGTTGTGTCGT-1
#> 7  spot259      0      4      4     G4 CCGTAGGGTTGTTTAC-1
#> 9  spot662      0      3      3     G3 TGATGTCAATTAAGTG-1
#> 6  spot212      0      2      2     G2 CAGAACTTAGCCCTCT-1
#> 10 spot699      0      2      2     G2 TTACACGATCTGCGAC-1
#> 1  spot000      0      0      0     G0 AAACAAGTATCTCCCA-1
head(df1[[2]]) # spotMutReadCount.csv
#>                 index groups
#> 8  CCTCCTGTTGTGTCGT-1     G5
#> 7  CCGTAGGGTTGTTTAC-1     G4
#> 9  TGATGTCAATTAAGTG-1     G3
#> 6  CAGAACTTAGCCCTCT-1     G2
#> 10 TTACACGATCTGCGAC-1     G2
#> 1  AAACAAGTATCTCCCA-1     G0
```

3.  nonZeRdCts()

``` r
d1 <- read.csv("./inst/extdata/spotBClster.csv", header = TRUE)
d2 <- read.csv("./inst/extdata/spotRdCtFinal.csv", header = TRUE)
d2 <- na.omit(d2)

df3 <- nonZeRdCts(df1 = d1, df2 = d2)
head(df3[[1]]) # NonZeroRdSpotIndex.csv
#>              barcode    spot   cluster TotalRDs array_row array_col RreadC
#> 1 CAGAACTTAGCCCTCT-1 spot212 Cluster 4     5019        53        85      0
#> 2 CCGTAGGGTTGTTTAC-1 spot259 Cluster 3    10758        60       104      0
#> 3 CCTCCTGTTGTGTCGT-1 spot266 Cluster 3    19768        62        96      0
#> 4 TGATGTCAATTAAGTG-1 spot662 Cluster 3     8780        59       103      0
#> 5 TTACACGATCTGCGAC-1 spot699 Cluster 4     3699        62       110      0
#>   MreadC Treads groups
#> 1      2      2     G2
#> 2      4      4     G4
#> 3      5      5     G5
#> 4      3      3     G3
#> 5      2      2     G2
head(df3[[2]]) # MutantSpotIndex.csv
#>              barcode    spot   cluster TotalRDs array_row array_col RreadC
#> 1 CAGAACTTAGCCCTCT-1 spot212 Cluster 4     5019        53        85      0
#> 2 CCGTAGGGTTGTTTAC-1 spot259 Cluster 3    10758        60       104      0
#> 3 CCTCCTGTTGTGTCGT-1 spot266 Cluster 3    19768        62        96      0
#> 4 TGATGTCAATTAAGTG-1 spot662 Cluster 3     8780        59       103      0
#> 5 TTACACGATCTGCGAC-1 spot699 Cluster 4     3699        62       110      0
#>   MreadC Treads groups
#> 1      2      2     G2
#> 2      4      4     G4
#> 3      5      5     G5
#> 4      3      3     G3
#> 5      2      2     G2
```

4.  spotSummary()

``` r
df <- read.csv("./inst/extdata/NonZeroRdSpotIndex.csv", header = TRUE)
path1 <-"./inst/extdata/mpileup/"

df4 <- spotSummary(df = df, path1 = path1)

head(df4) # NonZeroRdSpotIndex_Scored.csv
#>              barcode    spot   cluster TotalRDs array_row array_col RreadC
#> 1 TTACACGATCTGCGAC-1 spot699 Cluster 4     3699        62       110      0
#> 2 CAGAACTTAGCCCTCT-1 spot212 Cluster 4     5019        53        85      0
#> 3 CCGTAGGGTTGTTTAC-1 spot259 Cluster 3    10758        60       104      0
#> 4 TGATGTCAATTAAGTG-1 spot662 Cluster 3     8780        59       103      0
#> 5 CCTCCTGTTGTGTCGT-1 spot266 Cluster 3    19768        62        96      0
#>   MreadC Treads groups propM2Total propM2SumRM MoreThan1Mrd    Score
#> 1      2      2     G2    5.406867           1            1 6.406867
#> 2      2      2     G2    3.984858           1            1 4.984858
#> 3      4      4     G4    3.718163           1            1 4.718163
#> 4      3      3     G3    3.416856           1            1 4.416856
#> 5      5      5     G5    2.529340           1            1 3.529340
#>                                                  GenesWithMutRead
#> 1                        WNT4-1:22120068 C>T;WNT4-1:22120070 C>T;
#> 2                        WNT4-1:22120068 C>T;WNT4-1:22120070 C>T;
#> 3                        WNT4-1:22120068 C>T;WNT4-1:22120070 C>T;
#> 4                                           RFNG-17:82049969 G>A;
#> 5 WNT4-1:22120068 C>T;WNT4-1:22120070 C>T;CTNNBIP1-1:9850735 C>T;
```

Column meaning of df4: <br /> column 1: Spatial barcode <br /> column 2:
Spot name <br /> column 3: Cluster in Loupe.Loupe file <br /> column 4:
Total reads each spot has <br /> column 5: Row position in
position_lists_bc.csv <br /> column 6: Column position in
position_lists_bc.csv <br /> column 7: Reference read counts <br />
column 8: Mutant read counts <br /> column 9: Sum of Ref read counts and
Mut read counts (sum of col 7 and col 8) <br /> column 10: The number (
\< 10) in this column indicates how many mutant reads in a spot. For
example, G9 means there are 9 mutant reads in that spot. However, G10
means there are10 or more mutant reads in the spot. <br /> column 11:
Proportion of # of mutant reads to Total reads (col8/col4) <br /> column
12: Proportion of # of mutant reads to sum of ref and mut reads
(col8/col9) <br /> column 13: A binary column, when the spot has \> 1
mutant read, it is 1, otherwise, it is 0. <br /> column 14: A score
calculated based on the # of ref and mutant read counts. The higher the
score is, the more likely it is a tumor spot. The table is sorted by
this column. <br /> column 15: Details on which genes contain what
mutations in which spots. <br />

## II. Copy Number Variation Detection

------------------------------------------------------------------------

accStartCNR_CNS() and bulkCNVs() functions are used to generate
arm-level bulk CNVs plot. <br /> wtArmMedianOne(), CtArmGenes() and
cdt_filt_sort() functions are used in generating spot CNVs heatmap.
Kindly refer to the case study in our manuscript FigTableScripts
[here](https://github.com/limin321/stmut/blob/master/FigTableScripts/FigTables.md#spatial-spot-cnvs-figure-3)
<br />

Here are examples on how to run the following functions: <br /> 1.
wtArmMedianOne()

``` r
df <- read.csv("./inst/extdata/spot1_rep1.cnr", header = TRUE)
path1 <- read.csv("./inst/extdata/hg38_centromereSimple.bed", sep = "\t", header = FALSE) 

df5 <- wtArmMedianOne(data = cnr, centmere = centm)

head(df5[[1]])
#>   chromosome   start     end     gene       log2     depth     gc tx_length
#> 1          1 1001138 1014541    ISG15 -0.2183919 2.2842600 0.5643       788
#> 2          1 1216908 1232031     SDF4 -0.2183919 0.0000000 0.6135       604
#> 3          1 1373730 1375495 AURKAIP1 -0.2183919 0.1142860 0.6687       875
#> 4          1 1385711 1399328    CCNL2 -0.2183919 0.0538503 0.5137      1857
#> 5          1 1401908 1407313   MRPL20 -0.2183919 0.0000000 0.5309       314
#> 6          1 1541673 1574869    SSU72 -0.2183919 0.0239464 0.5387      4176
#>     weight
#> 1 0.600865
#> 2 0.814714
#> 3 0.792194
#> 4 0.728547
#> 5 0.798626
#> 6 0.799400
head(df5[[2]])
#>   chromosome p_Genes q_Genes  pArmEnds CM_row_pos
#> 1          1     100      83 145917714        100
#> 2          2      43      52  96184859         43
#> 3          3      32      46 100492619         32
#> 4          4      11      22  55395957         11
#> 5          5       8      47  73498408          8
#> 6          6      47      24  73515750         47
```

2.  CtArmGenes()

``` r
cdt <- read.table("./inst/extdata/cdt.cdt", header = TRUE)
d3 <- read.table("./inst/extdata/summary.txt", sep = "\t", header = TRUE) 

df6 <- CtArmGenes(cdt = cdt, data = d3)

head(df6)
#>   arms genes chromosome gene_row
#> 1   1p   100          1        1
#> 2   1q    83          1      101
#> 3   2p    43          2      184
#> 4   2q    52          2      227
#> 5   3p    32          3      279
#> 6   3q    46          3      311
```

3.  cdt_filt_sort()

``` r
library(dplyr)
#> Warning: package 'dplyr' was built under R version 4.1.2
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
cdt <- read.table("./inst/extdata/cdt.cdt", header = TRUE)
data4 <- read.csv("./inst/extdata/CtArmGenSummary.csv",  header = TRUE)
arm <- c("1p","3p","3q","4q","5q","8q","9q","10p","10q","11q","13p","13q", "20p","20q","21q","14q","17q")
d4 <- data4 %>% filter(arms %in% arm)
genes <- d4[,"genes"]
rs <- d4[,"gene_row"]
gainLoss <- c(1,-1,1,-1,-1,1,1,-1,-1,1,-1,-1,1,1,-1,1,1)

df7 <- cdt_filt_sort(cdt = cdt,genes = genes,gainLoss = gainLoss,rs=rs)
#> Warning in xtfrm.data.frame(x): cannot xtfrm data frames

head(df7)
#>      CLID                       NAME spot424_rep1_cluster1
#> 1 IMAGE:0    1:1001138-1014541:ISG15             0.7610323
#> 2 IMAGE:1     1:1216908-1232031:SDF4             0.7610323
#> 3 IMAGE:2 1:1373730-1375495:AURKAIP1             0.7610323
#> 4 IMAGE:3    1:1385711-1399328:CCNL2             0.7610323
#> 5 IMAGE:4   1:1401908-1407313:MRPL20             0.7610323
#> 6 IMAGE:5    1:1541673-1574869:SSU72             0.7610323
#>   spot325_rep1_cluster1 spot19_rep1_cluster1
#> 1           -0.01796844          -0.05211491
#> 2           -0.01796844          -0.05211491
#> 3           -0.01796844          -0.05211491
#> 4           -0.01796844          -0.05211491
#> 5           -0.01796844          -0.05211491
#> 6           -0.01796844          -0.05211491
```

4.  cnvs()

``` r
# data1: annotation data. A 2 column dataframe with column names as "spotName" and "annotate". Also, the annotate column must have "normal" cells or spots which will be used to center the data.
# cdt file is the weighted median cdt file generated by cnvkit.py.
# data2 is the output of wtArmMedianOne() function stored in the wtcnr/summary.txt

cnvs(data1 = data1, cdt = cdt0, data2 = data2, arm = arm, gainLoss = gainLoss) # the output is CNVs.cdt in your working directory.
```

## III. Allelic Imbalance

------------------------------------------------------------------------

accumStartPos() and bulkLOHplot() functions are for generating bulk
DNAseq allelic imbalance plots.

-   Generate ‘samtools mpileup’ input of counting major- and minor-
    reads per mutant of each spot.

``` r
# Tumor SNPs list
data1 <- read.table(file = "/Volumes/Bastian/Limin/Ji_data/Patient6/BulkDNASeq/LOH/MpileupOutput_TumorConverted.txt", sep = "\t",quote = "", header = TRUE)

# generate "samtools mpileup" input for counting major and minor alleles per mutant of each Spot
lohMpileupInput(data1 = data1) # the LOHmpileupInput.txt file will generate in your working dir
```

In our cases, the patient4_hg38_SNPs.txt and patient6_hg38_SNPs.txt
files, which can be found
[here](https://github.com/limin321/stmut/blob/master/VisualizingSomaticAlterations/DNAseqResourceFiles/),
are used to count the # of major and minor alleles of each spot in
patient4 and patient6.

-   Counting the # of majorAllele- and minorAllele- reads per mutant of
    each spot. The script Mpileup_RNA_alleImbalance.pl can be downloaded
    [here](https://github.com/limin321/stmut/tree/master/FigTableScripts)

``` bash
echo "perl ./Mpileup_RNA_alleImbalance.pl ./LOHmpileupInput.txt spot000/spot000.bam"

#> perl ./Mpileup_RNA_alleImbalance.pl ./LOHmpileupInput.txt spot000/spot000.bam
```

-   Generate a summary table of all spot major/minor allele counts of
    all spots.

``` r
files <- c("/Volumes/Bastian/Limin/Ji_data/Patient6/SpatialTranscriptomic/Rep1/LOH/allelicImbalance2/mpileupOutput/spot0001/MpileupOutput_RNA.txt","/Volumes/Bastian/Limin/Ji_data/Patient6/SpatialTranscriptomic/Rep1/LOH/allelicImbalance2/mpileupOutput/spot0002/MpileupOutput_RNA.txt")

x <- files[1]
y = match("spot0001",str_split_fixed(x,"/",15)) # 12

lohMajorAlleleCt(files = files, y=12)
```

The output is 2 csv files: SNPallMajorAlleleCount.csv and
SNPMajorAlleleCount.csv. The latter is used to generate Figures in the
manuscript. <br />

-   Scripts generating the allelic imbalance figures(Figure 4 and Figure
    S6) in the manuscript can be found
    [here](https://github.com/limin321/stmut/blob/master/FigTableScripts/FigTables.md#bulk-dnaseq-allelic-imbalance)
