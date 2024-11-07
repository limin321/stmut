
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

- Bash scripts displayed in `echo` command are for your reference when
  you run your own data. <br />
- This package relies on previously sequenced DNAseq data, for example,
  exome data. That is, you need to have your bulk CNVs, germline SNPs,
  and somatic mutations list ready before using this package. <br />
- Prepare the following 5 files from the spaceranger pipeline output:
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

- spotIndex generation: you can also run splitSpot() to generate an
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

- spotBam generation: the spot bam is generated as suggested by
  [10xGenomics subset-bam](https://github.com/10XGenomics/subset-bam)
  <br />

``` bash
echo "subset-bam_linux --bam possorted_genome_bam.bam --cell-barcodes spot000.txt --out-bam spot000.bam"
echo "samtools index spot000.bam"
#> subset-bam_linux --bam possorted_genome_bam.bam --cell-barcodes spot000.txt --out-bam spot000.bam
#> samtools index spot000.bam
```

- Count point mutations for each spot: we count the number of ref and
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

- [spaPointMutation](https://github.com/limin321/stmut/blob/master/FigTableScripts/FigTables.md#figure-1-figure-2a)
  creates a folder in your working directory including 8 files related
  to spot point mutations exploration. The AllSptTumPropsed.csv file
  contains a list of point mutations for visualization on the 10X Loupe
  Browser. The color scheme can be customized in the 10X Loupe Browser.
  The figures generated should be similar to Figure 1 in our manuscript.
  Make sure the format of your input files matches the examples provided
  by the package to ensure the smooth running of the codes. <br />

## II. Copy Number Variation Detection

------------------------------------------------------------------------

To call copy number variation from 10X spatial or single cell data. We
published 2 docker images, one built under ubuntu 20.04 (amd64), the
other in MacOS A pple M1 chip (arm64).

    Usage: /usr/local/bin/stmutcnv.sh cnv \
        --filteredFeatureCSV <value> \
        --clusterCSV <value> \
        --positionCSV <value> \
        --TotalReads <value> \
        --numSpots <value> \
        --group <value> \
        --annotate <value> \
        --arms <value> \
        --gainLoss <value> -\
        --pmtimes <value> \
        --clean <value>

### Details of each argument

| Argument \[default\] | Description                                                  |
|----------------------|--------------------------------------------------------------|
| `filteredFeatureCSV` | “filtered_feature_bc.csv” file                               |
| `clusterCSV`         | “Graph-Based.csv” file                                       |
| `positionCSV`        | “tissue_positions_list.csv” file                             |
| `numSpots` \[8\]     | \[opt\]Number of spots used for grouping                     |
| `TotalReads`\[1000\] | \[opt\]Number of reads or genes of a new spot after grouping |
| `group` \[gene\]     | \[opt\] One of ‘gene’,‘read’, ‘none’                         |
| `annotate`           | A two-column normal, tumor annotated csv file                |
| `arms`               | \[opt\] A list of arms; ex: 3p. 3q, 6p, 6q                   |
| `gainLoss`           | \[opt\] A list of 1, -1; ex: 1,-1,1,1                        |
| `pmtimes` \[100\]    | \[opt\] Permutation time when bulk DNA provided              |
| `clean` \[true\]     | \[opt\] clean intermediate files                             |

opt = optional

##### Pick your case type and refer to the example code to prepare yours. <br />

| case | group by | bulk tumor DNA-CNV data | Example |
|------|----------|-------------------------|---------|
| 1    | `gene`   | YES                     | \(1\)   |
| 2    | `gene`   | NO                      | \(2\)   |
| 3    | `read`   | YES                     | \(1\)   |
| 4    | `read`   | NO                      | \(2\)   |
| 5    | `none`   | YES                     | \(3\)   |
| 6    | `none`   | NO                      | \(4\)   |

## Ready to run your data.

### Step 1, download Docker or Singularity image.

Please pull the docker image from Docker Hub
[here](https://hub.docker.com/search?q=limin321) <br /> Ubutun:<br />
`docker pull limin321/stmutcnv_amd64:0.0.1` <br />
`singularity pull stmut.sig docker://limin321/stmutcnv_amd64:0.0.1`
<br />

Mac: <br /> `docker pull limin321/stmutcnv_arm64:0.0.1`<br />
`singularity pull stmut.sig docker://limin321/stmutcnv_arm64:0.0.1`
<br />

### Step 2, prepare the following 4 input files, better in the same working dir:

    1) filtered_feature_bc.csv, <br />     2) Graph-Based.csv, <br />
    3) tissue_positions_list.csv, <br />     4) annotate.csv. It should
look like below <br />

| cluster  | annotate |
|----------|----------|
| Cluster1 | tumor    |
| Cluster2 | normal   |
| Cluster3 | tumor    |
| Cluster4 | normal   |
| …        | …        |

The first 1) and third 3) inputs are standard Spaceranger outputs. The
second 2) is exported from the Loupe browser as shown above. The fourth
4) input is a two column csv file annotated by you for each cluster.
Please note the annoation “normal” and “tumor” should be little case.

### step 3: Infer CNV from 10X Spatial Data – 10X Platform

**Example 1**, group by ‘gene’ or ‘read’ ***with bulk tumor DNA data CNV
available***. Assuming the four inputs file are in the ‘inputs’ folder
inside ‘your_local_dir’

    docker run --rm -v <your_local_dir>:/home/stmut stmutcnv_arm64:0.0.1 bash /usr/local/bin/stmutcnv.sh cnv \
        --filteredFeatureCSV ./inputs/filtered_feature_bc.csv \
        --clusterCSV ./inputs/Graph-Based.csv \
        --positionCSV ./inputs/tissue_positions_list.csv \
        --TotalReads 1000 \
        --numSpots 8 \
        --group gene \ # change to 'read' if you want to group by read counts.
        --annotate ./inputs/annotate.csv \
        --arms 3p,6q,9q \
        --gainLoss -1,-1,-1 \
        --pmtimes 20

Expect output:

    analysis
    └── grouped_spots
        ├── BarcodeLegend.csv
        ├── cdt
        │   ├── CNVs_OrganizedByGEcluster_UMIcount.cdt
        │   ├── CNVs_OrganizedByGEcluster_UMIcount.pdf
        │   ├── CNVs_RankedBySimilarityToDNA.cdt
        │   ├── CNVs_RankedBySimilarityToDNA_CNVscoreHistogram.csv
        │   ├── CNVs_RankedBySimilarityToDNA_CNVscoreHistogram.pdf
        │   ├── CNVs_RankedBySimilarityToDNA_QQplot.pdf
        │   ├── CNVs_RankedbySimilaritytoDNA_Quintiles4Loupe.csv
        │   ├── CNVs_clustered.Rdata
        │   ├── CNVs_clustered_heatmap.pdf
        │   └── permutCNV_summ.csv
        └── histogram_genes_per_spot.png

    3 directories, 12 files

We provide three sets of outputs: <br />
**CNVs_OrganizedByGEcluster_UMIcount**: the three files are used to
generate Fig.4 included in our paper. <br />
**CNVs_RankedBySimilarityToDNA**: this set is optional. Only when you
provide bulk tumor CNVs data, these outputs will be generated. <br />
**CNVs_clustered**: we provide a dendrogram of CNVs info. The details
are saved in the .Rdata which you can extract by running the following
**R** codes. <br />

    load("./analysis_readgp/grouped_spots/cdt/CNVs_clustered.Rdata") # this will load the htp obj in R
    df1 <- read.table("./analysis/grouped_spots/cdt/CNVs_RankedBySimilarityToDNA.cdt", header = TRUE)
    data <- htp$carpet
    data1 <- cbind(df1[,1:2], data)
    write.table(data1, file = "./analysis/grouped_spots/cdt/CNVs_cluster.cdt", sep = "\t", row.names = FALSE)

**Example 2**, group by ‘gene’ or ‘read’ ***without bulk tumor DNA
data***. Assuming the four inputs file are in the ‘inputs’ folder inside
‘your_local_dir’

    docker run --rm -v <your_local_dir>:/home/stmut stmutcnv_arm64:0.0.1 bash /usr/local/bin/stmutcnv.sh cnv \
        --filteredFeatureCSV ./inputs/filtered_feature_bc.csv \
        --clusterCSV ./inputs/Graph-Based.csv \
        --positionCSV ./inputs/tissue_positions_list.csv \
        --TotalReads 1000 \
        --numSpots 8 \
        --group read \ # change to 'gene' if you want to group by gene counts.
        --annotate ./inputs/annotate.csv \

Expected outputs:

    analysis_grp_read_NObulk
    └── grouped_spots
        ├── BarcodeLegend.csv
        ├── cdt
        │   ├── CNVs_OrganizedByGEcluster_UMIcount.cdt
        │   ├── CNVs_OrganizedByGEcluster_UMIcount.pdf
        │   ├── CNVs_clustered.Rdata
        │   └── CNVs_clustered_heatmap.pdf
        └── histogram_genes_per_spot.png

    3 directories, 6 files

**Example 3**: no grouping spots is performed. Bulk tumor CNV data is
provided. Useful for single-cell data.

    docker run --rm -v <your_local_dir>:/home/stmut stmutcnv_arm64:0.0.1 bash /usr/local/bin/stmutcnv.sh cnv \
        --filteredFeatureCSV ./inputs/filtered_feature_bc.csv \
        --clusterCSV ./inputs/Graph-Based.csv \
        --positionCSV ./inputs/tissue_positions_list.csv \
        --TotalReads 1000 \
        --numSpots 8 \
        --group none \
        --annotate ./inputs/annotate.csv \
        --arms 3p,6q,9q \
        --gainLoss -1,-1,-1 \

*Expect outputs are the same as Example 1.*

**Example 4**: no grouping spots is performed. No Bulk tumor CNV data is
provided.

    docker run --rm -v <your_local_dir>:/home/stmut stmutcnv_arm64:0.0.1 bash /usr/local/bin/stmutcnv.sh cnv \
        --filteredFeatureCSV ./inputs/filtered_feature_bc.csv \
        --clusterCSV ./inputs/Graph-Based.csv \
        --positionCSV ./inputs/tissue_positions_list.csv \
        --TotalReads 1000 \
        --numSpots 8 \
        --group none \
        --annotate ./inputs/annotate.csv \

*Expected outputs are the same as Example 2.*

***One example of running in singularity:***

    singularity exec --bind <your_dir_to_mount>:/home/stmut --pwd /home/stmut <path/to>/stmut.sig bash stmutcnv.sh cnv \
        --filteredFeatureCSV ./inputs/filtered_feature_bc.csv \
        --clusterCSV ./inputs/Graph-Based.csv \
        --positionCSV ./inputs/tissue_positions_list.csv \
        --TotalReads 1000 --numSpots 8 \
        --group gene \
        --annotate ./inputs/annotate.csv

### StereoSeq Platform

If your data is from StereoSeq, you need to do some extra work. Here is
the step by step instructions.

###### First, run the following code to generate input files similar to 10X platform.

    docker run --rm  -v ./stmutCNVtest/scripts/:/home/stmut/ stmutcnv:latest bash /usr/local/bin/stmutcnv.sh gemconvert \
        --gemfile ./stereo/<chipID>.tissue.gem.gz \ # the output from stereo SAW pipeline.
        --binsize 200 \ # the bin_size, bin200 is assumed to similar size as 10X visium spot size.
        --outpath ./stereo/ # path to save the outputs.

It takes 3 arguments: the tissue.gem.gz; the bin_size, the output dir
you want to store the output.

Output 4 files. <br />     1) filtered_feature_bc.csv <br />     2)
graph_based.csv <br />     3) tissue_positions_list.csv <br />     4)
bin200_seurat.RDS <br /> The counterfeit barcodes in each file were
created to mimic the ones from 10X platform to keep consistent data
format for analysis. The real corresponding coordinates are stored in
the meta data of the rds. With this rds, you also need to perform
clustering analysis, and annotate each cluster so as to create a
annotate.csv file to run CNV analysis in the next step.

###### Second, infer CNV by following the same codes shown in the 10X Platform.

## III. Allelic Imbalance

------------------------------------------------------------------------

accumStartPos() and bulkLOHplot() functions are for generating bulk
DNAseq allelic imbalance plots.

- Generate ‘samtools mpileup’ input of counting major- and minor- reads
  per mutant of each spot.

``` r
# Tumor SNPs list
data1 <- read.table(file = "/Volumes/Bastian/Limin/Ji_data/Patient6/BulkDNASeq/LOH/MpileupOutput_TumorConverted.txt", sep = "\t",quote = "", header = TRUE)

# generate "samtools mpileup" input for counting major and minor alleles per mutant of each Spot
lohMpileupInput(data1 = data1) # the LOHmpileupInput.txt file will generate in your working dir
```

In our cases, the patient4_hg38_SNPs.txt and patient6_hg38_SNPs.txt
files, which can be found
[here](https://github.com/limin321/stmut/blob/master/VisualizingSomaticAlterations/DNAseqResourceFiles/),
are used to count the \# of major and minor alleles of each spot in
patient4 and patient6.

- Counting the \# of majorAllele- and minorAllele- reads per mutant of
  each spot. The script Mpileup_RNA_alleImbalance.pl can be downloaded
  [here](https://github.com/limin321/stmut/tree/master/FigTableScripts)

``` bash

echo "perl ./Mpileup_RNA_alleImbalance.pl ./LOHmpileupInput.txt spot000/spot000.bam"

#> perl ./Mpileup_RNA_alleImbalance.pl ./LOHmpileupInput.txt spot000/spot000.bam
```

- Generate a summary table of all spot major/minor allele counts of all
  spots.

``` r
files <- c("/Volumes/Bastian/Limin/Ji_data/Patient6/SpatialTranscriptomic/Rep1/LOH/allelicImbalance2/mpileupOutput/spot0001/MpileupOutput_RNA.txt","/Volumes/Bastian/Limin/Ji_data/Patient6/SpatialTranscriptomic/Rep1/LOH/allelicImbalance2/mpileupOutput/spot0002/MpileupOutput_RNA.txt")

x <- files[1]
y = match("spot0001",str_split_fixed(x,"/",15)) # 12

lohMajorAlleleCt(files = files, y=12)
```

The output is 2 csv files: SNPallMajorAlleleCount.csv and
SNPMajorAlleleCount.csv. The latter is used to generate Figures in the
manuscript. <br />

- Scripts generating the allelic imbalance figures(Figure 4 and Figure
  S6) in the manuscript can be found
  [here](https://github.com/limin321/stmut/blob/master/FigTableScripts/FigTables.md#bulk-dnaseq-allelic-imbalance)

``` r
sessionInfo()
#> R version 4.1.1 (2021-08-10)
#> Platform: x86_64-apple-darwin17.0 (64-bit)
#> Running under: macOS Big Sur 10.16
#> 
#> Matrix products: default
#> BLAS:   /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRblas.0.dylib
#> LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib
#> 
#> locale:
#> [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> loaded via a namespace (and not attached):
#>  [1] digest_0.6.29   lifecycle_1.0.3 magrittr_2.0.3  evaluate_0.16  
#>  [5] rlang_1.1.1     stringi_1.7.8   cli_3.4.1       rstudioapi_0.14
#>  [9] vctrs_0.6.2     rmarkdown_2.16  tools_4.1.1     stringr_1.5.0  
#> [13] glue_1.6.2      xfun_0.39       yaml_2.3.5      fastmap_1.1.0  
#> [17] compiler_4.1.1  htmltools_0.5.3 knitr_1.40
```
