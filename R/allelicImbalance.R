#' Set genomic start position to accumulated start position.
#'
#' @param dataf the df must have 2 cols, "chromosome" character show as 1,2,...,"X","Y",
#' the "start" column must be numeric.
#'
#' @return A data with accumulated start position, called in bulkLOHplot func.
#' @export
#'
#' @examples \dontrun{
#' loh <- read.table(system.file("extdata/", "p6LOH.txt", package = "stmut"),
#'   sep = "\t",
#'   header = TRUE
#' )
#' loh <- loh[, c(1, 2, 3, 10)]
#' colnames(loh) <- c("chromosome", "start", "end", "tumorshift")
#' loh$chromosome <- substr(loh$chromosome, 4, nchar(loh$chromosome))
#' df8 <- accumStartPos(loh)
#' }
.accumStartPos <- function(dataf) {
  if (!is.data.frame(dataf)) {
    flog.warn("the dataf must have 2 cols, 'chromosome' character show as 1,2,...,'X','Y', the 'start' column must be numeric.")
    stop("dataf must be a data frame.")
  }

  chromosome <- NULL
  chr <- unique(dataf$chromosome)[-1]
  new_ch <- dataf %>% filter(chromosome %in% "1")

  for (ch in chr) { # for loop to generate accumulated start position
    if (ch == "X") {
      val <- max(new_ch$start)
      valEnd <- max(new_ch$end)
      df <- dataf %>% filter(chromosome %in% ch)
      df$start <- as.numeric(df$start) + as.numeric(val)
      df$end <- as.numeric(df$end) + as.numeric(valEnd)
      df$chromosome <- "23"
      new_ch <- rbind(new_ch, df)
    } else if (ch == "Y") {
      val <- max(new_ch$start)
      valEnd <- max(new_ch$end)
      df <- dataf %>% filter(chromosome %in% ch)
      df$start <- as.numeric(df$start) + as.numeric(val)
      df$end <- as.numeric(df$end) + as.numeric(valEnd)
      df$chromosome <- "24"
      new_ch <- rbind(new_ch, df)
    } else if (ch == "MT") {
      val <- max(new_ch$start)
      valEnd <- max(new_ch$end)
      df <- dataf %>% filter(chromosome %in% ch)
      df$start <- as.numeric(df$start) + as.numeric(val)
      df$end <- as.numeric(df$end) + as.numeric(valEnd)
      df$chromosome <- "25"
      new_ch <- rbind(new_ch, df)
    } else {
      val <- max(new_ch$start)
      valEnd <- max(new_ch$end)
      df <- dataf %>% filter(chromosome %in% ch)
      df$start <- as.numeric(df$start) + as.numeric(val)
      df$end <- as.numeric(df$end) + as.numeric(valEnd)
      new_ch <- rbind(new_ch, df)
    }
  }
  new_ch$chromosome <- as.numeric(new_ch$chromosome) # this code is very important to plot the chromosomes in order.
  return(new_ch)
}

#' Allelic Imbalance or loss of heterozygocity plot.
#'
#' @param centmere hg38 centromere data
#' @param alle_imbal tumorshif/shift, which is the variant allele frequency deviated from 50% VAF
#'
#' @importFrom stringr str_split_fixed
#'
#' @return An arm-level LOH plot
#' @export
#'
#' @examples \dontrun{
#' centm <- read.csv(system.file("extdata/", "hg38_centromereSimple.bed", package = "stmut"),
#'   sep = "\t", header = FALSE
#' )
#' loh <- read.table(system.file("extdata/", "p6LOH.txt", package = "stmut"),
#'   sep = "\t",
#'   header = TRUE
#' )
#' loh <- loh[, c(1, 2, 3, 10)]
#' colnames(loh) <- c("chromosome", "start", "end", "tumorshift")
#' loh$chromosome <- substr(loh$chromosome, 4, nchar(loh$chromosome))
#' bulkLOHplot(centmere = centm, alle_imbal = loh)
#' }
bulkLOHplot <- function(centmere, alle_imbal) {
  if (!is.data.frame(centmere)) {
    stop("centmere must be a data frame!")
  }
  if (!is.data.frame(alle_imbal)) {
    stop("alle_imbal must be a data frame.")
  }

  chromosome <- NULL
  avg <- NULL
  start <- NULL
  median <- NULL
  colnames(centmere) <- c("chromosome", "start", "end")
  centmere$chromosome <- substr(centmere$chromosome, 4, nchar(centmere$chromosome))
  centmere[which(centmere$chromosome == "X"), 1] <- "23"
  centmere[which(centmere$chromosome == "Y"), 1] <- "24"
  centmere$chromosome <- as.numeric(centmere$chromosome)
  centmere <- centmere[order(centmere$chromosome), ] # adjust centmere data to be in chromosome order
  centmere$chromosome <- as.character(centmere$chromosome)

  # retrieve centromere row position of each chromosome of my sample
  pos_centM <- c()
  arm_median <- c()

  colIdx <- grep("chromosome", colnames(alle_imbal))

  alle_imbal[which(alle_imbal$chromosome == "X"), colIdx] <- "23" #
  alle_imbal[which(alle_imbal$chromosome == "Y"), colIdx] <- "24" # either deleting or keeping this code depends on male or female

  # create chromosome arm tumorshift median matrix
  armMedian <- data.frame(
    chromosome = character(),
    arm = character(),
    median = double(),
    start = double(),
    end = double(), stringsAsFactors = FALSE
  )

  chrs <- unique(alle_imbal$chromosome)
  for (ch in chrs) {
    d1 <- alle_imbal %>% filter(chromosome %in% ch)
    d1$start <- as.numeric(d1$start)
    d2 <- centmere %>% filter(chromosome %in% ch)
    for (pos in d1$start) {
      if (pos >= d2$start) {
        val <- which(d1$start == pos) # retrieve row position where the centromere should locate
        tumshft <- d1$tumorshift

        if (val == 1) {
          med <- round(median(tumshft[1:length(tumshft)]), 4)
          newR <- c(ch, "arm", med, min(d1$start), max(d1$end))
          armMedian[nrow(armMedian) + 1, ] <- newR
        } else {
          med1 <- round(median(tumshft[1:(val - 1)]), 4)
          newR1 <- c(ch, "arm1", med1, min(d1[c(1:(val - 1)), ]$start), max(d1[c(1:(val - 1)), ]$end))
          armMedian[nrow(armMedian) + 1, ] <- newR1
          med2 <- round(median(tumshft[val:length(tumshft)]), 4)
          newR2 <- c(ch, "arm2", med2, min(d1[c(val:length(tumshft)), ]$start), max(d1[c(val:length(tumshft)), ]$end))
          armMedian[nrow(armMedian) + 1, ] <- newR2
        }
        break
      }
    }
  }

  # adjust start position by calling accumStartPos
  alle_imbal[which(alle_imbal$chromosome == "23"), colIdx] <- "X" #
  alle_imbal[which(alle_imbal$chromosome == "24"), colIdx] <- "Y"

  data <- .accumStartPos(alle_imbal)
  df <- data

  # adjust arm start, end position fo armMedian data for plotting
  armMedian$chromosome <- as.numeric(armMedian$chromosome)
  armMedian$median <- as.numeric(armMedian$median)
  armMedian$start <- as.numeric(armMedian$start)
  armMedian$end <- as.numeric(armMedian$end)

  new_armD <- armMedian %>% filter(chromosome %in% 1)
  new_CM <- centmere %>% filter(chromosome %in% 1)

  nums <- unique(armMedian$chromosome)[-1]
  for (i in nums) {
    darM <- armMedian %>% filter(chromosome %in% i)
    data <- df %>% filter(chromosome %in% (i - 1))
    val1 <- max(data$start)
    val2 <- max(data$end)
    darM$start <- darM$start + val1
    darM$end <- darM$end + val2
    new_armD <- rbind(new_armD, darM)

    # adjust centromere start, end position for plotting
    CMdata <- centmere %>% filter(chromosome %in% i)
    CMdata$start <- CMdata$start + val1
    CMdata$end <- CMdata$end + val2
    new_CM <- rbind(new_CM, CMdata)
  }

  # preparing for plot
  if (length(unique(df$chromosome)) == 24) {
    chrom <- factor(c(paste0("chr", seq(22)), "chrX", "chrY"), levels = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"))
  } else {
    chrom <- factor(c(paste0("chr", seq(22)), "chrX"), levels = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX"))
  }

  pos <- df %>%
    group_by(chromosome) %>%
    summarize(avg = round(mean(start))) %>%
    pull(avg)
  t1 <- df %>%
    group_by(chromosome) %>%
    summarize(max = max(start)) %>%
    pull(max)

  ggplot(df) +
    annotate("point", x = df$start, y = df$tumorshift, size = 0.00002, color = "grey", pch = 0.2, alpha = 0.5) + # each row in data
    annotate("segment", x = 0, xend = max(df$start) * 1.02, y = 0, yend = 0, color = "black") + # x-axis
    annotate("segment", x = t1, xend = t1, y = 0, yend = 0.25, color = "black") +
    annotate("segment", x = new_armD$start, xend = new_armD$end, y = new_armD$median, yend = new_armD$median, color = "orange") + # add tumorshift median for each arm
    geom_vline(data = new_CM, mapping = aes(xintercept = start), linetype = 10, color = "black") + # centromere line
    scale_x_continuous(breaks = pos, labels = chrom, expand = c(0, 0)) +
    labs(x = "chromosome", y = "TumorShift") +
    ylim(c(0, 0.25)) + # if change ylim(), making sure also change the annotate() yend=?, otherwise, the chromosome boundaries will disappear
    theme_classic() + # use white color as background
    theme(axis.text.x = element_text(angle = 90))
}


#' Create samtools_Mpileup_input for counting spot allelic imbalance.
#'
#' @param data1 Tumor DNAseq germline SNPs ref/Alt counts
#'
#' @importFrom utils write.table
#'
#' @return A dataframe serves as the input for samtools mpileup
#' @export
#'
#' @examples \dontrun{
#' print("check user-guide")
#' }
lohMpileupInput <- function(data1) {
  if (!is.data.frame(data1)) {
    stop("data1 must be a data frame.")
  }

  data2 <- data.frame(matrix(nrow = 0, ncol = 13))
  colnames(data2) <- c(
    "Hugo_Symbol", "Chromosome", "Start_Position", "End_Position",
    "Variant_Classification", "Variant_Type", "Reference_Allele",
    "Tumor_Seq_Allele2", "dbSNP_RS", "dbSNP_Val_Status", "Genome_Change",
    "Protein_Change", "COSMIC_overlapping_mutations"
  )

  # add a major_alle column to the converted allelic imbalance dataframe
  data1$MajorAllele <- rep("dum", dim(data1)[1])
  for (i in 1:dim(data1)[1]) {
    if (data1[i, "TumorRef"] >= data1[i, "TumorAlt"]) {
      data1[i, 15] <- "TumorRef"
    } else {
      data1[i, 15] <- "TumorAlt"
    }
  }

  df1 <- data.frame(matrix("dummy", ncol = ncol(data2), nrow = nrow(data1)))
  colnames(df1) <- names(data2)
  colnames(df1) <- replace(names(df1), names(df1) == "Hugo_Symbol", "MajorAllele")
  df1$Chromosome <- data1$Chr
  df1$Start_Position <- data1$hg38_start
  df1$End_Position <- data1$hg38_end
  df1$Reference_Allele <- data1$Ref
  df1$Tumor_Seq_Allele2 <- data1$Mut
  df1$Chromosome <- as.character(df1$Chromosome)
  df1$Chromosome <- substr(df1$Chromosome, 4, nchar(df1$Chromosome))
  df1$MajorAllele <- data1$MajorAllele
  path <- getwd()
  # add one comment line in front of the data before saving it using `cat` command
  cat("##  Funcotator 4.2.0.0 | Date 20212019T112004 | Gencode 19 CANONICAL | Cosmic v84 | HGNC Nov302017 | Simple_Uniprot 2014_12 | dbSNP 9606_b150 \n", file = paste0(path, "/LOHmpileupInput.txt"), append = F)
  write.table(df1, file = paste0(path, "/LOHmpileupInput.txt"), append = TRUE, sep = "\t", row.names = FALSE, quote = FALSE)
}


#' Count the major allele and minor allele per mutant of each spot.
#'
#' @param x the path to the mpileup output of major- and minor- allele per mut of each spot
#' @param y a numeric value showing the sequential order of the spot name in the path x
#'
#' @return a dataframe of mutant major- and minor- allele counts. one row per mutant per spot
#' @export
#'
#' @examples \dontrun{
#' print("check user-guide")
#' }
SNPmajorAlleleCt <- function(x, y = 12) {
  if (!is.numeric(y)) {
    stop("y must be a numeric value!")
  }

  d1 <- read.delim(x, header = FALSE, sep = "\t", quote = "")
  d1[1, 9] <- "TumorRef"
  d1[1, 10] <- "TumorAlt"
  colnames(d1) <- d1[1, ]
  d1 <- d1[-1, ]

  df <- data.frame(matrix(nrow = nrow(d1), ncol = 6))
  CN <- c("Spot", "SNP", "Major_Count", "Minor_Count", "Total_Read_Count", "Fraction_Major")
  colnames(df) <- CN
  for (i in 1:nrow(d1)) {
    spot <- str_split_fixed(x, "/", 13)[y]
    SNP <- paste0("chr", d1[i, 2], ":", d1[i, 3], "-", d1[i, 5], "/", d1[i, 6])
    totalRDs <- as.numeric(d1[i, 9]) + as.numeric(d1[i, 10])
    majorAllele <- as.numeric(d1[i, d1[i, 1]])
    minorAllele <- totalRDs - majorAllele
    if (totalRDs > 0) {
      fractionMajor <- majorAllele / totalRDs
    } else {
      fractionMajor <- "NA"
    }

    df[i, 1] <- spot
    df[i, 2] <- SNP
    df[i, 3] <- majorAllele
    df[i, 4] <- minorAllele
    df[i, 5] <- totalRDs
    df[i, 6] <- fractionMajor
  }
  return(df)
}



#' Combine all spot major/minor allele counts into one table.
#'
#' @param files A list of pathes of spot mpileupOutput of major- and minor allele.
#' @param y The sequential order of the spot name is in the path
#'
#' @return a list of dataframe, the first is all spots SNPs major- and minor- count;
#' the other is the filtered version of SNPmajorAlleleCount table
#' @export
#'
#' @examples \dontrun{
#' print("check user guide")
#' }
lohMajorAlleleCt <- function(files, y = 12) {
  if (!is.numeric(y)) {
    stop("y must be a numeric value!")
  }

  Total_Read_Count <- NULL
  d3 <- parallel::mclapply(files, function(x) {
    data1 <- SNPmajorAlleleCt(x, y) # this takes 2 arguments; y=11 argument to make sure the spot name, need to test the previous code to figure out is 11 or 12 or ...
  }, mc.cores = 3)
  data <- data.frame(do.call("rbind", d3))
  colnames(data) <- c("Spot", "SNP", "Major_Count", "Minor_Count", "Total_Read_Count", "Fraction_Major")

  # filter out SNP with 0 reads
  datafilter <- data %>% filter(Total_Read_Count > 0)
  write.csv(data, file = "./SNPallMajorAlleleCount.csv", row.names = FALSE, quote = FALSE)
  write.csv(datafilter, file = "./SNPMajorAlleleCount.csv", row.names = FALSE, quote = FALSE)
}
