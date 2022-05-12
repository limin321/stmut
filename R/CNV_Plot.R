#' Get accumulated start position for cnr and cns file for arm-level CNV plotting.
#'
#' @param cnr the dataframe of cnr file from cnvkit output
#' @param seg the dataframe of cns file from cnvkit output
#' @export
#' @return A list of 2 dataframe, the new_ch for cnr file; and the newSeg for cns
#'
#'
#' @examples \dontrun{
#' cnr <- read.table(system.file("extdata/", "P6_T.deduplicated.realign.cnr",
#'   package = "stmut"
#' ), sep = "\t", header = TRUE)
#' cns <- read.table(system.file("extdata/", "P6_T.deduplicated.realign.cns",
#'   package = "stmut"
#' ), sep = "\t", header = TRUE)
#' df9 <- accStartCNR_CNS(cnr = cnr, seg = cns)
#' }
.accStartCNR_CNS <- function(cnr, seg) {
  chromosome <- NULL
  chr <- c(as.character(seq(22)), "X", "Y")[-1]
  new_ch <- cnr %>% filter(chromosome %in% "1")
  newSeg <- seg %>% filter(chromosome %in% "1")

  for (ch in chr) {
    if (ch == "X") {
      val <- max(new_ch$start)
      valEnd <- max(new_ch$end)
      df <- cnr %>% filter(chromosome %in% ch)
      df$start <- as.numeric(df$start) + as.numeric(val)
      df$end <- as.numeric(df$end) + as.numeric(valEnd)
      df$chromosome <- "23"
      new_ch <- rbind(new_ch, df)
      # change cns table gene start position
      sg <- seg %>% filter(chromosome %in% ch)
      sg$start <- as.numeric(sg$start) + as.numeric(val)
      sg$end <- as.numeric(sg$end) + as.numeric(valEnd)
      sg$chromosome <- "23"
      newSeg <- rbind(newSeg, sg)
    } else if (ch == "Y") {
      val <- max(new_ch$start)
      valEnd <- max(new_ch$end)
      df <- cnr %>% filter(chromosome %in% ch)
      df$start <- as.numeric(df$start) + as.numeric(val)
      df$end <- as.numeric(df$end) + as.numeric(valEnd)
      df$chromosome <- "24"
      new_ch <- rbind(new_ch, df)
      # change cns table gene start position
      sg <- seg %>% filter(chromosome %in% ch)
      sg$start <- as.numeric(sg$start) + as.numeric(val)
      sg$end <- as.numeric(sg$end) + as.numeric(valEnd)
      sg$chromosome <- "24"
      newSeg <- rbind(newSeg, sg)
    } else {
      val <- max(new_ch$start)
      valEnd <- max(new_ch$end)
      df <- cnr %>% filter(chromosome %in% ch)
      df$start <- as.numeric(df$start) + as.numeric(val)
      df$end <- as.numeric(df$end) + as.numeric(valEnd)
      new_ch <- rbind(new_ch, df)

      # change cns table gene start position
      sg <- seg %>% filter(chromosome %in% ch)
      sg$start <- as.numeric(sg$start) + as.numeric(val)
      sg$end <- as.numeric(sg$end) + as.numeric(valEnd)
      newSeg <- rbind(newSeg, sg)
    }
  }

  new_ch$chromosome <- as.numeric(new_ch$chromosome)
  newSeg$chromosome <- as.numeric(newSeg$chromosome)

  returnedList <- list(new_ch, newSeg)
  return(returnedList)
}

#' Arm-level copy number variations of bulk DNA-seq data.
#'
#' @param centmere the simplified 3-column hg38 centromere data
#' @param dcnr cnr file which is cnvkit output
#' @param dcns segment file which is cnvkit output
#'
#' @import ggplot2
#' @import dplyr
#'
#' @return An arm-level CNV plot, because the ref pool is male, so the plot also
#' has chrY data for female but the plot shows a deep deletion of chrY for female
#' @export
#'
#' @examples \dontrun{
#' centm <- read.csv(system.file("extdata/", "hg38_centromereSimple.bed", package = "stmut"),
#'   sep = "\t", header = FALSE
#' )
#' cnr <- read.table(system.file("extdata/", "P6_T.deduplicated.realign.cnr",
#'   package = "stmut"
#' ), sep = "\t", header = TRUE)
#' cns <- read.table(system.file("extdata/", "P6_T.deduplicated.realign.cns",
#'   package = "stmut"
#' ), sep = "\t", header = TRUE)
#' bulkCNVs(centmere = centm, dcnr = cnr, dcns = cns)
#' }
bulkCNVs <- function(centmere, dcnr, dcns) {
  if (!is.data.frame(centmere)) {
    stop("centmere must be a 3-column data frame.")
  }
  if (!is.data.frame(dcnr)) {
    stop("dcnr must be a dataframe. It is the output of cnvkit.")
  }
  if (!is.data.frame(dcns)) {
    stop("dcns must be a dataframe. It is the output of cnvkit")
  }

  chromosome <- NULL
  start <- NULL
  end <- NULL
  avg <- NULL
  colnames(centmere) <- c("chromosome", "start", "end")
  centmere$chromosome <- substr(centmere$chromosome, 4, nchar(centmere$chromosome))
  centmere[which(centmere$chromosome == "X"), 1] <- "23"
  centmere[which(centmere$chromosome == "Y"), 1] <- "24"
  centmere$chromosome <- as.numeric(centmere$chromosome)
  centmere <- centmere[order(centmere$chromosome), ]

  # get cnr data
  dcnr$chromosome <- substr(dcnr$chromosome, 4, nchar(dcnr$chromosome))
  dcns$chromosome <- substr(dcns$chromosome, 4, nchar(dcns$chromosome))

  r <- .accStartCNR_CNS(dcnr, dcns)
  df <- r[[1]]
  df2 <- r[[2]]

  # add previous chromosome max start position to centromere start position
  nums <- centmere$chromosome[-1]
  for (i in nums) {
    data <- df %>% filter(chromosome %in% (i - 1))
    val1 <- max(data$start)
    val2 <- max(data$end)
    centmere[which(centmere$chromosome == i), 2] <- centmere[which(centmere$chromosome == i), 2] + val1
    centmere[which(centmere$chromosome == i), 3] <- centmere[which(centmere$chromosome == i), 3] + val2
  }

  # calculate mean column to centmere data
  centmere <- centmere[order(centmere$chromosome), ]
  meanC <- rep(0, 24)
  centmere <- cbind(centmere, meanC)
  for (i in 1:nrow(centmere)) {
    centmere[i, 4] <- mean(c(centmere[i, 2], centmere[i, 3]))
  }

  x <- c(centmere$start, centmere$end)
  cenMean <- mean(x)

  # set chromosome names as factor for later use as X-axis.
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
    summarize(avg = max(start)) %>%
    pull(avg)

  ggplot(df) +
    annotate("point", x = df$start, y = df$log2, size = 0.00002, color = "grey", pch = 1, alpha = 0.2) +
    annotate("segment", x = df2$start, xend = df2$end, y = df2$log2, yend = df2$log2, color = "orange", size = 1) + # add segment layer to the plot.
    annotate("segment", x = 0, xend = max(df$start) * 1.02, y = 0, yend = 0, color = "black") +
    annotate("segment", x = t1, xend = t1, y = -1.5, yend = 1, color = "black") +
    geom_vline(data = centmere, mapping = aes(xintercept = end), linetype = 10, color = "black") + # centromere line
    scale_x_continuous(breaks = pos, labels = chrom, expand = c(0, 0)) +
    labs(x = "chromosome", y = "copy ratio (log2)") +
    ylim(c(-1.5, 1)) +
    theme_classic() + # use white color as background
    theme(axis.text.x = element_text(angle = 90))
}

#' Grouping spots from the same cluster to generate new spots to save spots with few reads.
#'
#' @param data A dataframe of spots to be grouped
#' @param Treads The min reads of the a new spot, by default Treads=5000
#' @param NumSpt NumSpt*NumSpt nearest spots will be included to determine which
#' spots to grouped together. By default, NumSpt=6
#'
#' @importFrom raster pointDistance
#'
#' @return A dataframe of new spots
#' @export
#'
#' @examples \dontrun{
#' clst1 <- read.csv(system.file("extdata/", "c1.csv", package = "stmut"), header = TRUE)
#' df10 <- groupSpots(clst1, Treads = 5000, NumSpt = 8)
#' }
groupSpots <- function(data, Treads = 5000, NumSpt = 8) {
  if (!is.data.frame(data)) {
    flog.error("data should be a dataframe and have 6 columns: 'barcode', 'cluster', 'spot', 'TotalRDs', 'array_row', 'array_col'.")
    stop("data must be a data frame.")
  }
  if (!is.numeric(Treads)) {
    stop("Treads must be numeric!")
  }
  if (!is.numeric(NumSpt)) {
    stop("NumSpt must be numeric!")
  }

  array_row <- NULL
  barcode <- NULL
  # to store newly combined group
  spt2 <- NumSpt * NumSpt
  d1 <- data # save to d1 for later use in merge(), the original data variable will be consumed.

  groups <- data.frame(
    barcode = character(),
    groups = character(),
    stringsAsFactors = FALSE
  )
  n <- 1 # initialize new group number
  if (sum(data$TotalRDs) < Treads) { # if total reads less than 5000, put all spots in one group
    g1 <- data.frame(barcode = data$barcode, groups = rep(paste0("group", n), dim(data)[1]))
    groups <- rbind(groups, g1)
    n <- n + 1
  } else { # if the total reads > given number, grouping spots

    for (j in 1:nrow(data)) { # loop each row to group spots
      if (dim(data)[1] == 0) {
        break
      } else if (dim(data)[1] == 1) {
        g1 <- data.frame(barcode = data$barcode, groups = paste0("group", n))
        groups <- rbind(groups, g1)
        break
      } else {
        newdata <- data.frame(
          barcode = character(), # empty dataframe to store tempor 36 spots for distance calculation
          reads = integer(),
          array_row = integer(),
          array_col = integer(),
          stringsAsFactors = FALSE
        )

        # first do QC to only keep 36 points for distance calculation by only keep 6 continuous x-axis spots, 6 continuous y-axis
        b <- unique(data$array_row)
        for (i in 1:length(b)) { # loop each row/x-axis
          nam <- data %>% filter(array_row %in% b[i]) # dataframe of the same row.
          a <- dim(nam)[1] # number of spots- spots having the same x-axis but different y-axis

          # each x position, only keep 6 spots with different y positions
          if (a <= NumSpt) {
            newdata <- rbind(newdata, nam)
          } else { # only keep 6 continuous spots
            nam <- nam[-c((NumSpt + 1):a), ] # remove spots too far away
            newdata <- rbind(newdata, nam)
          }
        }
        ########################
        # subset the first 36 spots
        if (nrow(newdata) >= spt2 & sum(newdata$TotalRDs) >= Treads) {
          newc4 <- newdata[c(1:spt2), ]
        } else if (nrow(newdata) >= spt2 & sum(newdata$TotalRDs) <= Treads) {
          stop("Some grouped spots have too few reads.Please assign a bigger number to variable NumSpt.")
          break
        } else if (nrow(newdata) < spt2 & sum(data$TotalRDs) >= Treads) {
          newc4 <- newdata
        } else {
          newc4 <- newdata
        }
        ############################
        dall <- c(0) # store distance of the 35 distance, because distance of the same spot is 0.
        p1 <- newc4[1, c("array_row", "array_col")]
        for (i in 2:nrow(newc4)) {
          p2 <- newc4[i, c("array_row", "array_col")]
          a <- pointDistance(p1, p2, lonlat = FALSE)
          dall <- append(dall, a)
        }
        newc4 <- cbind(newc4, dall, row.names = NULL)
        newc5 <- newc4[order(newc4$dall), ] # sort the df based on distance to the first spot
        # get fewest spots whose sum of reads is greater than 25000

        if (sum(newc5$TotalRDs) < Treads) {
          g1 <- data.frame(barcode = newc5$barcode, groups = rep(paste0("group", n), dim(newc5)[1])) # new group
          groups <- rbind(groups, g1)
          data <- data %>% filter(!barcode %in% g1$barcode) # remove the grouped spots and start loop again.
          n <- n + 1
        } else {
          for (i in 1:nrow(newc5)) {
            reads <- newc5$TotalRDs
            readT <- sum(reads[1:i])

            if (readT >= Treads) {
              k <- i # to filter each group
              break
            }
          }
          rmCode <- newc5$barcode[1:k] # spot barcodes grouped as a new spot
          g1 <- data.frame(barcode = rmCode, groups = rep(paste0("group", n), length(rmCode))) # new group
          groups <- rbind(groups, g1)

          data <- data %>% filter(!barcode %in% rmCode) # remove the grouped spots and start loop again.
          n <- n + 1
        }
      }
    }
  } # loop each row ends; else
  ##########
  final <- merge(d1, groups, by = "barcode")
  return(final)
}



# # version 2 of groupSpots function, taking 2 argument. For future use if the current
# # groupSpots keeps giving bugs.
# # The strategy of this version is easier,1) calculating the distance between
# # one spot and the rest spots; 2) sort the spots based on the distance;
# # 3) take the first N spots which the sum of spots-reads > 5k.
# groupSpots <- function(clst, Treads = 5000){
#   if (!is.data.frame(data)) {
#     flog.error("data should be a dataframe and have 6 columns: 'barcode', 'cluster', 'spot', 'TotalRDs', 'array_row', 'array_col'.")
#     stop("data must be a data frame.")
#   }
#   if (!is.numeric(Treads)) {
#     stop("Treads must be numeric!")
#   }
#
#   # initialize groups info
#   groups <- data.frame(
#     barcode = character(),
#     groups = character(),
#     stringsAsFactors = FALSE
#   )
#   n = 1 # initialize group numbers
#
#   for (x in 1:nrow(clst)){
#     clst <- clst[with(clst, order(array_row,array_col)),]
#
#     if (nrow(clst) == 0) {
#       break
#     } else if (nrow(clst) == 1) {
#       g1 <- data.frame(barcode = clst$barcode, groups = paste0("group", n))
#       groups <- rbind(groups, g1)
#
#       clst <- clst %>% filter(!barcode %in% g1$barcode) # remove the grouped spots and start loop again.
#       n <- n + 1
#       break
#     }else{
#       # calculate distance of the first and the rest
#       dall <- c(0) # store distance of the 35 distance, because distance of the same spot is 0.
#       p1 <- clst[1,c("array_row", "array_col")]
#       for (i in 2:nrow(clst)) {
#         p2 <- clst[i, c("array_row", "array_col")]
#         a <- pointDistance(p1, p2, lonlat = FALSE)
#         dall <- append(dall, a)
#       }
#
#       clst$dall <- dall
#       clst <- clst[order(dall),] # sort the data based on the distance from smallest to largest
#
#       # get the number of spots that meet your Reads requirements, by default is 5k reads/spot
#       for (i in 1:nrow(clst)) {
#         reads <- clst$TotalRDs
#         readT <- sum(reads[1:i])
#
#         if (readT >= Treads) {
#           k <- i # to filter each group, cutoff to keep rows in each group
#           break
#         }else{
#           k = nrow(clst)
#         }
#       }
#
#       rmCode <- clst$barcode[1:k] # spot barcodes grouped as a new spot
#       g1 <- data.frame("barcode" = rmCode, "groups" = rep(paste0("group", n), length(rmCode))) # new group
#       groups <- rbind(groups, g1)
#
#       clst <- clst %>% filter(!barcode %in% rmCode) # remove the grouped spots and start loop again.
#       n <- n + 1
#     }
#   }
#
#
# }








#' Generating new_spots_feature_bc.csv after calling groupSpots function.
#'
#' @param df1 dataframe of spots to be grouped
#' @param Treads min reads of new spot, default is 5000
#' @param NumSpt number of spots to be considered, default is 8
#' @param data Original filtered_feature_bc.csv
#'
#' @importFrom stats na.omit
#'
#' @return A list of 2 dataframe and a vector of representative barcodes of new spots.
#' One is grouping spots with grouping-name; the other is newspots(grouped_spots + non_group_spots) feature barcode dataframe
#' @export
#'
#' @examples \dontrun{
#' clst1 <- read.csv(system.file("extdata/", "c1.csv", package = "stmut"), header = FALSE)
#' c1BC <- read.csv(system.file("extdata/", "c1_feature_bc.csv", package = "stmut"), header = FALSE)
#' df11 <- newSptBC(df1 = clst1, Treads = 5000, NumSpt = 8, data = c1BC)
#' }
newSptBC <- function(df1, Treads = 5000, NumSpt = 8, data) {
  if (!is.data.frame(data)) {
    flog.error("data should be a dataframe and have 6 columns: 'barcode', 'cluster', 'spot', 'TotalRDs', 'array_row', 'array_col'.")
    stop("data must be a data frame.")
  }
  if (!is.data.frame(df1)) {
    stop("df1 must be a dataframe.")
  }
  if (!is.numeric(Treads)) {
    stop("Treads must be numeric!")
  }
  if (!is.numeric(NumSpt)) {
    stop("NumSpt must be numeric!")
  }

  cluster <- NULL
  uniqueName <- c()
  fdf <- data.frame(matrix(ncol = 6, nrow = 0)) # the final fdf has same dimension as df1 but adding group-name within each cluster
  name <- c("barcode", "array_row", "array_col", "cluster", "TotalRDs", "groups")
  colnames(fdf) <- name
  n <- length(unique(df1$cluster)) # number of clusters
  n <- unique(df1$cluster)
  for (i in n) { # loop clusters
    clst <- df1 %>% filter(cluster == i)
    clst <- clst[with(clst, order(array_row,array_col)),] # sort the data first by row, then by col.
    df <- groupSpots(data = clst, Treads = Treads, NumSpt = NumSpt) # call R function
    fdf <- rbind(fdf, df)
  }
  fdf1 <- fdf #

  # generate new new_feature_bc.csv
  fdf$barcode <- str_replace(fdf$barcode, "-1", ".1")
  keep <- data[, !(names(data) %in% c(fdf$barcode))] # not included for spotcombining
  combSpts <- data[, c(fdf$barcode)] # spots for grouping

  k <- 0 # number of new spots
  for (cl in unique(fdf$cluster)) {
    fdf2 <- fdf %>% filter(cluster == cl)
    len <- length(unique(fdf2$groups))
    k <- k + len
    for (gp in unique(fdf2$groups)) {
      d1 <- (fdf %>% filter(cluster == cl) %>% filter(groups == gp))

      if (dim(d1)[1] == 1) {
        d1 <- d1$barcode
        df1 <- combSpts[, c(d1)]
        nam1 <- d1[1]
        newcol <- df1
      } else {
        d1 <- d1$barcode
        df1 <- combSpts[, c(d1)]
        nam1 <- d1[1]
        newcol <- rowSums(df1)
      }
      keep <- cbind(keep, newcol)
      name <- colnames(keep)
      colnames(keep) <- replace(name, name == "newcol", nam1)
      uniqueName <- append(uniqueName, nam1)
    }
  }
  keep <- data.frame(keep)
  returnList <- list(fdf1, keep, uniqueName)
  return(returnList)
}

#' Calculate Chromosome Arm CNVs Weighted Median to Represent the CNVs of that arm.
#'
#' @param data The cnr dataframe of each spot.
#' @param centmere The centromere dataframe of hg38.
#'
#' @importFrom matrixStats weightedMedian
#'
#' @return A list of 2 dataframe, the weighted cnr file and a summary file counting
#' how many genes in each arm.
#' @export
#'
#' @examples \dontrun{
#' cnr <- read.table(system.file("extdata/", "spot1_rep1.cnr", package = "stmut"), header = TRUE)
#' centm <- read.csv(system.file("extdata/", "hg38_centromereSimple.bed", package = "stmut"),
#'   sep = "\t", header = FALSE
#' )
#' df5 <- wtArmMedianOne(data = cnr, centmere = centm)
#' }
wtArmMedianOne <- function(data, centmere) {
  if (!is.data.frame(data)) {
    stop("data must be a dataframe!")
  }
  if (!is.data.frame(centmere)) {
    stop("centmere must be a dataframe!")
  }

  chromosome <- NULL
  # when calculate rolling median, rolling through each arm within of each chromosome.
  dataNew <- data %>% filter(!chromosome %in% c("MT", "Y")) # remove MT, Y chromosome
  dataNew[which(dataNew$chromosome == "X"), 1] <- "23" # replace X with "23"

  chr <- unique(dataNew[, which(colnames(data) == "chromosome")])
  colnames(centmere) <- c("chromosome", "start", "end")
  centmere$chromosome <- substr(centmere$chromosome, 4, nchar(centmere$chromosome))
  centmere[which(centmere$chromosome == "X"), 1] <- "23"
  centmere[which(centmere$chromosome == "Y"), 1] <- "24"
  centmere$chromosome <- as.numeric(centmere$chromosome)
  centmere <- centmere[order(centmere$chromosome), ] # adjust centmere data to be in chromosome order
  centmere$chromosome <- as.character(centmere$chromosome)

  qArm <- c() # store only qArms

  df <- data.frame(
    characters = character(),
    Integers = integer(),
    Integers = integer(),
    Characters = character(),
    Doubles = double(),
    Doubles = double(),
    Doubles = double(),
    Integers = integer(),
    Doubles = double(),
    stringsAsFactors = FALSE
  )
  dfArms <- data.frame(matrix(nrow = length(chr), ncol = 5)) # store number genes of each arm for cdt sorting use
  colnames(dfArms) <- c("chromosome", "p_Genes", "q_Genes", "pArmEnds", "CM_row_pos")

  for (ch in chr) { # first loop over chromosome
    data1 <- dataNew %>% filter(chromosome %in% ch) # chr1
    ctm1 <- centmere %>% filter(chromosome %in% ch) # centromere of chr1
    nch <- length(data1[, which(colnames(data1) == "log2")])
    # to classify chromosomes into 2 groups, with and withour p-arms
    for (pos in data1$start) { # second loop over start pos of each chromosome to get break-points for centromere
      if (pos > ctm1$start) {
        val <- which(data1$start == pos) # row pos of p arm ends,q arm starts
        if (val == 1) {
          qArm <- append(qArm, ch)
        }
        dfArms[ch, "chromosome"] <- ch
        dfArms[ch, "p_Genes"] <- val - 1
        dfArms[ch, "q_Genes"] <- nrow(data1) - val + 1
        dfArms[ch, "pArmEnds"] <- pos
        dfArms[ch, "CM_row_pos"] <- val - 1
        break
      }
    } # end of second for loop

    if (ch %in% qArm) { # weighted-median for chromosomes without p-arms
      dat <- data1[c(1:nch), which(colnames(data1) == "log2")]
      w <- data1[c(1:nch), which(colnames(data1) == "weight")]
      m <- weightedMedian(dat, w)
      data1[which(data1$chromosome == ch), which(colnames(data1) == "log2")] <- m
      df <- rbind(df, data1)
    } else { # chromosome with both arms
      dat1 <- data1[c(1:(val - 1)), which(colnames(data1) == "log2")]
      w1 <- data1[c(1:(val - 1)), which(colnames(data1) == "weight")]
      m1 <- weightedMedian(dat1, w1)
      data1[c(1:(val - 1)), which(colnames(data1) == "log2")] <- m1
      dat2 <- data1[c(val:nch), which(colnames(data1) == "log2")]
      w2 <- data1[c(val:nch), which(colnames(data1) == "weight")]
      m2 <- weightedMedian(dat2, w2)
      data1[c(val:nch), which(colnames(data1) == "log2")] <- m2
      df <- rbind(df, data1)
    }
  }

  returnList <- list(df, dfArms)
  return(returnList)
} # function ending

#' Reformat Gene Summary of Each Arm for cdt_filt_sort function use in the next step.
#'
#' @param cdt the filtered cdt file.
#' @param data A summary file generated by wtArmMedianOne function.
#'
#' @return A dataframe prepared for sorting the cdt based on bulk CNVs information.
#' @export
#'
#' @examples \dontrun{
#' cdt <- read.table(system.file("extdata/", "cdt.cdt", package = "stmut"), header = TRUE)
#' data3 <- read.table(system.file("extdata/", "summary.txt", package = "stmut"),
#'   sep = "\t", header = TRUE
#' )
#' df6 <- CtArmGenes(cdt = cdt, data = data3)
#' }
CtArmGenes <- function(cdt, data) {
  if (!is.data.frame(cdt)) {
    stop("cdt must be data frame!")
  }
  if (!is.data.frame(data)) {
    stop("data must be a data frame!")
  }

  arms <- NULL
  chr <- NULL
  # reference hg38 centromere bed file
  chromosome <- seq(24)
  start <- c(122026459, 92188145, 90772458, 49712061, 46485900, 58553888, 58169653, 44033744, 43389635, 39686682, 51078348, 34769407, 16000000, 16000000, 17083673, 36311158, 22813679, 15460899, 24498980, 26436232, 10864560, 12954788, 58605579, 10316944)
  end <- c(124932724, 94090557, 93655574, 51743951, 50059807, 59829934, 61528020, 45877265, 45518558, 41593521, 54425074, 37185252, 18051248, 18173523, 19725254, 38265669, 26616164, 20861206, 27190874, 30038348, 12915808, 15054318, 62412542, 10544039)
  centmere <- data.frame(chromosome, start, end)
  centmere$chromosome <- as.character(centmere$chromosome)

  # add arm column to input cdt file
  df1 <- data.frame(cdt$NAME)
  colnames(df1) <- c("NAME")
  df1[c("chr", "pos", "gene")] <- str_split_fixed(df1$NAME, ":", 3)
  df1[c("start", "end")] <- str_split_fixed(df1$pos, "-", 2)
  df1 <- df1[, c("NAME", "chr", "start", "end")]
  df1$start <- as.numeric(df1$start)
  df1$end <- as.numeric(df1$end)
  df1$arms <- rep(0, nrow(df1))

  d0 <- data.frame(matrix(nrow = 0, ncol = 5))
  colnames(d0) <- c("NAME", "chr", "start", "end", "arms")
  for (ch in centmere$chromosome) {
    d1 <- df1 %>% filter(chr == ch)
    cent <- centmere %>% filter(chromosome == ch)
    start <- cent[1, "start"]
    end <- cent[1, "end"]
    for (loc in d1$end) {
      if (loc < start) {
        idx <- which(d1$end == loc)
        d1[idx, "arms"] <- paste0(ch, "p")
      } else {
        idx <- which(d1$end == loc)
        d1[idx, "arms"] <- paste0(ch, "q")
      }
    }
    d0 <- rbind(d0, d1)
  }
  d1 <- d0[, c(1, 5)]
  cdt <- merge(d1, cdt, by = "NAME", sort = FALSE)

  # Create each arm, genes summary file
  dfgene <- data.frame(matrix(nrow = 0, ncol = 4)) # store arm, chromosome,number of genes in each arm, arm_start_pos (gene_row)
  colnames(dfgene) <- c("arms", "genes", "chromosome", "gene_row")

  i <- 1 # count how many arms
  for (arm in unique(cdt$arms)) {
    d0 <- cdt %>% filter(arms == arm)
    ch <- substr(arm, 1, nchar(arm) - 1)
    numGene <- dim(d0)[1]
    rowN <- min(which(cdt$arms == arm))
    dfgene[i, "arms"] <- arm
    dfgene[i, "chromosome"] <- ch
    dfgene[i, "genes"] <- numGene
    dfgene[i, "gene_row"] <- rowN
    i <- i + 1
  }
  return(dfgene)
}


#' Rank the spots by their similarity to the DNA-seq copy number alterations.
#'
#' @param cdt cdt file to be ranked
#' @param genes Number of genes in each arm
#' @param gainLoss Gain or Loss info of each arm, 1 represent gain, -1 represent loss.
#' @param rs Representative row position of each arm.
#'
#' @return A spots dataframe sorted by CNVs
#' @export
#'
#' @examples \dontrun{
#' cdt <- read.table(system.file("extdata/", "cdt.cdt", package = "stmut"), header = TRUE)
#' data4 <- read.csv(system.file("extdata/", "CtArmGenSummary.csv", package = "stmut"), header = TRUE)
#' arm <- c(
#'   "1p", "3p", "3q", "4q", "5q", "8q", "9q", "10p", "10q", "11q", "13p", "13q",
#'   "20p", "20q", "21q", "14q", "17q"
#' )
#' d4 <- data4 %>% filter(arms %in% arm)
#' genes <- d4[, "genes"]
#' rs <- d4[, "gene_row"]
#' gainLoss <- c(1, -1, 1, -1, -1, 1, 1, -1, -1, 1, -1, -1, 1, 1, -1, 1, 1)
#' df7 <- cdt_filt_sort(cdt = cdt, genes = genes, gainLoss = gainLoss, rs = rs)
#' }
cdt_filt_sort <- function(cdt, genes, gainLoss, rs) {
  if (!is.data.frame(cdt)) {
    stop("cdt must be a data frame!")
  }

  name <- cdt[, c(1, 2)]
  cdt <- cdt[, -c(1, 2)]

  # convert cdt values to be numeric
  cdt[] <- lapply(cdt, as.numeric) # convert the dataframe to numeric

  wt <- c() # normalize number of genes by divide the max number of genes
  for (n in genes) {
    weit <- n / max(genes)
    wt <- append(wt, weit)
  }
  newR <- c() # newR created for sorting the dataframe
  nc <- dim(cdt)[2]
  for (i in c(1:nc)) {
    val <- 0
    for (k in 1:length(genes)) {
      val1 <- gainLoss[k] * (cdt[rs[k], i]) * wt[k] # val <- (-1)*(cdt[rs[1],i])*wt[1]+cdt[rs[2],i]*wt[2]+(-1)*(cdt[rs[3],i])*wt[3] # the sum of weighted values, 3q and 8q has CN gains, 3p and 8p has CN loss, so val = 3q*wt - 3p*wt +8q*wt-8p*wt
      val <- val + val1
    }
    newR <- append(newR, val)
  }

  cdt1 <- rbind(cdt, newR)
  cdt1 <- cdt1[, order(cdt1[nrow(cdt1), ], decreasing = FALSE)] # reorder by the added row, order() has warning, but works.
  cdt2 <- cdt1[-nrow(cdt1), ] # remove the added row
  cdt3 <- cbind(name, cdt2)
  return(cdt3)
}

#' Create grouped spots summary files.
#'
#' @param df1 spots summary file generated by sptBClstRds()
#' @param Treads minimum reads for new spot
#' @param NumSpt spots range for grouping
#' @param data filtered_feature_bc.csv from spaceranger.
#' @param locpath path to tissue_positions_list.csv file from spaceranger.
#'
#' @importFrom utils write.csv
#'
#' @return generate CNVs/spotIndex/, CNVs/txt/ folder; and grpGraph_based.csv,spotGroupingSummary.csv,
#' newSpotBcRds.csv, newSpotSummary.csv, tissue_positions_list.csv, newSptsGenExp.csv files.
#' @export
#'
#' @examples \dontrun{
#' print("check user guide")
#' }
spotsGrouped <- function(df1, Treads = 5000, NumSpt = 8, data, locpath) {
  V1 <- NULL
  barcode <- NULL
  dir.create("./CNVs/txt/", recursive = TRUE)
  dir.create("./CNVs/spotIndex/", recursive = TRUE)
  path <- getwd()
  path1 <- paste0(path, "/CNVs/spotIndex/")
  path2 <- paste0(path, "/CNVs/txt/")

  df <- newSptBC(df1, Treads, NumSpt, data)
  r21 <- df[[1]]
  r22 <- df[[2]]
  r23 <- df[[3]]
  if (nrow(df1) == nrow(r21)){
    write.csv(r21, file = paste0(path, "/CNVs/spotGroupingSummary.csv"), row.names = FALSE)
    write.csv(r22, file = paste0(path, "/CNVs/newSptsGenExp.csv"), row.names = FALSE)
  }else{
    stop("Grouping has issue. Please double check the parameters and inputs.")
  }


  for (i in 1:length(r23)) {
    name <- paste0("spot", i)
    bc <- r23[i]
    bc <- str_replace(bc, ".1", "-1")
    write.table(bc, file = paste0(path1, name, ".txt"), col.names = FALSE, quote = FALSE, row.names = FALSE) # new spot index file
    expSpt <- r22[, c(1, (i + 1))]
    write.table(expSpt, file = paste0(path2, name, ".txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE) # new spot gene expr txt file.
    i <- i + 1
  }

  v1 <- paste0("spot", seq(length(r23)))
  r21$barcode <- str_replace(r21$barcode, "-1", ".1")
  t2 <- r21 %>% filter(barcode %in% r23) # obtain representative cluster name
  t2$cluster <- gsub(" ", "", t2$cluster)
  data1 <- data.frame(r23, v1)
  colnames(data1) <- c("barcode", "spot")

  data2 <- merge(t2, data1, by = "barcode")
  data2$barcode <- str_replace(data2$barcode, ".1", "-1")
  graphbased <- data2[, c(1, 2)]
  write.csv(graphbased, file = paste0(path, "/CNVs/grpGraph_based.csv"), row.names = FALSE)

  loc <- read.csv(locpath, header = FALSE)
  loc$V1 <- str_replace(loc$V1, "-1", ".1")
  loc1 <- loc %>% filter(V1 %in% r23)
  loc1$V1 <- str_replace(loc1$V1, ".1", "-1")
  write.table(loc1, file = paste0(path, "/CNVs/tissue_positions_list.csv"), row.names = FALSE, col.names = FALSE, sep = ",", quote = FALSE)

  ## convert r22 to numeric
  # r22[] <- lapply(r22, as.numeric) # this will raise a warning bc the first column is ensembl gene names which cannot be converted to numeric.
  r22[] <- suppressWarnings(lapply(r22, as.numeric)) # suppress waring.


  files <- list.files(path = paste0(path, "/CNVs/spotIndex"), ".txt", full.names = TRUE, recursive = FALSE)
  p1 <- paste0(path, "/CNVs/tissue_positions_list.csv")
  df1 <- sptBClstRds(files = files, data1 = r22, data2 = graphbased, path = p1)
  write.csv(df1[[1]], file = paste0(path, "/CNVs/newSpotBcRds.csv"), row.names = FALSE)
  write.csv(df1[[2]], file = paste0(path, "/CNVs/newSpotSummary.csv"), row.names = FALSE)
}


#' Generate a cdt file of sorted CNVs data of spots based on DNAseq CNV signals.
#'
#' @param data1 new spots summary data
#' @param cdt CNVs data to be sorted.
#' @param data2 arm-based gene summary file, output of wtArmMedianOne()
#' @param arm a vector of arms with CNVs signals
#' @param gainLoss CNVs signals corresponding to arms, -1 is loss, +1 is gain
#' @param tumorclusters customer-defined tumor clusters.
#'
#' @importFrom stats median
#'
#' @return a cdt file for JavaTree vies or R visualization.
#' @export
#'
#' @examples \dontrun{
#' print("check user guide")
#' }
cnvsorted <- function(data1, cdt, data2, arm, gainLoss, tumorclusters) {
  stopifnot(is.data.frame(data1), is.data.frame(cdt), is.data.frame(data2))

  if (length(arm) != length(gainLoss)) {
    stop("arm and gainLoss should have the same length!")
  }

  arms <- NULL
  cdt1 <- cdt[-c(1:2), -c(1:4)]
  c1 <- data.frame(t(cdt1))
  c2 <- cbind(rownames(c1), c1)
  name <- names(c2)
  colnames(c2) <- replace(name, name == "rownames(c1)", "sptRep")
  c3 <- merge(data1, c2, by = "sptRep")

  c3 <- c3[order(c3$rep, c3$cluster, c3$TotalRDs), ] # sort the c3 by rep, cluster, TotalRDs
  rownames(c3) <- c3$final
  c4 <- c3[, -c(1:5)]
  c4[] <- lapply(c4, as.numeric)

  a <- ncol(c4) # how many rows in total
  c5 <- c4[rowSums(c4 == 0) / a < 0.05, ] # remove spots which has more 5% of genes/values are 0s.

  c6 <- c5[!grepl(tumorclusters, rownames(c5)), ] # remove tumor spots
  meds <- apply(c6, 2, median) # median of non-tumor values as ref
  c7 <- data.frame(apply(c5, 1, function(x) {
    x - meds
  })) # each spot - ref
  rows2 <- cdt[-c(1:2), c(2:3)]
  c8 <- cbind(rows2, c7)

  ArmGenSumm <- CtArmGenes(cdt = c8, data = data2)
  df1 <- ArmGenSumm %>% filter(arms %in% arm) # sort the spots based on the CNV signal found in the bulk data
  genes <- df1[, "genes"]
  rs <- df1[, "gene_row"]
  df <- cdt_filt_sort(cdt = c8, genes = genes, gainLoss = gainLoss, rs = rs)
  write.table(df, file = paste0(getwd(), "/CNVs_sorted.cdt"), sep = "\t", row.names = FALSE, quote = FALSE)
}



#' Generate filtered, centered cdt file for spot CNVs or cell CNVs.
#'
#' @param data1 A two-column csv file; the 2 column names are "spotName", "annotate".
#' In the "annotate" column, spots which are used as reference should be annotated as "normal". These spots are used to center the data.
#' @param cdt Original cdt file to be processed.
#' @param data2 Arm-based gene summary file, output of wtArmMedianOne()
#' @param data3 The spotSummary file, generated by sptBClstRds()
#' @param arm Optional. If you don't have DNAseq data. arm is not required. A vector of arms with CNVs signals
#' @param gainLoss Optional. If you don't have DNAseq data. gainLoss is not required. CNVs signals corresponding to arms, -1 is loss, +1 is gain
#'
#' @return A centered cdt file as final CNVs of spots or cells.
#' @export
#'

cnvs <- function(data1, cdt, data2, data3, arm=c(), gainLoss=c()) {
  stopifnot(is.data.frame(data1), is.data.frame(cdt), is.data.frame(data2))

  if (length(arm) != length(gainLoss)) {
    stop("arm and gainLoss should have the same length!")
  }

  evens <- function(x) subset(x, x %% 2 == 0) # only pick even numbers

  arms <- NULL
  spot <- NULL
  chr <- NULL
  cluster <- NULL
  cdt1 <- cdt[-c(1:2), -c(1:4)]

  # cap CNV values to [-1, 1]
  cdt1[] <- lapply(cdt1, as.numeric)
  cdt1[cdt1 > 1] <- 1
  cdt1[cdt1 < -1] <- -1

  # merge cdt data with annotation data
  c1 <- data.frame(t(cdt1))
  c2 <- cbind(rownames(c1), c1)
  name <- names(c2)
  colnames(c2) <- replace(name, name == "rownames(c1)", "spotName")
  c3 <- merge(data1, c2, by = "spotName")
  rownames(c3) <- c3[, "spotName"]

  # remove spots which has more 5% of genes/values are 0s.
  c4 <- c3[, -c(1, 2)]
  c4[] <- lapply(c4, as.numeric)
  c4 <- c4[rowSums(c4 == 0) / (ncol(c4) - 2) < 0.05, ]

  # calculate non-tumor/normal spots median
  nonT <- data1[data1$annotate == "normal", ]$spotName
  c5 <- c4[rownames(c4) %in% nonT, ]
  meds <- apply(c5, 2, median) # median of each gene/col of non-tumor/normal spots
  c6 <- data.frame(apply(c4, 1, function(x) {
    x - meds
  })) # each spot/row - normal_median
  rows2 <- cdt[-c(1:2), c(2:3)]
  c7 <- cbind(rows2, c6)

  # if NO DNA-seq data, NO NEED to run the sorting based on the bulk-CNV score. Instead, sort the data based on cluster and totalRDs.

  if (length(arm) ==0){
    dat2 <- rbind(names(c6), c6) # add spot-name as first row
    cdt2 <- data.frame(t(dat2))

    #change the first colname to be "spot".
    name <- colnames(cdt2)
    colnames(cdt2) <- replace(name, name == "X1","spot" )

    # spot summary file
    data3$cluster <- stringr::str_remove_all(data3$cluster," ")
    sum1 <- data3[, c("cluster", "spot", "TotalRDs")] %>% filter(spot %in% cdt2$spot)
    cdt3 <- merge(sum1, cdt2, by = "spot")
    cdt3$TotalRDs <- as.numeric(cdt3$TotalRDs)
    cdt4 <- cdt3[order(cdt3[,2], -cdt3[,3]),] # sort first by cluster column, then by TotalRDs column-decreasing.
    cdt5 <- data.frame(t(cdt4))
    colnames(cdt5) <- cdt5[1,]
    cdt6 <- cdt5[-c(1:3),]
    df1 <- cdt[-c(1:2),2:3]
    cdt7 <- cbind(df1, cdt6)
    write.table(cdt7, file = paste0(getwd(),"/cluster_TotalRDs.cdt"), sep = "\t", row.names = FALSE, col.names = TRUE)

    # 1.1 generating chromosome heatmap of sorted cluster_TotalRDs.cdt for annotating chr locations of cdt file.
    df2 <- df1
    df2[c("chr", "loc")] <- str_split_fixed(df2$NAME, ":",2)

    ### Generate chromosome-template cdt file----- For troubleshooting ChromosomeTemplate.png
    #   df0 <- df2
    #   df0$chr <- as.numeric(df0$chr)
    #   df0$value <- rep(0, nrow(df0))
    #   nums <- evens(1:length(unique(df0$chr)))
    #   chrs <- unique(df0$chr)[nums] # only pick even number clusters
    # for (i in 1:nrow(df0)){
    #   if(df0[i,"chr"] %in% chrs){
    #     df0[i,"value"] <- 1
    #   }
    # }
    #   df0 <- df0[,c("CLID", "NAME", "value")]
    # write.table(df0, file =paste0(getwd(), "/ChromosomeTemplate.cdt"), sep = "\t", row.names = FALSE, col.names = TRUE)

    df2 <- dplyr::count(df2,chr)
    df2$seq <- rep(1, nrow(df2))
    df2 <- df2[order(as.numeric(df2$chr), decreasing = FALSE),]
    df2$chr <- as.numeric(df2$chr)

    if (nrow(df2) %% 2 == 0){
      a = nrow(df2)/2
      colors <- rep(c("#5A5A5A", "#CCCCCC"), a)
    } else{
      a = (nrow(df2) - 1)/2
      colors <- c("#5A5A5A", rep(c("#CCCCCC","#5A5A5A"), a))
    }

    ####
    p1 <- paste0(getwd(), "/ChromosomeTemplate.png")
    ggplot(df2, aes(x = seq, y=n, fill = factor(chr))) +
      geom_bar(position = "stack", stat = "identity") +
      scale_fill_manual(values = colors)+
      theme_bw() +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()
      )
    ggsave(p1)

    # 1.2 generating cluster heatmaps of sorted cluster_TotalRDs.cdt
    NAME <- data.frame(t(cdt5[2,]))

    # ## clusterTemplate.cdt ----- For troubleshooting clusterTemplate image
    # image <- paste0("IMAGE:", seq(nrow(NAME))-1)
    # NAME$image <- image
    # NAME <- NAME[,c("image", "cluster")]
    # NAME$value <- rep(0, nrow(NAME))
    #
    # nums <- evens(1:length(unique(NAME$cluster)))
    # clst <- unique(NAME$cluster)[nums] # only pick even number clusters
    # for (i in 1:nrow(NAME)){
    #   if(NAME[i,"cluster"] %in% clst){
    #     NAME[i,"value"] <- 1
    #   }
    # }
    # write.table(NAME, file =paste0(getwd(), "/clusterTemplate.cdt"), sep = "\t", row.names = FALSE, col.names = TRUE)

    ## cluster1...clusterN heatmap
    dat1 <- dplyr::count(NAME, cluster)
    dat1$seq <- rep(1, nrow(dat1))
    if (nrow(dat1) %% 2 == 0){
      a = nrow(dat1)/2
      colors <- rep(c("#5A5A5A", "#CCCCCC"), a)
    } else{
      a = (nrow(dat1) - 1)/2
      colors <- c("#5A5A5A", rep(c("#CCCCCC","#5A5A5A"), a))
    }

    ##
    p2 <- paste0(getwd(), "/ClusterTemplate.png")
    ggplot(dat1, aes(x = seq, y=n, fill = cluster)) +
      geom_bar(position = "stack", stat = "identity") +
      scale_fill_manual(values = colors)+
      theme_bw() +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()
      )
    ggsave(p2)
  }else{ # sort the data using bulk CNV score
    ArmGenSumm <- CtArmGenes(cdt = c7, data = data2)
    df1 <- ArmGenSumm %>% filter(arms %in% arm) # sort the spots based on the CNV signal found in the bulk data
    genes <- df1[, "genes"]
    rs <- df1[, "gene_row"]
    df <- cdt_filt_sort(cdt = c7, genes = genes, gainLoss = gainLoss, rs = rs)
    write.table(df, file = paste0(getwd(), "/CNVs.cdt"), sep = "\t", row.names = FALSE, quote = FALSE)
  }
}


