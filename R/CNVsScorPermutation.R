#' Generate CNVs_Scores of 100 permutations.
#'
#' @param m number of times of permutation, by default is 100
#' @param cdt0 cdt file used to calculate CNVs
#' @param armGens dataframe of the output of CtArmGenes function
#' @param arm A vector of arm with gene gain and loss
#' @param gainLoss A vector of data showing gain(+1) and loss(-1)
#'
#' @return A dataframe of m times permutation values.
#' @export
#'
#' @examples \dontrun{
#' print("check user guider")
#' }
CNVpermut <- function(m = 100, cdt0, armGens, arm, gainLoss) {
  if (!is.data.frame(cdt0)) {
    stop("cdt0 must be a dataframe")
  }

  if (!is.data.frame(armGens)) {
    stop("armGens must be a data frame")
  }

  if (length(arm) != length(gainLoss)) {
    stop("arm and gainLoss should have the same length!")
  }

  arms <- NULL
  # bulk arm gain or loss info
  df1 <- armGens %>% filter(arms %in% arm)
  genes <- df1[, "genes"]
  rs <- df1[, "gene_row"]

  # spot cdt file
  rows2 <- cdt0[, c(1, 2)]
  cdt11 <- cdt0[, -c(1, 2)]
  cdt <- cdt11

  # create combined permutation file.
  permuts <- data.frame(matrix(nrow = dim(cdt)[2], ncol = 0))
  for (lif in 1:m) {
    df1 <- data.frame(t(apply(cdt, 1, sample))) # permute colums within each row
    df1 <- data.frame(apply(df1, 2, sample))
    df2 <- cbind(rows2, df1)
    # write.table(df2, file = paste0("/Volumes/Bastian/Limin/Ji_data/ST_paper/mergedP6/permutation/p6_permuted",lif,".cdt"), sep = "\t",row.names = FALSE, quote = FALSE)

    # create permuated scores
    cdt <- df2
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
    permuts <- cbind(permuts, newR)
  }
  colnames(permuts) <- c(paste0("permut", seq(m)))
  return(permuts)
}


#' Generate your sample CNV score dataframe.
#'
#' @param cdt0 sample cdt file
#' @param armGens dataframe of the output of CtArmGenes function
#' @param arm A vector of arm with gene gain and loss
#' @param gainLoss A vector of data showing gain(+1) and loss(-1)
#'
#' @return A dataframe of your sample CNV score
#' @export
#'
#' @examples \dontrun{
#' print("check user guider")
#' }
caseCNVscore <- function(cdt0, armGens, arm, gainLoss) {
  if (!is.data.frame(cdt0)) {
    stop("cdt0 must be a dataframe")
  }

  if (!is.data.frame(armGens)) {
    stop("armGens must be a data frame")
  }

  if (length(arm) != length(gainLoss)) {
    stop("arm and gainLoss should have the same length!")
  }


  arms <- NULL
  df1 <- armGens %>% filter(arms %in% arm)
  genes <- df1[, "genes"]
  rs <- df1[, "gene_row"]

  name <- cdt0[, c(1, 2)]
  cdt <- cdt0[, -c(1, 2)]
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

  df <- data.frame(colnames(cdt), newR)
  df <- df[order(df$newR, decreasing = TRUE), ]
  df$Rank <- seq(nrow(df))
  colnames(df) <- c("Spot", "CNVScore", "Rank")
  df <- df[, c("Rank", "Spot", "CNVScore")]
  return(df)
}


#' Generate patient scores and permutation scores and rank them from high to low.
#'
#' @param sptScor This is the output of caseCNVscore function.
#' @param permuts This is the output of CNVpermut function.
#'
#' @return A list of 2 dataframes. The first is the sample_score_summary,
#' the second is the permuts score summary.
#' @export
#'
#' @examples \dontrun{
#' print("check user guider")
#' }
permutScore <- function(sptScor, permuts) {
  if (!is.data.frame(sptScor)) {
    stop("sptScor must be a dataframe.")
  }

  if (!is.data.frame(permuts)) {
    stop("permuts must be a dataframe.")
  }


  CNVscor <- NULL
  FDR <- NULL
  Fraction_at_that_scor_or_lower_obs <- NULL
  Fraction_at_that_scor_or_lower_permut <- NULL
  Fraction_at_this_diff_or_higher_observed <- NULL
  Fraction_at_this_diff_or_higher_permut <- NULL
  Fraction_at_this_scor_or_higher_obs <- NULL
  Fraction_at_this_scor_or_higher_permut <- NULL
  Rank <- NULL
  Rank_reverse <- NULL
  arms <- NULL
  diff_mean <- NULL
  meanDiffRank <- NULL
  spot <- NULL


  # case spotScore
  colnames(sptScor) <- c("Rank", "spot", "CNVscor")
  sptScor$Fraction_at_this_scor_or_higher_obs <- sptScor$Rank / dim(sptScor)[1] # add Fraction_at_this_scor_or_higher_obs
  sptScor <- sptScor[order(sptScor$CNVscor, decreasing = FALSE), ]
  sptScor$Rank_reverse <- seq(dim(sptScor)[1])
  sptScor$Fraction_at_that_scor_or_lower_obs <- sptScor$Rank_reverse / dim(sptScor)[1] # add Fraction_at_that_scor_or_lower_obs
  sptScor <- sptScor[order(sptScor$CNVscor, decreasing = TRUE), ]
  sptScor$diff_mean <- abs(sptScor$CNVscor - mean(sptScor$CNVscor)) # absolute value of CNVscor mean difference
  sptScor <- sptScor[order(sptScor$diff_mean, decreasing = TRUE), ]
  sptScor$meanDiffRank <- seq(dim(sptScor)[1])
  sptScor$Fraction_at_this_diff_or_higher_observed <- sptScor$meanDiffRank / dim(sptScor)[1]
  sptScor <- sptScor[order(sptScor$Rank), ]

  # permuts scor
  p1 <- c()
  for (i in 1:dim(permuts)[2]) {
    p1 <- append(p1, permuts[, i])
  }
  r <- length(p1)
  permD <- data.frame(seq(r), sort(p1, decreasing = TRUE))
  colnames(permD) <- c("Rank", "PermutScore")

  nr <- dim(permD)[1]
  permD$Fraction_at_this_scor_or_higher_permut <- permD$Rank / nr
  permD <- permD[order(permD$PermutScore), ]
  permD$Rank_reverse <- seq(nr)
  permD$Fraction_at_that_scor_or_lower_permut <- permD$Rank_reverse / nr
  permD$diff_Mean <- abs(permD$PermutScore - mean(permD$PermutScore))
  permD <- permD[order(permD$diff_Mean, decreasing = TRUE), ]
  permD$meanDiffRank <- seq(nr)
  permD$Fraction_at_this_diff_or_higher_permut <- permD$meanDiffRank / nr
  permD <- permD[order(permD$Rank), ]

  # summary table
  sptScor$Fraction_at_this_scor_or_higher_permut <- rep(0, dim(sptScor)[1])
  sptScor$Fraction_at_that_scor_or_lower_permut <- rep(0, dim(sptScor)[1])
  sptScor$Fraction_at_this_diff_or_higher_permut <- rep(0, dim(sptScor)[1])

  numb <- c()
  for (v in 1:dim(sptScor)[1]) {
    for (i in 1:nr) {
      if (permD[i, "PermutScore"] < sptScor[v, "CNVscor"]) {
        numb <- append(numb, i)
        val1 <- permD[i, "Fraction_at_this_scor_or_higher_permut"]
        val2 <- permD[i, "Fraction_at_that_scor_or_lower_permut"]
        val3 <- permD[i, "Fraction_at_this_diff_or_higher_permut"]
        sptScor[v, "Fraction_at_this_scor_or_higher_permut"] <- val1
        sptScor[v, "Fraction_at_that_scor_or_lower_permut"] <- val2
        sptScor[v, "Fraction_at_this_diff_or_higher_permut"] <- val3
        break
      }
    }
  }

  # Add FDR column
  sptScor$FDR <- sptScor$Fraction_at_this_scor_or_higher_permut / sptScor$Fraction_at_this_scor_or_higher_obs
  sptScor1 <- dplyr::select(sptScor, Rank, Rank_reverse, spot, CNVscor, Fraction_at_this_scor_or_higher_obs, Fraction_at_that_scor_or_lower_obs, Fraction_at_this_diff_or_higher_observed, diff_mean, meanDiffRank, Fraction_at_this_scor_or_higher_permut, Fraction_at_that_scor_or_lower_permut, Fraction_at_this_diff_or_higher_permut, FDR)

  returnList <- list(sptScor1, permD)
  return(returnList)
}



#' Generate tumor groups based on FDR rate for Loupe visualization.
#'
#' @param a CNVscor cutoff less than a for Non_tumor
#' @param b CNVscor cutoff for potential_tumor and Tossup ( b=< Potential_tumor_CNVscor < c;
#' a =< Tossup_CNVscor_cutoff < b)
#' @param c CNVscor cutoff greater than c for Tumor spots
#' @param data1 Sample/case CNVscore summary, the output of permutScore function
#' @param data2 Sample grouped spots representative summary
#' @param data3 Grouped_spots summary, which spots into the same group, which
#' should be given the same color for loupe visualization
#' @param data4 All spots summary
#'
#' @return A dataframe of FDR-based tumor groups, used for Loupe visualization
#' @export
#'
#' @examples \dontrun{
#' print("check user guider")
#' }
fdrTumor <- function(a = 0, b = 1.82, c = 4.4, data1, data2, data3, data4) {
  if (!is.numeric(a)) {
    stop("a must be numeric!")
  }

  if (!is.numeric(b)) {
    stop("b must be numeric! ")
  }

  if (!is.numeric(c)) {
    stop("c must be numeric")
  }

  stopifnot(is.data.frame(data1), is.data.frame(data2), is.data.frame(data3), is.data.frame(data4))

  spot <- NULL
  final <- NULL
  barcode <- NULL
  Group <- NULL
  n <- nrow(data1)
  data1$Group <- rep(0, n)
  for (i in 1:n) {
    if (data1[i, "CNVscor"] < a) {
      data1[i, "Group"] <- "Non_tumor"
    } else if (data1[i, "CNVscor"] > c) {
      data1[i, "Group"] <- "Tumor"
    } else if (data1[i, "CNVscor"] < c & data1[i, "CNVscor"] >= b) {
      data1[i, "Group"] <- "Potential_Tumor"
    } else if (data1[i, "CNVscor"] < b & data1[i, "CNVscor"] >= a) {
      data1[i, "Group"] <- "Tossup"
    }
  }

  # visualization
  data1 <- data1[, c("spot", "Group")]
  ## patient6 rep2
  r2 <- data1 %>% filter(spot %in% data2$spot)
  r22 <- data2 %>% filter(spot %in% r2$spot)
  r22$barcode <- str_replace(r22$barcode, ".1", "-1")
  colnames(r2) <- c("spot", "Group")
  p6rep2 <- merge(r22, r2, by = "spot")
  p6rep2 <- p6rep2[, c("barcode", "Group")]

  # including spots that are being represented by these spots
  reprt2 <- data3 %>% filter(barcode %in% p6rep2$barcode)
  sub <- p6rep2 %>% filter(barcode %in% reprt2$barcode)
  reprt2 <- merge(reprt2, sub, by = "barcode")
  data3$Group <- rep(0, dim(data3)[1])

  # assign group name to spots
  for (i in 1:dim(data3)[1]) {
    for (j in 1:dim(reprt2)[1]) {
      if (data3[i, "cluster"] == reprt2[j, "cluster"] & data3[i, "groups"] == reprt2[j, "groups"]) {
        data3[i, "Group"] <- reprt2[j, "Group"]
      }
    }
  }

  d2 <- data3 %>%
    filter(!Group == 0) %>%
    dplyr::select(barcode, Group)
  d3 <- rbind(p6rep2, d2)

  # read in all barcodes
  d41 <- data4 %>% filter(!barcode %in% d3$barcode)
  d5 <- data.frame(d41$barcode, rep("out", dim(d41)[1]))
  colnames(d5) <- c("barcode", "Group")
  d6 <- distinct(rbind(d3, d5))
  return(d6)
}
