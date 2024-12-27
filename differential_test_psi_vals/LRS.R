library("data.table")
library("dplyr")

LRS <- function(incl, total, group) {
    # step 1: calc MLE calses
    p0 <- sum(incl) / sum(total)
    p1 <- sum(incl[group == 1]) / sum(total[group == 1])
    p2 <- sum(incl[group == 2]) / sum(total[group == 2])
    
    # step 2: get log likelihood
    log_likelihood_reduced <- sum(dbinom(incl, total, p0, log = TRUE))
    log_likelihood_full <- sum(dbinom(incl[group == 1], total[group == 1], p1, log = TRUE)) + sum(dbinom(incl[group == 2], total[group == 2], p2, log = TRUE))
    
    # step 3: get LRS
    LRS_value <- -2 * (log_likelihood_reduced - log_likelihood_full)

    # chi test
    df <- 1 # since 2 - 1
    p_value <- pchisq(LRS_value, df = df, lower.tail = FALSE)
    
    return(list(
       p0 = p0,
       p1 = p1,
       p2 = p2,
       llreduced = log_likelihood_reduced,
       llfull = log_likelihood_full,
       lrs = LRS_value,
       pvalue = p_value
     ))
}

diff.splicing <- function(psi.files, group) {
  all_data <- lapply(seq_along(psi.files), function(i) {
    file <- psi.files[i]
    data <- read.table(file, header = TRUE, sep = "\t")
    data$group <- group[i]
    return(data)
  })
  
  combined_data <- bind_rows(all_data)

  # only look at unique gene-exon combinations
  unique_features <- unique(combined_data[c("gene", "exon")])
  
  results <- apply(unique_features, 1, function(feature) {
    subset_data <- combined_data[combined_data$gene == feature["gene"] &
                                combined_data$exon == feature["exon"], ]
    
    lrs_result <- LRS(subset_data$num_incl_reads,
                     subset_data$num_total_reads,
                     subset_data$group)
    
    return(c(
      gene = feature["gene"],
      exon = feature["exon"],
      p0 = lrs_result$p0,
      p1 = lrs_result$p1,
      p2 = lrs_result$p2,
      llreduced = lrs_result$llreduced,
      llfull = lrs_result$llfull,
      lrs = lrs_result$lrs,
      pvalue = lrs_result$pvalue
    ))
  })
  
  results_df <- as.data.frame(t(results))
  colnames(results_df)[1:2] <- c("gene", "exon")
  results_df$pvalue <- as.numeric(as.character(results_df$pvalue))
  results_df$padj <- p.adjust(results_df$pvalue, method = "BH")
  
  return(results_df)
}

# Testing
# load("./diff_psi_test.RData")
# psi.files = c("./../psi/sample1.psi", "./../psi/sample2.psi", "./../psi/sample3.psi", "./../psi/sample4.psi", "./../psi/sample5.psi", "./../psi/sample6.psi", "./../psi/sample7.psi", "./../psi/sample8.psi", "./../psi/sample9.psi", "./../psi/sample10.psi")
# res <- diff.splicing(psi.files, group)
# head(res)
# library(data.table)
# res  <- data.table(res)
# res[res[, gene == "ENSG00000001167.10"]]
# write.table(res, file="output.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
# res
