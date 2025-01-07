#'
#' Simple script to load data in DEXSeq and perform
#' Differential Exon Usage analysis
#'

# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
#
# BiocManager::install("DEXSeq")

# BiocManager::install("GenomicRanges")
# BiocManager::install("ROCR")

library(DEXSeq)
library(GenomicRanges)
library(rtracklayer)
library(ROCR)
library(data.table)
load("./differential_exons.RData")

countFiles = list.files("dexseq", full=T)
names(countFiles) <- gsub("_dexseq.txt", "", basename(countFiles))

## prepare the sample annotation
group <- rep(1,length(countFiles))
# sample names ending in 6..10 are in group 2
group[grepl("[06-9]$", names(countFiles))] <- 2
sampleTable = data.frame(condition=factor(group))
rownames(sampleTable) = names(countFiles)

## load the data
dxd = DEXSeqDataSetFromHTSeq(
  countFiles,
  sampleData=sampleTable,
  design= ~ sample + exon + condition:exon,
  flattenedfile = "./../../BamFiles/annotation_dexseq_b37.gff" )

## estimate the size factors
dxd = estimateSizeFactors( dxd ) 
## fit function to compute the variance in dependence on the mean
dxd = estimateDispersions( dxd ) 
## check the fit of the variance (dispersion) model
plotDispEsts( dxd ) 
## run the test for differential exon usage
dxd = testForDEU( dxd ) 
## obtain the results
dxr1 = DEXSeqResults( dxd )

res  <- as.data.frame(dxr1)

# flatten df (due to list)
res_flat <- as.data.frame(lapply(res, function(col) {
  if (is.list(col)) {
    sapply(col, toString)  # Convert list elements to strings
  } else {
    col
  }
}))


# get gtf
gtf <- rtracklayer::import("../../BamFiles/annotation_b37.gtf")
exons_gtf <- gtf[gtf$type == "exon"]
true_exons  <- differential.skipped
psi_results  <- fread("../output.txt")

# First, create GRanges object from DEXSeq results
dexseq_ranges <- GRanges(
  seqnames = res$genomicData.seqnames,
  ranges = IRanges(
    start = res$genomicData.start,
    end = res$genomicData.end
  ),
  strand = res$genomicData.strand,
  gene_id = res$groupID
)

head(dexseq_ranges)
# GRanges object with 6 ranges and 1 metadata column:
#       seqnames            ranges strand |           gene_id
#          <Rle>         <IRanges>  <Rle> |       <character>
#   [1]        X 99839799-99840063      + | ENSG00000000005.5
#   [2]        X 99840228-99840359      + | ENSG00000000005.5
#   [3]        X 99848621-99848891      + | ENSG00000000005.5
#   [4]        X 99848892-99849032      + | ENSG00000000005.5
#   [5]        X 99849258-99849359      + | ENSG00000000005.5
#   [6]        X 99852501-99852528      + | ENSG00000000005.5
#   -------
#   seqinfo: 24 sequences from an unspecified genome; no seqlengths

psi_results[, c("start", "end") := tstrsplit(exon, "-", type.convert = TRUE)]
psi_results[, `:=`(
  start = as.numeric(start),
  end = as.numeric(end)
)]


# Match gene IDs in psi_results with those in exons_gtf and extract chromosome information
exons_gtf_seqnames <- as.character(seqnames(exons_gtf))
psi_results$chromosome <- exons_gtf_seqnames[match(psi_results$gene, exons_gtf$gene_id)]

# Now create the GRanges object with the extracted chromosome info
psi_ranges <- GRanges(
  seqnames = psi_results$chromosome,  # Use the matched chromosome info
  ranges = IRanges(start = psi_results$start, end = psi_results$end),
  gene_id = psi_results$gene
)

# Check the resulting GRanges object
head(psi_ranges)
# GRanges object with 6 ranges and 1 metadata column:
#       seqnames            ranges strand |           gene_id
#          <Rle>         <IRanges>  <Rle> |       <character>
#   [1]       10 13275735-13275801      * | ENSG00000165623.5
#   [2]       16 70346512-70346561      * | ENSG00000260537.1
#   [3]       16 70348805-70348859      * | ENSG00000260537.1
#   [4]        7 16503941-16504013      * | ENSG00000171243.7
#   [5]       15 30700391-30700499      * | ENSG00000186399.8
#   [6]       15 30700605-30700692      * | ENSG00000186399.8
#   -------
#   seqinfo: 24 sequences from an unspecified genome; no seqlengths

# Find overlaps between DEXSeq results and annotation

dexseq_overlaps <- findOverlaps(dexseq_ranges, exons_gtf)

dexseq_matched <- data.frame(
  dexseq_id = dexseq_ranges$gene_id[queryHits(dexseq_overlaps)],
  gtf_gene_id = exons_gtf$gene_id[subjectHits(dexseq_overlaps)],
  gtf_exon_start = start(exons_gtf)[subjectHits(dexseq_overlaps)],
  gtf_exon_end = end(exons_gtf)[subjectHits(dexseq_overlaps)]
)
head(dexseq_matched)
#           dexseq_id       gtf_gene_id gtf_exon_start gtf_exon_end
# 1 ENSG00000000005.5 ENSG00000000005.5       99839799     99840063
# 2 ENSG00000000005.5 ENSG00000000005.5       99840228     99840359
# 3 ENSG00000000005.5 ENSG00000000005.5       99848621     99849032
# 4 ENSG00000000005.5 ENSG00000000005.5       99848892     99849032
# 5 ENSG00000000005.5 ENSG00000000005.5       99848621     99849032
# 6 ENSG00000000005.5 ENSG00000000005.5       99849258     99849359
# Find overlaps between PSI results and annotation

psi_overlaps <- findOverlaps(psi_ranges, exons_gtf)
psi_matched <- data.frame(
  psi_id = psi_ranges$gene_id[queryHits(psi_overlaps)],
  gtf_gene_id = exons_gtf$gene_id[subjectHits(psi_overlaps)],
  gtf_exon_start = start(exons_gtf)[subjectHits(psi_overlaps)],
  gtf_exon_end = end(exons_gtf)[subjectHits(psi_overlaps)]
)

head(psi_matched)
#              psi_id       gtf_gene_id gtf_exon_start gtf_exon_end
# 1 ENSG00000165623.5 ENSG00000165623.5       13275735     13275800
# 2 ENSG00000260537.1 ENSG00000260537.1       70346512     70346560
# 3 ENSG00000260537.1 ENSG00000260537.1       70348805     70348858
# 4 ENSG00000171243.7 ENSG00000171243.7       16503941     16504012
# 5 ENSG00000186399.8 ENSG00000186399.8       30700391     30700498
# 6 ENSG00000186399.8 ENSG00000186399.8       30700605     30700691


# update coords
res$genomicData.start[queryHits(dexseq_overlaps)] <- dexseq_matched$gtf_exon_start
res$genomicData.end[queryHits(dexseq_overlaps)] <- dexseq_matched$gtf_exon_end

head(res)

psi_results$start[queryHits(psi_overlaps)] <- psi_matched$gtf_exon_start
psi_results$end[queryHits(psi_overlaps)] <- psi_matched$gtf_exon_end

# Compare significant results
# For DEXSeq, consider results with padj < 0.05 as significant
dexseq_sig <- res[!is.na(res$padj) & res$padj < 0.05, ]
head(dexseq_sig)
dexseq_sig_ranges <- GRanges(
  seqnames = dexseq_sig$genomicData.seqnames,
  ranges = IRanges(
    start = dexseq_sig$genomicData.start,
    end = dexseq_sig$genomicData.end
  )
)

# For PSI results, consider results with padj < 0.05 as significant
psi_sig <- psi_results[psi_results$padj < 0.05, ]
psi_sig
psi_sig_ranges <- GRanges(
  seqnames = psi_sig$chromosome,
  ranges = IRanges(
    start = psi_sig$start,
    end = psi_sig$end
  )
)


# Calculate performance with correct ground truth
calculate_performance_corrected <- function(pred_ranges, true_ranges) {
  overlaps <- findOverlaps(pred_ranges, true_ranges)
  
  # Calculate metrics
  true_positives <- length(unique(queryHits(overlaps)))
  false_positives <- length(pred_ranges) - true_positives
  false_negatives <- length(true_ranges) - true_positives
  
  # Calculate precision, recall, and F1 score
  precision <- true_positives / (true_positives + false_positives)
  recall <- true_positives / (true_positives + false_negatives)
  f1_score <- 2 * (precision * recall) / (precision + recall)
  
  return(list(
    precision = precision,
    recall = recall,
    f1_score = f1_score,
    true_positives = true_positives,
    false_positives = false_positives,
    false_negatives = false_negatives,
    total_predictions = length(pred_ranges),
    total_true = length(true_ranges)
  ))
}

# Print dataset sizes for verification
print(paste("Number of true differential exons:", length(true_exons)))
print(paste("Number of DEXSeq significant exons:", length(dexseq_sig_ranges)))
print(paste("Number of PSI significant exons:", length(psi_sig_ranges)))

# Calculate performance metrics using true differential exons
dexseq_performance <- calculate_performance_corrected(dexseq_sig_ranges, true_exons)
psi_performance <- calculate_performance_corrected(psi_sig_ranges, true_exons)

# Create detailed summary table
performance_summary <- data.frame(
  Metric = c("True Positives", "False Positives", "False Negatives", 
             "Precision", "Recall", "F1 Score",
             "Total Predictions", "Total True Positives"),
  DEXSeq = c(dexseq_performance$true_positives,
             dexseq_performance$false_positives,
             dexseq_performance$false_negatives,
             dexseq_performance$precision,
             dexseq_performance$recall,
             dexseq_performance$f1_score,
             dexseq_performance$total_predictions,
             dexseq_performance$total_true),
  PSI = c(psi_performance$true_positives,
          psi_performance$false_positives,
          psi_performance$false_negatives,
          psi_performance$precision,
          psi_performance$recall,
          psi_performance$f1_score,
          psi_performance$total_predictions,
          psi_performance$total_true)
)

# Print summary table with formatted numbers
print("Performance Summary:")
print(performance_summary)

# Create overlap analysis
dexseq_true_overlaps <- findOverlaps(dexseq_sig_ranges, true_exons)
psi_true_overlaps <- findOverlaps(psi_sig_ranges, true_exons)

# Print overlap details
print("\nDetailed Overlap Analysis:")
print(paste("DEXSeq overlaps with true exons:", length(unique(queryHits(dexseq_true_overlaps)))))
print(paste("PSI overlaps with true exons:", length(unique(queryHits(psi_true_overlaps)))))

# Let's also see which true differential exons were detected by both methods
dexseq_found <- subjectHits(dexseq_true_overlaps)
psi_found <- subjectHits(psi_true_overlaps)
common_found <- intersect(dexseq_found, psi_found)

print(paste("\nExons found by both methods:", length(common_found)))
print(paste("Exons found only by DEXSeq:", length(setdiff(dexseq_found, psi_found))))
print(paste("Exons found only by PSI:", length(setdiff(psi_found, dexseq_found))))


# First, we need to create prediction objects for ROCR
# For each exon in our analysis, we need:
# 1. The prediction scores (1 - p-value, so higher values = more likely to be differential)
# 2. The true labels (1 if in true_exons, 0 if not)

prepare_rocr_data <- function(ranges, pvals, true_ranges) {
  # Remove NA values
  valid_idx <- !is.na(pvals)
  ranges_clean <- ranges[valid_idx]
  pvals_clean <- pvals[valid_idx]
  
  # Find all overlaps with true ranges
  overlaps <- findOverlaps(ranges_clean, true_ranges)
  
  # Create labels vector (0 by default)
  labels <- rep(0, length(ranges_clean))
  # Set 1 for those that overlap with true ranges
  labels[queryHits(overlaps)] <- 1
  
  # Transform p-values to prediction scores (1 - pvalue)
  predictions <- 1 - pvals_clean
  
  return(list(predictions = predictions, labels = labels))
}

# Prepare data for DEXSeq
dexseq_data <- prepare_rocr_data(
  dexseq_ranges,
  res$pvalue,
  true_exons
)

# Prepare data for PSI
psi_data <- prepare_rocr_data(
  psi_ranges,
  psi_results$pvalue,
  true_exons
)

head(psi_data)
head(dexseq_data)

# Create ROCR prediction objects
dexseq_pred <- prediction(dexseq_data$predictions, dexseq_data$labels)
psi_pred <- prediction(psi_data$predictions, psi_data$labels)

# Calculate performance metrics
dexseq_perf <- performance(dexseq_pred, "tpr", "fpr")
psi_perf <- performance(psi_pred, "tpr", "fpr")

# Calculate AUC
dexseq_auc <- performance(dexseq_pred, "auc")@y.values[[1]]
psi_auc <- performance(psi_pred, "auc")@y.values[[1]]

# Create plot
pdf("roc_comparison.pdf", width = 7, height = 7)  # Open PDF device

# Set up the plot
plot(dexseq_perf, 
     col = "blue", 
     lwd = 2,
     main = "ROC Curve Comparison",
     xlab = "False Positive Rate (1 - Specificity)",
     ylab = "True Positive Rate (Sensitivity)")

# Add PSI curve
plot(psi_perf, 
     col = "red", 
     lwd = 2,
     add = TRUE)

# Add diagonal reference line
abline(0, 1, lty = 2, col = "gray")

# Add legend with AUC values
legend("bottomright",
       legend = c(
         paste("DEXSeq (AUC =", round(dexseq_auc, 3), ")"),
         paste("PSI (AUC =", round(psi_auc, 3), ")")
       ),
       col = c("blue", "red"),
       lwd = 2)

dev.off()  # Close PDF device

calculate_metrics_at_threshold <- function(pred_ranges, true_ranges, adj_pvals, threshold = 0.05) {
  
  # Remove NA values from adj_pvals, pred_ranges, and true_ranges
  valid_idx <- !is.na(adj_pvals)
  
  # Filter adj_pvals, pred_ranges, and true_ranges based on valid indices
  adj_pvals_clean <- adj_pvals[valid_idx]
  pred_ranges_clean <- pred_ranges[valid_idx]
  
  # Get significant results based on the threshold
  sig_idx <- adj_pvals_clean < threshold
  sig_ranges <- pred_ranges_clean[sig_idx]
  
  # Find overlaps between significant predicted ranges and true ranges
  overlaps <- findOverlaps(sig_ranges, true_ranges)
  
  # Calculate true positives, false positives, and false negatives
  true_positives <- length(unique(queryHits(overlaps)))  # Unique queryHits are true positives
  false_positives <- length(sig_ranges) - true_positives  # Remaining significant predicted ranges are false positives
  false_negatives <- length(true_ranges) - true_positives  # Remaining true ranges are false negatives
  
  # Calculate FDR (False Discovery Rate) and Sensitivity
  fdr <- false_positives / (false_positives + true_positives)
  sensitivity <- true_positives / (true_positives + false_negatives)
  
  return(list(
    FDR = fdr,
    Sensitivity = sensitivity
  ))
}

# Calculate metrics for both methods
dexseq_metrics <- calculate_metrics_at_threshold(dexseq_ranges, true_exons, res$padj)
psi_metrics <- calculate_metrics_at_threshold(psi_ranges, true_exons, psi_results$padj)

# Create summary table
metrics_summary <- data.frame(
  Method = c("DEXSeq", "PSI"),
  FDR = c(dexseq_metrics$FDR, psi_metrics$FDR),
  Sensitivity = c(dexseq_metrics$Sensitivity, psi_metrics$Sensitivity)
)

# Print results
print("Performance metrics at 5% adjusted p-value threshold:")
print(metrics_summary)

print("\nAUC values:")
print(paste("DEXSeq AUC:", round(dexseq_auc, 3)))
print(paste("PSI AUC:", round(psi_auc, 3)))
