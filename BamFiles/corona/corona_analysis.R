library("dplyr")
library("data.table")
library("DESeq2")
library("ggplot2")
library("pheatmap")


dt_star  <- fread("./gene.counts.star")
counts  <- as.matrix(dt_star, , rownames = "Geneid")
colnames(counts) <- sub(".bam$", "", colnames(counts))
meta_data <- fread("./sample.list")

meta_data <- meta_data %>%
  mutate(
    time = sub(".*_(\\d+h)$", "\\1", condition),
    condition = sub("_(\\d+h)$", "", condition)
  )
meta_data <- meta_data %>%
  mutate(
    condition = factor(condition),
    time = factor(time)
)

with(meta_data,
       table(condition, time))



dds_m = DESeqDataSetFromMatrix(
  countData = counts,
  colData   = meta_data,
  design    = ~ condition)

dds  <- DESeq(dds_m)

res = results(dds)
res[order(res$padj), ] |> head()


ggplot(as(res, "data.frame"), aes(x = pvalue)) +
  geom_histogram(binwidth = 0.01, fill = "Royalblue", boundary = 0)


plotMA(dds, ylim = c( -2, 2))


pas_rlog = rlogTransformation(dds)
plotPCA(pas_rlog, intgroup=c("condition", "time")) + coord_fixed()

select = order(rowMeans(assay(pas_rlog)), decreasing = TRUE)[1:30]
pheatmap( assay(pas_rlog)[select, ],
     scale = "row",
     annotation_col = as.data.frame(
        colData(pas_rlog)[, c("condition", "time")] ))


ddsTwoFactor = dds_m
design(ddsTwoFactor) = formula(~ time + condition)
ddsTwoFactor = DESeq(ddsTwoFactor)

res2 = results(ddsTwoFactor)
head(res2, n = 3)

trsf = function(x) ifelse(is.na(x), 0, (-log10(x)) ^ (1/6))
ggplot(tibble(pOne = res$pvalue,
              pTwo = res2$pvalue),
    aes(x = trsf(pOne), y = trsf(pTwo))) +
    geom_hex(bins = 75) + coord_fixed() +
    xlab("Single factor analysis (condition)") +
    ylab("Two factor analysis (type + condition)") +
    geom_abline(col = "orange")
