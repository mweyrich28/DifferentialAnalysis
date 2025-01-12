library(data.table)
library(patchwork)
library(ggplot2)
psi_values  <- fread("./output_chr.tsv")
#                    gene                exon        p0        p1        p2 llreduced    llfull         lrs      pvalue    chr       padj
#                  <char>              <char>     <num>     <num>     <num>     <num>     <num>       <num>       <num> <char>      <num>
#   1:  ENSG00000111701.6     7807201-7807229 0.6141732 0.5555556 0.6575342 -19.29236 -18.61300  1.35871339 0.243760586     12 0.49152973
#   2:  ENSG00000102802.5   31495698-31495984 0.7100671 0.6919060 0.7292818 -25.80034 -25.16816  1.26437271 0.260825016     13 0.50849919
#   3:  ENSG00000197579.3   32550772-32550967 0.7792969 0.7511521 0.8000000 -23.13918 -22.27742  1.72351950 0.189240638      9 0.41031615
#   4: ENSG00000144460.10 226378087-226378389 0.6558237 0.6225000 0.6799277 -33.28296 -31.59294  3.38002404 0.065991092      2 0.18011686
#   5: ENSG00000144460.10 226446657-226447752 0.8376238 0.8193780 0.8505068 -33.74650 -32.01408  3.46484053 0.062686552      2 0.17313429
#  ---                                                                                                                                   
# 228:  ENSG00000177459.6   99101304-99102258 0.8553001 0.8380567 0.8675624 -30.78575 -29.27479  3.02192264 0.082146025      8 0.21413346
# 229:  ENSG00000148835.9 105145083-105145248 0.6554455 0.6596639 0.6516854 -28.68912 -28.67138  0.03547648 0.850600720     10 0.95333028
# 230:  ENSG00000142634.8   15754394-15754506 0.7123519 0.6807692 0.7371601 -23.13804 -22.01201  2.25205820 0.133436819      1 0.30957342
# 231:  ENSG00000250741.2   18745111-18745386 0.7815951 0.7282913 0.8231441 -30.30504 -25.05322 10.50363099 0.001191402      2 0.01454764
# 232:  ENSG00000250741.2   18768260-18768440 0.7453988 0.7075812 0.7733333 -25.62633 -23.82307  3.60652109 0.057553396      2 0.16681090



# load("./dexseq/differential_exons.RData")
# true_exons  <- differential.skipped

q <- ggplot(psi_values, aes(p0)) +
  geom_density(aes(fill = "p0"), color = "black", alpha = 0.8) + 
  scale_fill_manual(values = c("p0" = "blue"), name = "Typ", labels = c("p0")) +
  labs(
    title = "p0 Verteilung",
    x = "Wahrscheinlichkeit",
    y = "Häufigkeit"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 16, face = "bold"), 
    legend.text = element_text(size = 14) 
  )

psi_values_melted  <- melt(psi_values, 
                           measure.vars = c("p1", "p2"),
                           variable.name = "Typ",
                           value.name = "prob"
)
psi_values_melted

p <- ggplot(psi_values_melted, aes(x = prob, fill = Typ)) + 
  geom_density( color = "black", alpha = 0.8) + 
  labs(
    title = "p1, p2 Verteilung",
    x = "Wahrscheinlichkeit",
    y = "Häufigkeit"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 16, face = "bold"), 
    legend.text = element_text(size = 14) 
  )
p / q
ggsave("../das_report/plots/prob.png", plot = , width = 14, height = 8, dpi = 300)


chromosome_levels <- c(as.character(1:22), "X", "Y", "MT")

ggplot(psi_values, aes(x = chr, y = neg_log10_padj, fill = as.factor(as.integer(as.factor(chr)) %% 2))) +
  geom_jitter(width = 0.2, alpha = 0.7, size = 4, shape = 21, color = "black") + 
  geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed") + 
  scale_fill_manual(
    values = c("blue", "lightblue"), 
    guide = "none" 
  ) +
  labs(
    title = "Korrigierte P-Werte nach Chromosom",
    x = "Chromosom",
    y = expression(-log[10]~"(Padj)")
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 25), 
    axis.title.x = element_text(size = 21), 
    axis.title.y = element_text(size = 21), 
    axis.text.x = element_text(size = 18, angle = 45, hjust = 1), 
    axis.text.y = element_text(size = 18)
  )


p <- ggplot(psi_values, aes(x = chr, y = neg_log10_padj, fill = as.factor(as.integer(as.factor(chr)) %% 2))) +
  geom_jitter(width = 0.2, alpha = 0.8, size = 3, shape = 21, color = "black") + 
  geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed") + 
  scale_fill_manual(
    values = c("blue", "lightblue"), 
    guide = "none" 
  ) +
  labs(
    title = "Korrigierte P-Werte nach Chromosom",
    x = "Chromosom",
    y = expression(-log[10]~"(Padj)")
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 25), 
    axis.title.x = element_text(size = 21), 
    axis.title.y = element_text(size = 21), 
    axis.text.x = element_text(size = 18, angle = 45, hjust = 1), 
    axis.text.y = element_text(size = 18)
  )

# Save the plot as a PNG file
ggsave("../das_report/plots/neg_log10_padj_by_chromosome.png", plot = p, width = 14, height = 8, dpi = 300)



library(ggplot2)
data <- data.frame(
  Sample = c("bam_1", "bam_2", "bam_3", "bam_4", "bam_5", 
             "bam_6", "bam_7", "bam_8", "bam_9", "bam_10"),
  Time = c(1.193, 1.164, 1.167, 1.144, 1.183, 1.316, 1.307, 1.210, 1.223, 1.216)
)


# Create the bar plot
bar_plot <- ggplot(data, aes(x = Sample, y = Time)) +
  geom_bar(stat = "identity", fill = "blue", color = "black", alpha = 0.8) + 
  labs(
    title = "Durchschnittliche Laufzeit pro Datei",
    x = "Datei",
    y = "Zeit in Sekunden" 
  ) +
  scale_x_discrete(labels = data$Sample) + 
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 25), 
    axis.title.x = element_text(size = 21), 
    axis.title.y = element_text(size = 21), 
    axis.text.x = element_text(size = 18, angle = 45, hjust = 1), 
    axis.text.y = element_text(size = 18)

  )

print(bar_plot)

# Display the plot
ggsave("../das_report/plots/time.png", plot =bar_plot , width = 14, height = 8, dpi = 300)


data <- data.table(
    A =  c(150,385,334,487,1100,183,233,2031,399,286,227,315,210,155,2088,450,211,401,202,2330,132,1018,670,254,513,679,368,473,1025,80,398,511,293,240,1703,568,266,230,830,544,664,144,1145,222,213,396,142,8,250,605,380,899,346,396,142,256,340,350,181,2762,399,424,216,189,522,511,167,223,228,183,234,257,698,385,393,1278,292,154,786,204,156,278,539,179,1019,1073,432,150,1140,139,359,38,848,219,257,222,368,167,176,343,1041,171,227,294,240,343,716,1107,611,227,232,437,1134,196,240,137,387,250,865,168,602,871,175,560,175,1033,249,346,300,224,666,372,571,647,220,190,604,224,622,287,273,783,259,346,184,780,667,1052,143,253,308,326,184,829,795,305,566,781,430,223,260,524,198,228,185,238,680,198,188,431,247,160,372,261,179,769,269,743,231,755,626,257,321,289,123,167,286,321,212,126,2108,299,722),
    B = c( 106,273,179,355,507 ,87,139,644,179,157,136,181,112,106,749,187,147,267,128,670,89,332,273,141,204,272,154,295,316,41,221,282,125,151,485,238,176,152,421,291,431,83,342,115,122,178,83,6,163,281,168,298,174,144,88,138,241,207,121,734,148,201,153,125,215,198,98,145,113,104,129,143,241,187,166,375,160,92,450,123,95,174,210,117,507,324,246,84,435,104,191,27,399,164,155,98,170,102,110,158,323,116,118,160,141,199,273,517,211,156,139,144,336,127,150,85,257,188,341,122,224,325,94,210,101,316,154,198,158,129,292,151,199,253,136,127,304,129,236,163,126,347,136,178,131,423,244,313,83,146,146,186,108,231,389,181,331,254,256,145,137,249,131,169,114,128,245,122,116,297,154,98,191,140,137,319,163,348,154,358,219,160,221,136,84,113,183,137,142,90,862,151,279),
    Gene = c("ENSG00000165623.5", "ENSG00000260537.1", "ENSG00000171243.7", "ENSG00000186399.8", "ENSG00000111885.5", "ENSG00000197579.3", "ENSG00000108242.8", "ENSG00000155886.7", "ENSG00000164615.3", "ENSG00000171136.6", "ENSG00000176428.5", "ENSG00000257529.1", "ENSG00000140379.7", "ENSG00000248672.1", "ENSG00000151012.9", "ENSG00000137745.7", "ENSG00000154589.2", "ENSG00000268083.1", "ENSG00000170293.4", "ENSG00000108256.4", "ENSG00000163600.8", "ENSG00000187902.7", "ENSG00000174808.7", "ENSG00000160339.11", "ENSG00000102802.5", "ENSG00000124469.6", "ENSG00000177459.6", "ENSG00000250741.2", "ENSG00000081052.10", "ENSG00000198035.9", "ENSG00000100147.9", "ENSG00000087128.5", "ENSG00000196844.4", "ENSG00000249624.4", "ENSG00000162924.9", "ENSG00000142634.8", "ENSG00000172426.11", "ENSG00000108255.3", "ENSG00000109674.3", "ENSG00000144659.6", "ENSG00000244255.1", "ENSG00000186583.7", "ENSG00000151233.6", "ENSG00000107731.8", "ENSG00000204149.5", "ENSG00000188916.4", "ENSG00000155282.7", "ENSG00000198946.3", "ENSG00000042304.6", "ENSG00000148110.11", "ENSG00000178690.2", "ENSG00000173674.6", "ENSG00000099284.9", "ENSG00000127528.5", "ENSG00000170835.10", "ENSG00000147223.5", "ENSG00000196860.3", "ENSG00000196407.7", "ENSG00000102245.3", "ENSG00000165186.9", "ENSG00000078579.8", "ENSG00000120688.7", "ENSG00000142273.6", "ENSG00000088325.11", "ENSG00000139890.5", "ENSG00000138463.8", "ENSG00000101335.5", "ENSG00000170509.7", "ENSG00000172292.10", "ENSG00000101280.6", "ENSG00000204520.8", "ENSG00000130224.10", "ENSG00000102580.10", "ENSG00000174903.10", "ENSG00000101213.5", "ENSG00000130338.8", "ENSG00000158481.8", "ENSG00000106780.7", "ENSG00000204645.4", "ENSG00000064995.12", "ENSG00000148835.9", "ENSG00000100373.5", "ENSG00000013016.10", "ENSG00000182271.8", "ENSG00000115850.5", "ENSG00000158195.6", "ENSG00000147378.7", "ENSG00000146013.6", "ENSG00000144460.10", "ENSG00000152592.9", "ENSG00000163421.4", "ENSG00000111701.6", "ENSG00000256646.3", "ENSG00000188747.4", "ENSG00000166143.5", "ENSG00000117479.8", "ENSG00000124789.7", "ENSG00000021574.7", "ENSG00000053438.7", "ENSG00000229544.6", "ENSG00000173705.4", "ENSG00000143839.12", "ENSG00000148943.7", "ENSG00000077935.12", "ENSG00000116962.10", "ENSG00000188219.10", "ENSG00000001167.10", "ENSG00000163624.5", "ENSG00000183783.6", "ENSG00000137161.12", "ENSG00000168930.9", "ENSG00000244509.3", "ENSG00000163840.5", "ENSG00000135226.12", "ENSG00000057149.10", "ENSG00000198553.4", "ENSG00000079689.9", "ENSG00000230031.6", "ENSG00000105887.10", "ENSG00000117560.6", "ENSG00000101282.4", "ENSG00000165661.11", "ENSG00000189099.7", "ENSG00000150361.7", "ENSG00000152454.3", "ENSG00000197953.5", "ENSG00000197646.6", "ENSG00000165264.6", "ENSG00000112077.11", "ENSG00000125831.5", "ENSG00000198865.5", "ENSG00000185053.8", "ENSG00000147588.6", "ENSG00000112175.6", "ENSG00000240403.1", "ENSG00000081248.6", "ENSG00000187151.3", "ENSG00000124664.6", "ENSG00000146197.7", "ENSG00000105668.3", "ENSG00000170162.9", "ENSG00000056277.11", "ENSG00000268182.1", "ENSG00000135069.9", "ENSG00000196800.2", "ENSG00000123454.6", "ENSG00000213221.4", "ENSG00000196715.5", "ENSG00000165506.10", "ENSG00000125498.15", "ENSG00000188883.4", "ENSG00000104147.4", "ENSG00000163513.13", "ENSG00000128342.4", "ENSG00000198483.8", "ENSG00000100053.5", "ENSG00000169340.5", "ENSG00000163961.4", "ENSG00000162851.6", "ENSG00000148572.10", "ENSG00000259112.1", "ENSG00000132881.7", "ENSG00000242875.2", "ENSG00000116218.8", "ENSG00000128524.4", "ENSG00000101856.8", "ENSG00000177324.9", "ENSG00000102390.6", "ENSG00000101161.6", "ENSG00000080618.9", "ENSG00000259900.1", "ENSG00000158109.10", "ENSG00000270149.1", "ENSG00000125869.5", "ENSG00000133105.3", "ENSG00000086544.2", "ENSG00000104918.4", "ENSG00000214491.4", "ENSG00000167230.6", "ENSG00000214102.3", "ENSG00000123901.4", "ENSG00000101057.11", "ENSG00000272968.1", "ENSG00000169116.7", "ENSG00000251283.1", "ENSG00000183324.6", "ENSG00000265264.1", "ENSG00000204033.5", "ENSG00000265590.5", "ENSG00000159212.8", "ENSG00000100077.10", "ENSG00000114487.5", "ENSG00000149679.7")
)

mean(data$A)
mean(data$B)

data_m  <- melt(data,
                id.vars = "Gene",
                measure.vars=c("A", "B"),
                variable.name = "Typ",
                value.name = "Anzahl"
)

data_m
#                    Gene    Typ Anzahl
#                  <char> <fctr>  <num>
#   1:  ENSG00000165623.5      A    150
#   2:  ENSG00000260537.1      A    385
#   3:  ENSG00000171243.7      A    334
#   4:  ENSG00000186399.8      A    487
#   5:  ENSG00000111885.5      A   1100
#  ---                                 
# 382:  ENSG00000265590.5      B    142
# 383:  ENSG00000159212.8      B     90
# 384: ENSG00000100077.10      B    862
# 385:  ENSG00000114487.5      B    151
# 386:  ENSG00000149679.7      B    279

ggplot(data_m, aes(y = Typ, x =Anzahl, fill = Typ)) + 
  geom_boxplot(alpha = 0.8, color = "black",  outlier.size = 2) + 
  labs(
    title = "Verteilung der Anzahl an Knoten pro Gen pro Intervallbaum",
    x = "Anzahl",
    y = "Intervallbaum"
  ) +
  scale_fill_manual(values = c("A" = "steelblue", "B" = "tomato")) + 
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 25), 
    axis.title.x = element_text(size = 21), 
    axis.title.y = element_text(size = 21), 
    axis.text.x = element_text(size = 18, hjust = 1), 
    axis.text.y = element_text(size = 18),
    legend.title = element_text(size = 16, face = "bold"), 
    legend.text = element_text(size = 14) 
  )

  library(dplyr)
medians <- data_m %>%
  group_by(Typ) %>%
  summarise(median_value = median(Anzahl))

# Create the boxplot with median labels
ggplot(data_m, aes(y = Typ, x = Anzahl, fill = Typ)) + 
  geom_boxplot(alpha = 0.8, color = "black", outlier.size = 2) + 
  labs(
    title = "Verteilung der Anzahl an Knoten pro Gen pro Intervallbaum",
    x = "Anzahl",
    y = "Intervallbaum"
  ) +
  scale_fill_manual(values = c("A" = "steelblue", "B" = "tomato")) + 
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 25), 
    axis.title.x = element_text(size = 21), 
    axis.title.y = element_text(size = 21), 
    axis.text.x = element_text(size = 18, hjust = 1), 
    axis.text.y = element_text(size = 18),
    legend.title = element_text(size = 16, face = "bold"), 
    legend.text = element_text(size = 14)
  ) +
  # Add median labels above the boxes
  geom_text(data = medians, aes(x = median_value, y = Typ, label = round(median_value, 2)), 
            color = "black", size = 6, vjust = -7.5)

ggsave("../das_report/plots/trees.png", plot = , width = 14, height = 8, dpi = 300)


max(data$A)
#          A     B
#      <num> <num>
#   1:   150   106
#   2:   385   273
#   3:   334   179
#   4:   487   355
#   5:  1100   507
#  ---            
# 189:   212   142
# 190:   126    90
# 191:  2108   862
# 192:   299   151
# 193:   722   279

data[, Gen := .I]

# Reshape the data to long format for ggplot
data_long <- melt(data, id.vars = "RowID", measure.vars = c("A", "B"), variable.name = "Typ", value.name = "Anzahl")

ggplot(data_long, aes(x = RowID, y = Anzahl, fill = Typ)) +
  geom_bar(stat = "identity", position = "identity", alpha = 0.8) +  # Overlaid bars with transparency
  theme_minimal() +
  labs(title = "Anzahl an Intervallen pro Gen", x = "Gen", y = "Anzahl") +
  scale_fill_manual(values = c("A" = "steelblue", "B" = "tomato")) + 
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 25), 
    axis.title.x = element_text(size = 21), 
    axis.title.y = element_text(size = 21), 
    axis.text.x = element_text(size = 18, hjust = 1), 
    axis.text.y = element_text(size = 18),
    legend.title = element_text(size = 16, face = "bold"), 
    legend.text = element_text(size = 14)
  )

ggsave("../das_report/plots/genetrees.png", plot = , width = 14, height = 8, dpi = 300)
