
filePath <- "~/Google Drive/MyCode/RCode/AD_nanoString/data/nTgIL10degs.txt"
gene_list <- read.table(filePath, header = TRUE, sep = "\t")

require(ggplot2)
## Highlight genes that have an absolute fold change > 2 and a p-value < Bonferroni cut-off
isSig <- gene_list$Q.value < 0.05
isFC <- abs(gene_list$Log.ratio) > 2
sigLevels <- numeric(length(gene_list$Log.ratio))
sigLevels[isSig] <- 1
sigLevels[isSig & isFC] <- 2
## gene_list$threshold <- as.factor(abs(gene_list$Log.ratio) > 2 & gene_list$Q.value < 0.05)
gene_list$threshold <- as.factor(sigLevels)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
               "#D55E00", "#CC79A7")



## Construct the plot object
g <- ggplot(data = gene_list, aes(x = Log.ratio, y = -log10(P.value), colour = threshold)) +
  geom_point(alpha = 1, size = 2) +
  theme(legend.position = "none") +
  xlim(c(-3, 10)) + ylim(c(0, 8)) +
  xlab("log2 fold change") + ylab("-log10 p-value") +
  scale_colour_manual(values=cbPalette)
g
## Add labels
interestGenes <- c("IL10", "APOE", "Ccl8", "CCL5", "Ccl2", "FCGR3A")
geneIdx <- as.character(gene_list$IPA.Name) %in% interestGenes
genes <- gene_list[geneIdx, ]
g <- g + geom_text(aes(x = genes$Log.ratio-0.4, y = -log10(genes$P.value)+0.25,
                  label = genes$IPA.Name, size = 1.2), colour="black")
s
ggsave(filename = "nTgVolcano.eps", plot = g,
       path = "~/Google Drive/MyCode/RCode/AD_nanoString/results/",
       width = 8.2, height = 8.2, units = "cm", dpi = 300)