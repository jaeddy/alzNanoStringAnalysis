filePath <- "~/Google Drive/MyCode/MatlabCode/Expression Analysis/AD Expression/ipaCanonPaths.csv"

canonPaths <- read.csv(filePath, header = FALSE, stringsAsFactors = FALSE, na.strings = "")
pathNames <- canonPaths$V1
pathGenes <- canonPaths[ ,-1]
totalGenes <- length(unique(unlist(pathGenes)))-1

targetPath <- pathNames[13]


targetPathIdx <- pathNames %in% targetPath
targetPathGenes <- pathGenes[targetPathIdx, ]
targetPathGenes <- targetPathGenes[!is.na(targetPathGenes)]

overlapMat <- data.frame(matrix(nrow = length(pathNames), ncol = 4))
colnames(overlapMat) <- c("pathway", "overlap", "size", "p-value")

for(i in 1:length(pathNames)){
  queryPath <- pathNames[i]
  queryPathGenes <- pathGenes[pathNames %in% queryPath, ]
  queryPathGenes <- queryPathGenes[!is.na(queryPathGenes)]
  
  overlap <- sum(queryPathGenes %in% targetPathGenes)
  total <- length(queryPathGenes)
  
  pval <- phyper(overlap-1, length(targetPathGenes), totalGenes, total,
                 lower.tail = FALSE)
  overlapMat[i, ] <-  c(queryPath, overlap, total, pval)
  
}
