require("NanoStringNorm")
data("NanoString")

library("ggplot2")

# import NanoString data from tab-separated file
NanoString.mRNA <- read.table("AD_nanoString_data.txt", fill = TRUE, header = TRUE, as.is = TRUE)
str(NanoString.mRNA)

# append tag to duplicated gene names
duplicateNames <- duplicated(NanoString.mRNA$Name)
NanoString.mRNA$Name[duplicateNames] <- paste(NanoString.mRNA$Name[duplicateNames], "-2", sep = "")

# extract sample subsets
annoCols <- colnames(NanoString.mRNA)[c(1:3)]
sampleIDs <- colnames(NanoString.mRNA)[-c(1:3)]
tgIL10Samples <- sampleIDs[c(1:3,7:9)]
tgCtrlSamples <- sampleIDs[c(4:6,10:12)]
wtIL10Samples <- sampleIDs[c(13,17:19)]
wtCtrlSamples <- sampleIDs[c(14:16,20:21)]

tgNS.mRNA <- subset(x = NanoString.mRNA, 
                    select = colnames(NanoString.mRNA) %in% c(annoCols, tgIL10Samples, tgCtrlSamples))
wtNS.mRNA <- subset(x = NanoString.mRNA, 
                    select = colnames(NanoString.mRNA) %in% c(annoCols, wtIL10Samples, wtCtrlSamples))

# build a trait indicator vector for Tg samples
tg.names <- names(tgNS.mRNA)[-c(1:3)]
IL10vsCTRL <- (tg.names %in% tgIL10Samples) + 1 # value 2 indicates IL10 class
trait.tg <- data.frame(
  row.names = tg.names,
  IL10vsCTRL = IL10vsCTRL)

# build a trait indicator vector for WT samples
wt.names <- names(wtNS.mRNA)[-c(1:3)]
IL10vsCTRL <- (wt.names %in% wtIL10Samples) + 1 # value 2 indicates IL10 class
trait.wt <- data.frame(
  row.names = wt.names,
  IL10vsCTRL = IL10vsCTRL)

# define biological replicates for Tg samples
tg.bioReps <- rep("", times = length(tg.names))
tg.bioReps[tg.names %in% tgIL10Samples] <- "Tg.IL10"
tg.bioReps[tg.names %in% tgCtrlSamples] <- "Tg.CTRL"

# define biological replicates for WT samples
wt.bioReps <- rep("", times = length(wt.names))
wt.bioReps[wt.names %in% wtIL10Samples] <- "WT.IL10"
wt.bioReps[wt.names %in% wtCtrlSamples] <- "WT.CTRL"

# compare normalization methods for Tg samples
norm.comp.results.test <- norm.comp(
    x = tgNS.mRNA,
    replicates = tg.bioReps,
    CodeCount.methods = c("none", "sum", "geo.mean"),
    Background.methods = c("none", "mean", "mean.2sd", "max"),
    SampleContent.methods = c("none", "housekeeping.sum", "housekeeping.geo.mean",
                              "total.sum", "low.cv.geo.mean", 
                              "top.mean", "top.geo.mean"),
    OtherNorm.methods = "none",
    histogram = FALSE,
    verbose = TRUE,
    icc.method = "anova")

p <- ggplot(data = norm.comp.results.test, aes(method, icc.anova.results))
p + geom_point()


# compare normalization methods for WT samples
# norm.comp.results.test <- norm.comp(
#   x = wtNS.mRNA,
#   replicates = wt.bioReps,
#   CodeCount.methods = "none",
#   Background.methods = "none",
#   SampleContent.methods = c("none", "housekeeping.sum", "housekeeping.geo.mean",
#                             "top.mean", "top.geo.mean"),
#   OtherNorm.methods = "none",
#   verbose = FALSE)

# recommended normalization for Tg samples + diff exp
tgNS.mRNA.norm <- NanoStringNorm(
  x = tgNS.mRNA,
  anno = NA,
  CodeCount = "geo.mean",
  Background = "mean.2sd",
  SampleContent = "housekeeping.geo.mean",
  round.values = TRUE,
  take.log = TRUE,
  traits = trait.tg)

# recommended normalization for WT samples + diff exp
wtNS.mRNA.norm <- NanoStringNorm(
  x = wtNS.mRNA,
  anno = NA,
  CodeCount = "geo.mean",
  Background = "mean.2sd",
  SampleContent = "housekeeping.geo.mean",
  round.values = TRUE,
  take.log = TRUE,
  traits = trait.wt)

# collect norm outputs
tgNormStats <- tgNS.mRNA.norm$gene.summary.stats.norm
wtNormStats <- wtNS.mRNA.norm$gene.summary.stats.norm

# plot differential expression in Tg samples
pdf("TgIL10_NanoStringNorm_results.pdf")
Plot.NanoStringNorm(
  x = tgNS.mRNA.norm,
  label.best.guess = TRUE,
  plot.type = c("volcano"))
dev.off()

# plot differential expression in WT samples
pdf("WTIL10_NanoStringNorm_results.pdf")
Plot.NanoStringNorm(
  x = wtNS.mRNA.norm,
  label.best.guess = TRUE,
  plot.type = c("volcano"))
dev.off()
