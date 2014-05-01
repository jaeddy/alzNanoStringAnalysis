require("NanoStringNorm")
data("NanoString")

library("ggplot2")

# import NanoString data from tab-separated file
NanoString.mRNA <- read.table("./data/AD_nanoString_data.txt", fill = TRUE, header = TRUE, as.is = TRUE)
str(NanoString.mRNA)

# append tag to duplicated gene names
duplicateNames <- duplicated(NanoString.mRNA$Name)
NanoString.mRNA$Name[duplicateNames] <- paste(NanoString.mRNA$Name[duplicateNames], "-2", sep = "")

# define housekeeping genes
NanoString.mRNA[NanoString.mNRA$Name %in% c("Cltc", "GAPDH", "Hprt1", "Pgk1", "Tubb5"), "Code.Class"] <- "Housekeeping"

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

# check normalization options to find best match for previous data
CodeCountOpts <- c("none", "sum", "geo.mean")
BackgroundOpts <- c("none", "mean", "mean.2sd", "max")

optParams <- c()
for (i in length(CodeCountOpts)) {
  for (j in length(BackgroundOpts)) {
    tgNS.mRNA.norm <- NanoStringNorm(
      x = tgNS.mRNA,
      anno = NA,
      CodeCount = CodeCountOpts[i],
      Background = BackgroundOpts[j],
      SampleContent = "housekeeping.geo.mean",
      round.values = FALSE,
      take.log = FALSE,
      verbose = FALSE)
    
    testVal <- tgNS.mRNA.norm$normalized.data$X20140224_il10_brains_custom2_01.RCC[1]
    if (abs(testVal - 444.17) < 10) {
      optParams <- c(optParams, CodeCountOpts[i], BackgroundOpts[j])
      optParams
    }
  }
}
  

# compare normalization methods for Tg samples
tg.norm.comp.results.test <- norm.comp(
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

# compare normalization methods for WT samples
wt.norm.comp.results.test <- norm.comp(
  x = wtNS.mRNA,
  replicates = wt.bioReps,
  CodeCount.methods = c("none", "sum", "geo.mean"),
  Background.methods = c("none", "mean", "mean.2sd", "max"),
  SampleContent.methods = c("none", "housekeeping.sum", "housekeeping.geo.mean",
                            "total.sum", "low.cv.geo.mean", 
                            "top.mean", "top.geo.mean"),
  OtherNorm.methods = "none",
  histogram = FALSE,
  verbose = TRUE,
  icc.method = "anova")

# plot ICC vs. sample content correction method for Tg data
p <- ggplot(data = tg.norm.comp.results.test, aes(SampleContent.method, icc.anova.results))
p + geom_jitter(aes(color = factor(CodeCount.method),
                    shape = factor(Background.method)),
               size = 3)

# plot ICC vs. background correction method for top.geo.mean subset of Tg data
p2 <- ggplot(data = subset(tg.norm.comp.results.test, SampleContent.method == "top.geo.mean"),
             aes(Background.method, icc.anova.results))
p2 + geom_jitter(aes(color = factor(CodeCount.method),
                     size = factor(cv.bio2tech.ratio)))

# plot ICC vs. sample content correction method for WT data
p <- ggplot(data = wt.norm.comp.results.test, aes(SampleContent.method, icc.anova.results))
p + geom_jitter(aes(color = factor(CodeCount.method),
                    shape = factor(Background.method)),
                size = 3)

# plot ICC vs. background correction method for top.geo.mean subset of Tg data
p2 <- ggplot(data = subset(wt.norm.comp.results.test, SampleContent.method == "top.geo.mean"),
             aes(Background.method, icc.anova.results))
p2 + geom_jitter(aes(color = factor(CodeCount.method),
                     size = factor(cv.bio2tech.ratio)))

# plot code count correction vs. background correction method scaled according to ICC
p3 <- ggplot(data = subset(wt.norm.comp.results.test, SampleContent.method == "top.geo.mean"),
             aes(Background.method, CodeCount.method))
p3 + geom_point(aes(size = factor(icc.anova.results),
                    alpha = factor(cv.bio2tech.ratio)))

# recommended normalization for Tg samples + diff exp (sum, max, top.geo.mean)
tgNS.mRNA.norm <- NanoStringNorm(
  x = tgNS.mRNA,
  anno = NA,
  CodeCount = "geo.mean",
  Background = "max",
  SampleContent = "housekeeping.geo.mean",
  round.values = FALSE,
  take.log = FALSE,
  traits = trait.tg)

tgNS.mRNA.norm$normalized.data$X20140224_il10_brains_custom2_01.RCC[1]


# recommended normalization for WT samples + diff exp
wtNS.mRNA.norm <- NanoStringNorm(
  x = wtNS.mRNA,
  anno = NA,
  CodeCount = "sum",
  Background = "max",
  SampleContent = "top.geo.mean",
  round.values = TRUE,
  take.log = TRUE,
  traits = trait.wt)

# collect norm outputs
tgNormStats <- tgNS.mRNA.norm$gene.summary.stats.norm
tgNormData <- tgNS.mRNA.norm$normalized.data;
wtNormStats <- wtNS.mRNA.norm$gene.summary.stats.norm
wtNormStats <- wtNS.mRNA.norm$normalized.data;

# plot differential expression in Tg samples
pdf("./results/TgIL10_NanoStringNorm_results.pdf")
Plot.NanoStringNorm(
  x = tgNS.mRNA.norm,
  label.best.guess = TRUE,
  plot.type = c("cv", "mean.sd", "RNA.estimates", "volcano", "missing",
                "norm.factors", "positive.controls"))
dev.off()

# plot differential expression in WT samples
pdf("./results/WTIL10_NanoStringNorm_results.pdf")
Plot.NanoStringNorm(
  x = wtNS.mRNA.norm,
  label.best.guess = TRUE,
  plot.type = c("cv", "mean.sd", "RNA.estimates", "volcano", "missing",
                "norm.factors", "positive.controls"))
dev.off()
