library(viridis)
library(pheatmap)
library(data.table)
library(ggplot2)
library(magrittr)
library(stringr)
library(edgeR)
library(limorhyde)
library(circular)
library(cetcolor)
library(export)
library(here)

load(here("out/limorhyde/dgeLists.RData"))
model <- glmQLFit(dgeLists$circadian)

calcExpected <- function(condition, model, genes, timeVariable = "ZT_", nPoints = 1000, timeRange = c(0, 24))
{
  dt <- model$coefficients %>% as.data.table(keep.rownames = TRUE)
  setnames(dt, "rn", "ENSEMBL")
  setkey(dt, "ENSEMBL")

  sinVal <- paste0(condition, ":", timeVariable, "sin")
  cosVal <- paste0(condition, ":", timeVariable, "cos")

  zt <- data.table(ZT = seq(timeRange[1], timeRange[2], length.out = nPoints))
  zt[, c("ZT_cos", "ZT_sin"):=limorhyde(zt, 'ZT', period = 24)]

  dt <- dt[genes, c("ENSEMBL", condition, sinVal, cosVal), with = FALSE]
  setnames(dt, c("ENSEMBL", "grp", "sin", "cos"))
  calcVals <- function(cosVal, sinVal, grp)
  {
    zt$ZT_cos * cosVal + zt$ZT_sin * sinVal + grp
  }
  out <- dt[, .(expression = calcVals(cos, sin, grp), time = zt$ZT, condition = condition), by = "ENSEMBL"]
  out[, expression:=(expression + log(1e6))/log(2)]
  out
}

makeHeatmap <- function(genes, value.var, show_rownames = FALSE, clusterRows = TRUE)
{
  predictedValues <- lapply(conditions,
                            calcExpected,
                            model = model,
                            genes = genes) %>%
    do.call(what = "rbind")
  predictedParam <- predictedValues[, .(peak = x <- time[which.max(expression)],
                                        trough = y <- time[which.min(expression)],
                                        amplitude = expression[time == x] - expression[time == y],
                                        meanExpression = mean(expression)
  ),
  by = c("ENSEMBL", "condition")]
  symbol2ens <- structure(names(genes), names = genes)
  predictedParam[, SYMBOL:=symbol2ens[ENSEMBL]]
  predictedParam[, condition:=str_remove_all(condition, "Tissue|Genotype")]

  phasesOfInterest <- dcast(predictedParam[!is.na(SYMBOL), ],
                            SYMBOL ~ condition,
                            value.var = value.var,
                            fun.aggregate = mean)


  phasesOfInterestMatrix <- as.matrix(phasesOfInterest[, !"SYMBOL"])
  rownames(phasesOfInterestMatrix) <- phasesOfInterest$SYMBOL
  colnames(phasesOfInterestMatrix) <- str_replace(colnames(phasesOfInterestMatrix), ":", ", ")
  if (!clusterRows)
  {
    phasesOfInterestMatrix <- phasesOfInterestMatrix[names(genes), ]
  }

  if (value.var == "peak")
  {
    hmCols <- cet_pal("c2", n = 24)
    legend = FALSE
    distFn <- . %>% dist.circular(method = "angularseparation") #geodesic
    phasesOfInterestMatrix <- (phasesOfInterestMatrix * 2 * pi / 24)
  } else if (value.var == "amplitude")
  {
    hmCols <- viridis_pal()(10)
    legend = FALSE
    distFn <- . %>% dist(method = "euclidean")
    phasesOfInterestMatrix <- t(scale(t(phasesOfInterestMatrix)))
  }

  myDendro <- function(x, order)
  {
    x %>%
      distFn %>%
      hclust(method = "ward.D2") %>%
      as.dendrogram %>%
      reorder(order, agglo.FUN = max) %>%
      rev %>%
      as.hclust
  }
  if (clusterRows)
  {
    rowDendro_goi <- myDendro(phasesOfInterestMatrix, rownames(phasesOfInterestMatrix) %>% seq_along)
  } else
  {
    rowDendro_goi <- FALSE
  }
  #colDendro_goi <- myDendro(t(phasesOfInterestMatrix), colnames(phasesOfInterestMatrix) %>% seq_along)


  pheatmap(phasesOfInterestMatrix,
           show_rownames = show_rownames,
           show_colnames = TRUE,
           color = hmCols,
           cluster_rows = rowDendro_goi,
           cluster_cols = FALSE,
           treeheight_row = 25,
           legend = legend,
           border_color = NA,
           labels_col = "",
           scale = "none",
           fontface = "italic")
}

### See what the names of the goups are:
# The non-circadian parts are:
conditions <- c("TissueBrown:GenotypeControl",
                "TissueWhite:GenotypeControl",
                "TissueBrown:GenotypeKO",
                "TissueWhite:GenotypeKO")

coreClockSymbol <- c("Clock", "Arntl", "Npas2", "Per1", "Per2", "Per3", "Cry1", "Cry2", "Nr1d1", "Nr1d2", "Rora", "Rorb", "Rorc", "Dbp")
load(here("out/limorhyde/limorhydeEdgerResults.RData"))
forConversion <- limorhydeEdgerResults$CircadianBrownControl
coreClock <- forConversion[coreClockSymbol, ENSEMBL, on = "SYMBOL"]
names(coreClock) <- coreClockSymbol

cmToInches <- 0.393701
png(here("out/heatmaps/S1L.png"), height = 10.95*cmToInches, width = 6.55*cmToInches, units = "in", res = 300)
makeHeatmap(coreClock, value.var = "peak", show_rownames = TRUE, clusterRows = FALSE)
dev.off()

png(here("out/heatmaps/1N.png"), height = 10.95*cmToInches, width = 6.55*cmToInches, units = "in", res = 300)
makeHeatmap(coreClock, value.var = "amplitude", show_rownames = TRUE, clusterRows = FALSE)
dev.off()

### Heatmaps of all genes
allGenes <- lapply(limorhydeEdgerResults[1:4], extract, FDR < 0.01, ENSEMBL) %>% Reduce(f = "union")
names(allGenes) <- forConversion[allGenes, SYMBOL, on = "ENSEMBL"]

png(here("out/heatmaps/2E.png"), width = 8*cmToInches, height = 8*cmToInches, units = "in", res = 300)
makeHeatmap(allGenes, value.var = "peak")
dev.off()

png(here("out/heatmaps/2D.png"), height = 8*cmToInches, width = 8*cmToInches, units = "in", res = 300)
makeHeatmap(allGenes, value.var = "amplitude")
dev.off()
