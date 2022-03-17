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
library(here)

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

makeCircPlot <- function(goi)
{
  conditions <- c("TissueBrown:GenotypeControl",
                  "TissueWhite:GenotypeControl",
                  "TissueBrown:GenotypeKO",
                  "TissueWhite:GenotypeKO")

  # Calculate expected values
  predictedValues <- lapply(conditions, calcExpected, model = model, genes = goi) %>%
    do.call(what = "rbind")

  ### Visualize it
  predictedValues[, c("Tissue", "Genotype"):=tstrsplit(condition, ":")]
  ens2symbol <- structure(names(goi), names = goi)
  predictedValues[, SYMBOL:=ens2symbol[ENSEMBL]]
  predictedValues[, Genotype:=str_remove(Genotype, "Genotype")]
  predictedValues[, Tissue:=str_remove(Tissue, "Tissue")]

  realValues <- cpm(dgeLists$circadian, log = TRUE)[goi, , drop = FALSE] %>% reshape2::melt() %>%
    data.table()
  setnames(realValues, c("ENSEMBL", "id", "expression"))
  setkey(realValues, id)

  colDat <- metadata
  setkey(colDat, id)

  pDReal <- colDat[realValues]
  pDReal[, time:=as.integer(sub("ZT", "", Timepoint))]
  pDReal[, SYMBOL:=ens2symbol[ENSEMBL]]
  pDReal[Tissue == "Brown" & Genotype == "KO", Genotype := "FANKO BAT"]
  pDReal[Tissue == "White" & Genotype == "KO", Genotype := "FANKO eWAT"]

  colScale <- scale_color_manual(values = c("Control" = "#606060", "FANKO BAT" = "#A00000", "FANKO eWAT" = "#94641F"), #values = c("Control" = "#1f78b4", "KO" = "#33a02c"),
                                 name = element_blank())

  predictedValues[Tissue == "Brown" & Genotype == "KO", Genotype := "FANKO BAT"]
  predictedValues[Tissue == "White" & Genotype == "KO", Genotype := "FANKO eWAT"]

  ggplot(predictedValues, aes(x = time, y = expression, colour = Genotype)) +
    geom_line() +
    geom_point(data = pDReal) +
    facet_grid(Tissue ~ SYMBOL) +
    colScale +
    scale_x_continuous(name = NULL, limits = c(0, 24), expand = c(0, 0)) +
    scale_y_continuous(name = element_blank()) +
    theme_bw(base_size = 9) +
    guides(colour = "none") +
    theme(strip.background.y = element_blank(), strip.text.y = element_blank())
}

load(here("out/limorhyde/limorhydeEdgerResults.RData"))
coreClockSymbol <- c("Clock", "Arntl", "Npas2", "Per1", "Per2", "Per3", "Cry1", "Cry2", "Nr1d1", "Nr1d2", "Rora", "Rorb", "Rorc", "Dbp")
forConversion <- limorhydeEdgerResults$CircadianBrownControl
coreClock <- forConversion[coreClockSymbol, ENSEMBL, on = "SYMBOL"]
names(coreClock) <- coreClockSymbol

load(here("out/limorhyde/dgeLists.RData"))
model <- glmQLFit(dgeLists$circadian)

load(here("data/metadata.rda"))

cmToInches <- 0.393701
makeAndSavePlots <- function(gene)
{
  p <- makeCircPlot(gene)
  ggsave(plot = p, filename = here(paste0("out/expression-plots/", names(gene), "_expression.png")), height = 7.5, width = 3.5, units = "cm")
}
for (i in seq_along(coreClock)){
  makeAndSavePlots(coreClock[i])
}
