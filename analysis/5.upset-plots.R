library(data.table)
library(ggplot2)
library(magrittr)
library(stringr)
library(UpSetR)
library(here)
library(export)
load(here("out/limorhyde/limorhydeEdgerResults.RData"))

#####
# upsetPlot
contrastsOfInterest <- c("CircadianBrownControl",
                         "CircadianBrownKO",
                         "CircadianWhiteControl",
                         "CircadianWhiteKO")

calculateOverlaps <- function(listOfContrasts, fdrCutoff)
{
  genes <- listOfContrasts[[1]]$ENSEMBL
  highlightOverlaps <- function(x)
  {
    genes %in% x[FDR < fdrCutoff, ENSEMBL] %>% as.integer
  }
  overlapData <- lapply(listOfContrasts, highlightOverlaps) %>%
    as.data.table()
  overlapData[, ENSEMBL:=genes]
  overlapData[rowSums(overlapData[, !"ENSEMBL"]) > 0]
}

upsetPlotData <- calculateOverlaps(
  listOfContrasts = limorhydeEdgerResults[contrastsOfInterest],
  fdrCutoff = 0.01
)

setnames(upsetPlotData,
         c("CircadianBrownControl",
           "CircadianBrownKO",
           "CircadianWhiteControl",
           "CircadianWhiteKO"),
         c("WT_BAT",
           "KO_BAT",
           "WT_eWAT",
           "KO_eWAT")
)

intersectionsOfInterest = list(
  list("KO_eWAT"),
  list("WT_eWAT"),
  list("KO_BAT"),
  list("WT_BAT"),
  list("WT_BAT", "KO_BAT"),
  list("WT_eWAT", "KO_eWAT"),
  list("WT_BAT", "WT_eWAT"),
  list("KO_BAT", "KO_eWAT"),
  list("KO_BAT", "KO_eWAT", "WT_BAT", "WT_eWAT")
)

p <- upset(upsetPlotData,
           nsets = 4,
           intersections = intersectionsOfInterest,
           keep.order = TRUE
)
p

cmToInches <- 0.393701
png(here("out/upset-plot/upsetPlot.png"), width = 10.95*cmToInches, height = 6.55*cmToInches, units = "in", res = 300)
p
dev.off()

graph2ppt(x = p, file = here("out/upset-plot/upsetPlot.pptx"))
# It is difficult to remove the underscore from the name. Do it in powerpoint instead
