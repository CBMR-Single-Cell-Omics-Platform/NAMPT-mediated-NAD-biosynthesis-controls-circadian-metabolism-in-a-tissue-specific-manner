library(data.table)
library(ggplot2)
library(magrittr)
library(stringr)
library(edgeR)
library(SummarizedExperiment)
library(here)
library(limorhyde)
library(openxlsx)

load(here("data/counts.rda"))
load(here("data/metadata.rda"))

prepareForEdgeR <- function(counts, metadata, differential)
{
  y <- DGEList(counts[, -c(1:8)], genes = counts[, c(1:8)], samples = metadata, group = metadata$group)
  rownames(y) <- y$genes$ENSEMBL
  metadata <- copy(metadata)
  metadata[, ZT:=str_sub(Timepoint, 3) %>% as.integer]
  limo <- limorhyde(metadata$ZT, 'ZT_', period = 24)
  metadata <- cbind(metadata, limo)
  metadata[, Genotype:=factor(Genotype, levels = c("Control", "KO"))]

  design <- model.matrix(~ 0 + Tissue:Genotype + Tissue:Genotype:(ZT_cos + ZT_sin), data=metadata)
  if (differential)
  {
    design[, "TissueBrown:GenotypeControl:ZT_cos"] <- design[, "TissueBrown:GenotypeControl:ZT_cos"] + design[, "TissueBrown:GenotypeKO:ZT_cos"]
    design[, "TissueBrown:GenotypeControl:ZT_sin"] <- design[, "TissueBrown:GenotypeControl:ZT_sin"] + design[, "TissueBrown:GenotypeKO:ZT_sin"]
    design[, "TissueWhite:GenotypeControl:ZT_cos"] <- design[, "TissueWhite:GenotypeControl:ZT_cos"] + design[, "TissueWhite:GenotypeKO:ZT_cos"]
    design[, "TissueWhite:GenotypeControl:ZT_sin"] <- design[, "TissueWhite:GenotypeControl:ZT_sin"] + design[, "TissueWhite:GenotypeKO:ZT_sin"]

    colnames(design)[c(7,8,11,12)] <- paste0(colnames(design)[c(7,8,11,12)], "_differential")
  }
  colnames(design)
  is.fullrank(design)

  y <- y[filterByExpr(y, design = design), , keep.lib.sizes=FALSE]
  y <- calcNormFactors(y)
  y <- estimateDisp(y, design)
  y
}

edgerTester <- function(x)
{
  fit <- glmQLFit(x)

  helper <- . %>%
    topTags(n = Inf, p.value = 1) %>%
    extract2("table") %>%
    as.data.table
  list(
    BrownControl = glmQLFTest(glmfit = fit, coef = c(5, 9)) %>% helper,
    WhiteControl = glmQLFTest(glmfit = fit, coef = c(6, 10)) %>% helper,
    BrownKO = glmQLFTest(glmfit = fit, coef = c(7, 11)) %>% helper,
    WhiteKO = glmQLFTest(glmfit = fit, coef = c(8, 12)) %>% helper,
    BrownStaticKO = glmQLFTest(glmfit = fit, contrast = c(-1,0,1,0,0,0,0,0,0,0,0,0)) %>% helper,
    WhiteStaticKO = glmQLFTest(glmfit = fit, contrast = c(0,-1,0,1,0,0,0,0,0,0,0,0)) %>% helper
  )
}

dgeLists <- list(
  circadian = prepareForEdgeR(counts, metadata, differential = FALSE),
  differentialCircadian = prepareForEdgeR(counts, metadata, differential = TRUE)
)

edgerRes <- lapply(dgeLists, edgerTester)
limorhydeEdgerResults <- list(
  CircadianBrownControl = edgerRes$circadian$BrownControl,
  CircadianWhiteControl = edgerRes$circadian$WhiteControl,
  CircadianBrownKO = edgerRes$circadian$BrownKO,
  CircadianWhiteKO = edgerRes$circadian$WhiteKO,
  CircadianBrownKODifferential = edgerRes$differentialCircadian$BrownKO,
  CircadianWhiteKODifferential = edgerRes$differentialCircadian$WhiteKO,
  StaticBrownKO = edgerRes$circadian$BrownStaticKO,
  StaticWhiteKO = edgerRes$circadian$WhiteStaticKO
)

save(limorhydeEdgerResults, file = here("out/limorhyde/limorhydeEdgerResults.RData"))
write.xlsx(limorhydeEdgerResults, file = here("out/limorhyde/limorhydeEdgerResults.xlsx"))

save(dgeLists, file = here("out/limorhyde/dgeLists.RData"))

expression <- cpmByGroup(dgeLists$circadian)
expression_sample <- cpm(dgeLists$circadian, normalized.lib.sizes = TRUE)

expression <- cbind(expression, dgeLists$circadian$genes[rownames(expression), c("SYMBOL", "GENENAME")])
setDT(expression, keep.rownames = TRUE)
setnames(expression, "rn", "ENSEMBL")

expression_sample <- cbind(expression_sample, dgeLists$circadian$genes[rownames(expression_sample), c("SYMBOL", "GENENAME")])
setDT(expression_sample, keep.rownames = TRUE)
setnames(expression_sample, "rn", "ENSEMBL")

write.xlsx(expression, file = here("out/expressions/expressions.xlsx"), asTable = TRUE)
write.xlsx(expression_sample, file = here("out/expressions/expressions_sample.xlsx"), asTable = TRUE)
