library(magrittr)
library(here)
library(data.table)
library(stringr)
library(topGO)
library(openxlsx)
load(file = here("out/limorhyde/limorhydeEdgerResults.RData"))

names(limorhydeEdgerResults) <- str_remove(names(limorhydeEdgerResults), "Circadian")

allGenes <- factor(rep(1, nrow(limorhydeEdgerResults[[1]])),
                   levels = c(0, 1))
names(allGenes) <- limorhydeEdgerResults[[1]]$ENSEMBL

goData <- new("topGOdata",
              description = "Simple session",
              ontology = "BP",
              allGenes = allGenes,
              nodeSize = 10,
              annot = annFUN.org,
              mapping = "org.Mm.eg.db",
              ID = "ensembl"
)


prepareGeneSets <- function(cond1, cond2, resultsList){
  getLists <- function(x, y, xname, yname){
    out <- list(
      setdiff(x, y),
      setdiff(y, x),
      intersect(x, y)
    )

    names(out) <- c(xname, yname, paste(xname, yname, sep = "_and_"))
    out
  }
  fdrCutoff <- 0.01
  geneLists <- getLists(resultsList[[cond1]][FDR < fdrCutoff, ENSEMBL],
                        resultsList[[cond2]][FDR < fdrCutoff, ENSEMBL],
                        cond1,
                        cond2)

  geneLists
}


tester <- function(geneList, testName, goData){
  allGeneList <- geneScore(goData, use.names = TRUE)
  newAllGenes <- names(allGeneList) %in% geneList %>%
    as.numeric %>%
    factor(levels = c(0, 1))
  names(newAllGenes) <- names(allGeneList)
  goData <- updateGenes(goData, geneList = newAllGenes)
  description(goData) <- testName
  list(goRes = runTest(goData, algorithm = "weight01", statistic = "fisher"),
       goData = goData,
       genes = geneList
  )
}


formatOutput <- function(go_list){
  goRes <- go_list[["goRes"]]
  goData <- go_list[["goData"]]
  signif_genes <- go_list[["genes"]]
  nGo <- length(score(goRes))
  out <- GenTable(goData, weight01_PValue = goRes, topNodes = nGo, numChar = 10^9)

  # The P-values from weight01 are not independent, and does not work with
  # FDR correction. The manual suggests treating them as already corrected.
  # The P-value distribution seems to support this notion.
  #out[, FDR:=p.adjust(weight01_PValue, method = "fdr")]

  all_terms <- genesInTerm(goData, out$GO.ID)
  all_ens <- Reduce(union, all_terms)
  ens2symbol <- select(org.Mm.eg.db, all_ens, columns = c("ENSEMBL", "SYMBOL"), keytype = "ENSEMBL")
  setDT(ens2symbol)
  ens2symbol <- ens2symbol[, .(SYMBOL = SYMBOL[[1]]), by = "ENSEMBL"] # Keep only one ENSEMBL to SYMBOL
  setkey(ens2symbol, ENSEMBL)
  ens2symbol <- ens2symbol[signif_genes, ]
  setkey(ens2symbol, ENSEMBL)
  genes_in_term <- lapply(all_terms, \(x) ens2symbol[x, paste0(SYMBOL, collapse = ", "), nomatch = FALSE])
  setDT(out)
  out[, genes := NA_character_]
  setkey(out, GO.ID)
  out[names(genes_in_term), genes := unlist(genes_in_term)]
  out[, weight01_PValue := as.numeric(weight01_PValue)]
  out <- out[order(weight01_PValue, decreasing = FALSE), ]
  out
}


####### Gene ontologies of overlaps
#############################################
goFromEdgeR <- function(cond1, cond2, edgerResults, goData){
  sets <- prepareGeneSets(cond1,
                          cond2,
                          edgerResults
  )
  goResults <- mapply(tester, sets, names(sets), MoreArgs = list(goData = goData), SIMPLIFY = FALSE)
  lapply(goResults, formatOutput)
}
brownRes <- goFromEdgeR("BrownControl", "BrownKO", edgerResults = limorhydeEdgerResults, goData = goData)
whiteRes <- goFromEdgeR("WhiteControl", "WhiteKO", edgerResults = limorhydeEdgerResults, goData = goData)
write.xlsx(c(brownRes, whiteRes), file = here("out/ontologies/Supplemental table 2.xlsx"))


####### Gene ontologies of loss of rythmicity
##############################################
cutoff <- 0.01
loss_in_brown <- setdiff(limorhydeEdgerResults$BrownControl[FDR < cutoff, ENSEMBL],
                         limorhydeEdgerResults$BrownKO[FDR < cutoff, ENSEMBL])
loss_in_white <- setdiff(limorhydeEdgerResults$WhiteControl[FDR < cutoff, ENSEMBL],
                         limorhydeEdgerResults$WhiteKO[FDR < cutoff, ENSEMBL])
loss_in_both <- intersect(loss_in_brown, loss_in_white)

lossGoResults <- tester(loss_in_both,
                        "Loss in both",
                        goData = goData)
lossGoResults <- formatOutput(lossGoResults)
write.xlsx(lossGoResults, file = here("out/ontologies/Supplemental table 1.xlsx"))
