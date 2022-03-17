library(openxlsx)
library(data.table)
library(here)
library(magrittr)
library(stringr)

#### First dataset

load(here("data/metabolite_data1.rda"))

modelBuilder <- function(metabolite, dataset){
  model <- aov(abundance ~ (Genotype*Timepoint), dataset[metabolite])
  summ <- summary(model)[[1]]
  names <- c("GenotypeKO", "TimepointNight", "Interaction")
  structure(c(model$coefficients[-1], summ[-4, 5]), names = paste0(names, rep(c("_Estimate", "_PValue"), each = 3)))
}

metabolites <- unique(metabolite_data1$Label)
names(metabolites) <- metabolites

resultsPValues <- lapply(metabolites, . %>% modelBuilder(dataset = metabolite_data1)) %>%
  do.call(what = "rbind") %>%
  as.data.table(keep.rownames = TRUE)

setnames(resultsPValues, c("rn"), c("Label"))
resultsPValues[, GenotypeKO_FDR:=p.adjust(GenotypeKO_PValue, method = "fdr")]
resultsPValues[, TimepointNight_FDR:=p.adjust(TimepointNight_PValue, method = "fdr")]
resultsPValues[, Interaction_FDR:=p.adjust(Interaction_PValue, method = "fdr")]

out <- list(
  KO_effect = resultsPValues[, .(Label, GenotypeKO_Estimate, GenotypeKO_PValue, GenotypeKO_FDR)],
  Night_effect = resultsPValues[, .(Label, TimepointNight_Estimate, TimepointNight_PValue, TimepointNight_FDR)],
  Interaction_effect = resultsPValues[, .(Label, Interaction_Estimate, Interaction_PValue, Interaction_FDR)]
)
setkey(out$KO_effect, GenotypeKO_PValue)
setkey(out$Night_effect, TimepointNight_PValue)
setkey(out$Interaction_effect, Interaction_PValue)

write.xlsx(out, file = here("out/metabolites/metabolites_first.xlsx"))
rm(metabolite_data1, out, resultsPValues, metabolites, modelBuilder)

#### Second dataset
load(here("data/metabolite_data2.rda"))

model_builder <- function(dataset, model){
  model <- aov(model, dataset)
  summ <- summary(model)[[1]]

  model_coefs <- coef(model)[-1]
  pvals <- summ[-nrow(summ), 5]

  col_names <- names(model_coefs) |>
    str_replace_all(":", "_") |>
    str_remove_all("Genotype|Timepoint|Diet")

  out <- as.list(c(model_coefs, pvals))
  names(out) <- paste0(col_names,
                       rep(c("_Estimate", "_PValue"),
                           each = length(col_names)))
  out
}

three_way <- metabolite_data2[, model_builder(.SD, abundance ~ (Genotype*Timepoint*Diet)),
                             by = "Label",
                             .SDcols = c("abundance", "Genotype", "Timepoint", "Diet")]

two_way_chow <- metabolite_data2["Chow", model_builder(.SD, abundance ~ (Genotype*Timepoint)),
                                by = "Label",
                                .SDcols = c("abundance", "Genotype", "Timepoint")]

two_way_hfd <- metabolite_data2["HFD", model_builder(.SD, abundance ~ (Genotype*Timepoint)),
                               by = "Label",
                               .SDcols = c("abundance", "Genotype", "Timepoint")]


formatter <- function(x) {
  sheets <- str_remove(str_subset(colnames(x), "Estimate"), "_Estimate")
  helper <- function(y) {
    out <- x[, c("Label", paste0(y, c("_Estimate", "_PValue"))), with = FALSE]
    setnames(out, str_remove(colnames(out), paste0(y, "_")))
    out[, FDR := p.adjust(PValue, "fdr")]
    out[order(PValue, decreasing = FALSE)]
  }
  names(sheets) <- sheets
  lapply(sheets, helper)
}

writexl::write_xlsx(formatter(three_way), path = here("out/metabolites/metabolites_second_three_way.xlsx"))
writexl::write_xlsx(formatter(two_way_chow), path = here("out/metabolites/metabolites_second_two_way_chow.xlsx"))
writexl::write_xlsx(formatter(two_way_hfd), path = here("out/metabolites/metabolites_second_two_way_hfd.xlsx"))
