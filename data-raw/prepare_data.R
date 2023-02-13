#library(org.Mm.eg.db)
library(data.table)
library(here)
library(magrittr)
library(readxl)
library(usethis)

counts <- fread(here("data-raw/featureCounts.csv.gz"))
use_data(counts, overwrite = TRUE)

#### Construct metadata
metadata <- read_xlsx(path = here("data-raw/metadata.xlsx"))
setDT(metadata)
use_data(metadata, overwrite = TRUE)

#### Prepare metabolites
prepare_bat_normalized <- function() {
  file <- here("data-raw/BAT26062020_both metaboanalyst auto scaled and log2.csv")
  abundance <- fread(file, skip = 1)

  out <- suppressWarnings(melt.data.table(abundance,
                                          id.vars = "Label",
                                          variable.factor = FALSE))
  setnames(out, c("Label", "Sample", "abundance"))
  out[, c("Genotype", "Timepoint"):=tstrsplit(Sample, "_")]
  out[, Genotype:=factor(Genotype, levels = c("WT", "KO"))]
  out[, Timepoint:=factor(Timepoint, levels = c("day", "night"))]
  setkey(out, "Label")
  out[]
}

format_bat_ik <- function(){
  file <- here("data-raw/BAT normalized for IK.csv")
  sample_names <- t(fread(file, nrows = 2))
  sample_names[1,1] <- "Label"
  metadata <- data.table(sample_names[-1, ])
  setnames(metadata, c("Sample", "tmp"))
  metadata[, c("Genotype", "Diet", "Timepoint") := tstrsplit(tmp, "_")]
  metadata[, Genotype := factor(fifelse(Genotype == "p", "KO", "WT"),
                                levels = c("WT", "KO"))]
  metadata[, Diet := factor(fifelse(Diet == "F", "HFD", "Chow"),
                            levels = c("Chow", "HFD"))]
  metadata[, Timepoint := factor(fifelse(Timepoint == "N", "Noon", "Midnight"),
                                 levels = c("Noon", "Midnight"))]
  metadata[, tmp := NULL]

  metabolite_data <- fread(file, skip = 1)
  setnames(metabolite_data, sample_names[,1])
  metabolite_data <- metabolite_data[, lapply(.SD, as.numeric), by = "Label"]

  metabolite_data <- melt(metabolite_data,
                          id.vars = "Label",
                          variable.name = "Sample",
                          value.name = "abundance")
  metabolite_data <- merge.data.table(metabolite_data, metadata, by = "Sample")
  setkey(metabolite_data, Diet)
  metabolite_data[]
}

format_ewat <- function() {
  file <- here("data-raw/eWat pos neg 26062020 log2.csv")
  sample_names <- t(fread(file, nrows = 2))
  sample_names[1,1] <- "Label"
  metadata <- data.table(sample_names[-1, ])
  setnames(metadata, c("Sample", "tmp"))

  metadata[, c("Genotype", "Timepoint") := tstrsplit(tmp, "_")]
  metadata[, Genotype := factor(Genotype, levels = c("WT", "KO"))]
  metadata[, Timepoint := factor(Timepoint, levels = c("day", "night"))]
  metadata[, tmp := NULL]

  metabolite_data <- fread(file, skip = 1)
  setnames(metabolite_data, sample_names[,1])
  metabolite_data <- metabolite_data[, lapply(.SD, as.numeric), by = "Label"]

  metabolite_data <- melt(metabolite_data,
                          id.vars = "Label",
                          variable.name = "Sample",
                          value.name = "abundance")
  metabolite_data <- merge.data.table(metabolite_data, metadata, by = "Sample")
  setkey(metabolite_data, "Label")
  metabolite_data[]
}

bat_fig3 <- prepare_bat_normalized()
bat_fig4 <- format_bat_ik()
ewat_fig4 <- format_ewat()

use_data(bat_fig3, bat_fig4, ewat_fig4, overwrite = TRUE)
