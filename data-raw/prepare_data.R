library(org.Mm.eg.db)
library(data.table)
library(here)
library(magrittr)

counts <- fread(here::here("data-raw/featureCounts.csv.gz"))
usethis::use_data(counts, overwrite = TRUE)

#### Construct metadata
metadata <- readxl::read_xlsx(path = here::here("data-raw/metadata.xlsx"))
setDT(metadata)
usethis::use_data(metadata, overwrite = TRUE)

#### Prepare metabolites
formatData <- function(file){
  abundance <- fread(file, skip = 1)

  out <- melt(abundance, id.vars = "Label")
  setnames(out, c("Label", "Sample", "abundance"))
  out[, c("Genotype", "Timepoint"):=tstrsplit(Sample, "_")]
  out[, Genotype:=factor(Genotype, levels = c("WT", "KO"))]
  out[, Timepoint:=factor(Timepoint, levels = c("day", "night"))]
  setkey(out, "Label")
  out
}
metabolite_data1 <- formatData(here("data-raw/BAT25102019_both metaboanalyst auto scaled and log2.csv"))

metabolite_data1[, abundance:=trimws(abundance) %>% as.numeric]
usethis::use_data(metabolite_data1, overwrite = TRUE)

format_data <- function(file){
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
  metabolite_data
}

metabolite_data2 <- format_data(here("data-raw/BAT normalized.csv"))
usethis::use_data(metabolite_data2, overwrite = TRUE)
