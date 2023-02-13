# NAMPT-dependent NAD+ biosynthesis controls circadian metabolism in a tissue-specific manner:

Molecular clocks in the periphery coordinate tissue-specific daily biorhythms by integrating input from the hypothalamic master clock and intracellular metabolic signals. One such key metabolic signal is the cellular concentration of NAD+, which oscillates along with its biosynthetic enzyme, nicotinamide phosphoribosyltransferase (NAMPT). NAD+ levels feed back into the clock to influence rhythmicity of biological functions, yet whether this metabolic fine-tuning occurs ubiquitously across cell types and is a core clock feature is unknown. Here we show that NAMPT-dependent control over the molecular clock varies substantially between tissues. Brown adipose tissue (BAT) requires NAMPT to sustain the amplitude of the core clock, whereas rhythmicity in white adipose tissue (WAT) is only moderately dependent on NAD+ biosynthesis and the skeletal muscle clock is completely refractory to loss of NAMPT. In BAT and WAT, NAMPT differentially orchestrates oscillation of clock-controlled gene networks and the diurnality of metabolite levels. NAMPT coordinates the rhythmicity of TCA cycle intermediates in BAT, but not WAT, and loss of NAD+ abolishes these oscillations similarly to high fat diet (HFD)-induced circadian disruption. Adipose NAMPT depletion also improved the ability of animals to defend body temperature during cold stress but this is independent of time-of-day. Thus, our findings reveal that peripheral molecular clocks and metabolic biorhythms are shaped in a highly tissue-specific manner by NAMPT-dependent NAD+ synthesis.

## Folders and files

Having a standard folder and file structure is one step to being more reproducible.
Here is an explanation of some of the general contents of this project. The
project directory is generally structured with the following folders:

- `analysis/`: Contains R scripts for figure generation and differential gene
expression analysis.
- `data/`: Contains RData files of the processed data from data-raw.
- `data-raw/`: Contains raw data files as well a script for processing it.
- `out/`: Contains output files, sorted by type.
- `.git`: The Git version control history and related files. Used *only* by Git,
do not manually edit.

The base folder contains a few files:

- `.gitignore` tells [Git](https://git-scm.com/) to ignore certain files from
being tracked and prevents them from entering the version control history.
- `DESCRIPTION` is a standard file that includes metadata about your project, in
a machine readable format for others to obtain information on about your
project. It provides a description of what the project does and most importantly
what R packages your project uses on.
- The `project_171022.Rproj` file dictates that the directory is a RStudio
project. Open the project by opening this file.

All folders have their own README files inside with more specific details.

## Installing project software dependencies
The following packages are used in this project:
- cetcolor
- circular
- data.table
- edgeR
- eulerr
- export
- ggplot2
- here
- limorhyde
- magrittr
- openxlsx
- pheatmap
- readxl
- stringr
- SummarizedExperiment
- topGO
- UpSetR
- usethis
- viridis

## GEO Accession

The raw fastq files used in the project can be accessed at this GEO link:
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE221550
