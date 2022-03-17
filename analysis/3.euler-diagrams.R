library(eulerr)
library(here)
### Figure 2
# D
png(filename = here("out/euler-diagrams/2D.png"), width = 15, height = 15, units = "cm", res = 300)
plot(euler(combinations = c("eWAT" = 1295,
                            "BAT" = 2231,
                            "eWAT&BAT" = 590
)),
legend = FALSE,
labels = FALSE,
edges = FALSE,
fills = c("eWAT" = "#A0A0A4",
          "BAT" = "#606060"
))
dev.off()

# E
png(filename = here("out/euler-diagrams/2E.png"), width = 15, height = 15, units = "cm", res = 300)
plot(euler(combinations = c("WT" = 2500,
                            "BAT" = 259,
                            "BAT&WT" = 321
)),
legend = FALSE,
labels = FALSE,
edges = FALSE,
fills = c("WT" = "#606060",
          "BAT" = "#A00000"
))
dev.off()

# F
png(filename = here("out/euler-diagrams/2F.png"), width = 15, height = 15, units = "cm", res = 300)
plot(euler(combinations = c("WT" = 1460,
                            "eWAT" = 954,
                            "eWAT&WT" = 425
)),
legend = FALSE,
labels = FALSE,
edges = FALSE,
fills = c("WT" = "#606060",
          "eWAT" = "#94641F"
))
dev.off()


### Figure 3
# A
png(filename = here("out/euler-diagrams/3A.png"), width = 15, height = 15, units = "cm", res = 300)
plot(euler(combinations = c("genotype" = 199,
                            "time" = 119,
                            "interaction" = 4,
                            "genotype&time" = 138,
                            "interaction&genotype" = 24,
                            "interaction&time" = 32,
                            "interaction&genotype&time" = 33
), control = list(extraopt_threshold = 0), shape = "ellipse"),
legend = FALSE,
labels = FALSE,
edges = FALSE,
fills = c("genotype" = "#A00000",
          "time" = "#60856C",
          "interaction" = "#FFA040"))
dev.off()

# F
png(filename = here("out/euler-diagrams/3F.png"), width = 15, height = 15, units = "cm", res = 300)
plot(euler(combinations = c("genotype" = 140,
                            "time" = 18,
                            "genotype&time" = 4
)),
legend = FALSE,
labels = FALSE,
edges = FALSE,
fills = c("genotype" = "#94641F",
          "time" = "#50856C"
))
dev.off()

### Figure 4
# A
png(filename = here("out/euler-diagrams/4A.png"), width = 15, height = 15, units = "cm", res = 300)
plot(euler(combinations = c("genotype" = 229,
                            "time" = 72,
                            "HFD" = 114,
                            "genotype&time" = 105,
                            "HFD&genotype" = 95,
                            "HFD&time" = 31,
                            "HFD&genotype&time" = 43
), control = list(extraopt_threshold = 0), shape = "ellipse"),
legend = FALSE,
labels = FALSE,
edges = FALSE,
fills = c("genotype" = "#A00000",
          "time" = "#60856C",
          "HFD" = "#0F99B2"))
dev.off()

# E
png(filename = here("out/euler-diagrams/4E.png"), width = 15, height = 15, units = "cm", res = 300)
plot(euler(combinations = c("genotype" = 183,
                            "time" = 167,
                            "interaction" = 1,
                            "genotype&time" = 62,
                            "interaction&genotype" = 6,
                            "interaction&time" = 3,
                            "interaction&genotype&time" = 11
), control = list(extraopt_threshold = 0), shape = "ellipse"),
legend = FALSE,
labels = FALSE,
edges = FALSE,
fills = c("genotype" = "#A00000",
          "time" = "#60856C",
          "interaction" = "#FFA040"))
dev.off()

# F
png(filename = here("out/euler-diagrams/4F.png"), width = 15, height = 15, units = "cm", res = 300)
plot(euler(c("a" = 341)), legend = FALSE, labels = FALSE, edges = FALSE, fills = "#A00000")
dev.off()
