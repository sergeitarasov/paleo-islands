# Open the paleo-islands RStudio project. Run this script from its root.
library(readxl)

# Modern-island calibration data ----
islands <- as.data.frame(read_excel("data/richness/island-richness-area.xlsx"))

# Full, unsmoothed richness distribution at exactly 9 Ma ----
richness <- readRDS("smaps/summaries/posterior-aic-saved-asr/richness-through-time.rds")
i.nesovin <- which(abs(richness$nesovin$times - 9) < 1e-8)
i.nesosisyphus <- which(abs(richness$nesosisyphus$times - 9) < 1e-8)
i.combined <- which(abs(richness$combined$summary$time - 9) < 1e-8)
stopifnot(length(i.nesovin) == 1, length(i.nesosisyphus) == 1, length(i.combined) == 1)

counts <- data.frame(
  map = seq_len(ncol(richness$combined$counts)),
  age_ma = 9,
  Nesovinsonia = as.integer(richness$nesovin$counts[i.nesovin, ]),
  Nesosisyphus = as.integer(richness$nesosisyphus$counts[i.nesosisyphus, ]),
  Combined = as.integer(richness$combined$counts[i.combined, ])
)
# The upstream script pairs maps by model and posterior row before summing.
stopifnot(all(counts$Combined == counts$Nesovinsonia + counts$Nesosisyphus))

# Summarize, without truncating the tails or smoothing ----
summary <- data.frame()
for (clade in c("Nesovinsonia", "Nesosisyphus", "Combined")) {
  n <- counts[[clade]]
  summary <- rbind(summary, data.frame(
    clade = clade, age_ma = 9, n_maps = length(n), mean = mean(n),
    median = median(n), lower = quantile(n, .025, type = 7),
    upper = quantile(n, .975, type = 7), minimum = min(n), maximum = max(n)
  ))
}
frequency <- as.data.frame(table(counts$Combined))
names(frequency) <- c("n_species", "n_maps")
frequency$n_species <- as.integer(as.character(frequency$n_species))
frequency$weight <- frequency$n_maps / nrow(counts)

# Derived inputs; original workbook and map summaries stay in place ----
dir.create("area-estimation/results/inputs", recursive = TRUE, showWarnings = FALSE)
write.csv(islands, "area-estimation/results/inputs/island-calibration.csv", row.names = FALSE)
write.csv(counts, "area-estimation/results/inputs/richness-9ma-draws.csv", row.names = FALSE)
write.csv(summary, "area-estimation/results/inputs/richness-9ma-summary.csv", row.names = FALSE)
write.csv(frequency, "area-estimation/results/inputs/richness-9ma-frequency.csv", row.names = FALSE)

# Relative source paths and checksums identify the input versions.
sources <- c("data/richness/island-richness-area.xlsx",
  "smaps/summaries/posterior-aic-saved-asr/richness-through-time.rds",
  "smaps/lineages-through-time-simple.R")
provenance <- data.frame(source = sources, md5 = unname(tools::md5sum(sources)))
write.csv(provenance, "area-estimation/results/inputs/source-files.csv", row.names = FALSE)
summary
