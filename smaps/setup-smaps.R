# Setup ----
# Run from the paleo-islands project root in RStudio.
library(divextra)
library(diversitree)
library(ape)
library(phytools)

N.MAPS <- 2000             # Full design: maps PER clade, also defines pilot draws.
run.name <- "posterior-aic-saved-asr"
if (exists("mode") && mode == "pilot") run.name <- "posterior-aic-saved-asr-pilot"
maps.dir <- file.path("smaps/maps",run.name)
run.dir <- file.path("smaps/runs",run.name)
dir.create(maps.dir, recursive = TRUE, showWarnings = FALSE)
dir.create(run.dir, recursive = TRUE, showWarnings = FALSE)

# Data and fitted models: the same inputs used by MCMC ----
phy <- readRDS("mcmc/inputs/phy-271.rds")
geo.regions <- readRDS("mcmc/inputs/geo_regions-271.rds")
tip.states <- c(Afr = 1L, OP = 2L, Mada = 3L, Aus = 4L, Maur = 6L)
tip.states <- tip.states[as.character(geo.regions)]
names(tip.states) <- names(geo.regions)
mle.gr6 <- readRDS("mcmc/inputs/AfrSunken16_gr-6_r-5.rds")
mle.gr7 <- readRDS("mcmc/inputs/AfrSunken16_gr-7_r-1.rds")

state.labels <- c("A", "E", "M", "U", "S", "R", "A.E", "A.M", "A.U",
  "A.S", "A.R", "E.M", "E.U", "E.S", "E.R", "S.R")
colors <- setNames(rep("grey80", 16), state.labels)
colors[c("A", "S", "R", "A.S", "A.R", "S.R")] <- c(
  "#0072B2", "#E69F00", "#009E73", "#CC79A7", "#D55E00", "#56B4E9")

# Complete posterior vectors; discard 1,000 warmup sweeps from each chain ----
posterior.gr6 <- NULL
posterior.gr7 <- NULL
for (chain in 1:4) {
  x <- readRDS(paste0("mcmc/results/roihu-5000-20261023/gr6-chain-", chain, ".rds"))
  posterior.gr6 <- rbind(posterior.gr6, x$rates[1001:5000, ])
  x <- readRDS(paste0("mcmc/results/roihu-5000-20261023/gr7-chain-", chain, ".rds"))
  posterior.gr7 <- rbind(posterior.gr7, x$rates[1001:5000, ])
}

# AIC weights saved with the ASR reconstructions ----
model.weights <- readRDS("asr/results/model-weights.rds")
print(model.weights)
N.GR6 <- round(N.MAPS * model.weights$weight[model.weights$model == "gr6"])
N.GR7 <- N.MAPS - N.GR6

# Full-tree likelihood and model-specific parameter expansion ----
lik <- make.classe.td(phy, tip.states, k = 16, n.epoch = 2, strict = FALSE,
  sampling.f = c(.039, .027, .446, .008, 1e-6, .833, rep(1e-6, 10)),
  control = list(backend = "gslode"))
pars <- read_yaml_pars_td("mcmc/inputs/backbone-Afr-Sunken-16.yml")
groups <- read_par_groups_yaml("mcmc/inputs/groups-Afr-Sunken-16.yml")
pars.gr6 <- regroup_parameters(pars, groups[[6]])
pars.gr7 <- regroup_parameters(pars, groups[[7]])
lik.gr6 <- constrain(lik, formulae = make_constraints_sse_td(pars.gr6))
lik.gr7 <- constrain(lik, formulae = make_constraints_sse_td(pars.gr7))

# Saved clades and fixed ML-ASR root priors ----
sisyphini.gr6 <- readRDS("asr/results/sisyphini.gr6.rds")
sisyphini.gr7 <- readRDS("asr/results/sisyphini.gr7.rds")
nesovin.gr6 <- readRDS("asr/results/nesovin.gr6.rds")
nesovin.gr7 <- readRDS("asr/results/nesovin.gr7.rds")
sisyphini <- sisyphini.gr6$tree
nesovin <- nesovin.gr6$tree
tips.nesovin <- nesovin$tip.label[startsWith(nesovin$tip.label, "Nesovin")]
tips.nesosisyphus <- sisyphini$tip.label[startsWith(sisyphini$tip.label, "Nesosisyphus")]
stopifnot(length(tips.nesovin) == 1, length(tips.nesosisyphus) == 4)

# Identify the root by topology; reorder states explicitly.
root.sisyphini <- setdiff(sisyphini$edge[,1], sisyphini$edge[,2])
root.nesovin <- setdiff(nesovin$edge[,1], nesovin$edge[,2])
root.nesosisyphus.gr6 <- sisyphini.gr6$node.probs[state.labels, as.character(root.sisyphini)]
root.nesosisyphus.gr7 <- sisyphini.gr7$node.probs[state.labels, as.character(root.sisyphini)]
root.nesovin.gr6 <- nesovin.gr6$node.probs[state.labels, as.character(root.nesovin)]
root.nesovin.gr7 <- nesovin.gr7$node.probs[state.labels, as.character(root.nesovin)]

# Treat each subclade as a separate tree with its fixed ML-ASR vector as prior.
# The mapper updates that prior using the subclade likelihood.
# This two-stage plug-in analysis reuses descendant evidence; it is not an
# exact joint full-tree posterior. Root-prior uncertainty is not propagated.
saveRDS(list(nesovin.gr6 = root.nesovin.gr6, nesovin.gr7 = root.nesovin.gr7,
  nesosisyphus.gr6 = root.nesosisyphus.gr6, nesosisyphus.gr7 = root.nesosisyphus.gr7),
  file.path(run.dir, "fixed-root-priors.rds"))
write.csv(model.weights, file.path(run.dir, "model-weights.csv"), row.names = FALSE)

# Choose complete posterior rows before starting any maps ----
# Pair clades by model and posterior row for combined richness.
set.seed(10)
rows.nesosisyphus.gr6 <- sample.int(nrow(posterior.gr6), N.GR6, replace = TRUE)
rows.nesosisyphus.gr7 <- sample.int(nrow(posterior.gr7), N.GR7, replace = TRUE)
rows.nesovin.gr6 <- rows.nesosisyphus.gr6
rows.nesovin.gr7 <- rows.nesosisyphus.gr7
saveRDS(list(nesosisyphus.gr6 = rows.nesosisyphus.gr6,
  nesosisyphus.gr7 = rows.nesosisyphus.gr7,
  nesovin.gr6 = rows.nesovin.gr6, nesovin.gr7 = rows.nesovin.gr7),
  file.path(run.dir, "posterior-rows.rds"))

