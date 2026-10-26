# Reference ML example. Active posterior workflow: bash smaps/run-smaps-local.sh
# Note ----

## Use make-smaps-parallel.R

# Setup and input data ----
# Run from the divextra-sunken project root.
devtools::load_all("/Users/taravser/GitHub/divextra")
library(ape)
library(phytools)

# Args ----
results.dir <- "asr/results/smaps"

# Read ----
model.weights <- readRDS("asr/results/model-weights.rds")
mle.gr6 <- readRDS("runs/rds/AfrSunken16_gr-6_r-1.rds")
mle.gr7 <- readRDS("runs/rds/AfrSunken16_gr-7_r-1.rds")

sisyphini.gr6 <- readRDS("asr/results/sisyphini.gr6.rds")
sisyphini.gr7 <- readRDS("asr/results/sisyphini.gr7.rds")

nesovin.gr6 <- readRDS("asr/results/nesovin.gr6.rds")
nesovin.gr7 <- readRDS("asr/results/nesovin.gr7.rds")




# State labels and plot colours ----
state.labels <- c("A", "E", "M", "U", "S", "R", "A.E", "A.M", "A.U",
                  "A.S", "A.R", "E.M", "E.U", "E.S", "E.R", "S.R")
colors <- setNames(rep("grey80", length(state.labels)), state.labels)
colors[c("A", "S", "R", "A.S", "A.R", "S.R")] <- c(
  "#0072B2", "#E69F00", "#009E73", "#CC79A7", "#D55E00", "#56B4E9"
)
#



# Model weights ----
# simulate proportional to weights
print(model.weights)
N.MAPS <- 2000
N.GR6 <- round(N.MAPS * model.weights$weight[
  model.weights$model == "gr6"
])
N.GR7 <- N.MAPS - N.GR6
N.GR7
N.GR6


# Sisyphini ----

x <- sisyphini.gr6
plot(x$tree)
NSIM <- N.GR6
state.labels <- rownames(x$node.probs)
tips <- x$tree$tip.label[startsWith(x$tree$tip.label, "Nesosisyphus")]
root.p = x$node.probs[, as.character(length(x$tree$tip.label) + 1)]

timing <- system.time(
  smaps.nesosisyphus.gr6 <- make.simmap.classe.td(
    tree = x$tree,
    pars = mle.gr6$par.full,
    tip.weights = x$tip.states,
    k = 16,
    state.labels = state.labels,
    backbone.method = "tree-conditioned",
    root.prior = root.p,
    side.lineage.mode = "extinct_clades",
    side.lineage.tips = tips,
    nsim = NSIM,
    seed = NULL
  )
)

timing
timing["elapsed"]


x <- sisyphini.gr7
plot(x$tree)
NSIM <- N.GR7
state.labels <- rownames(x$node.probs)
tips <- x$tree$tip.label[startsWith(x$tree$tip.label, "Nesosisyphus")]
root.p = x$node.probs[, as.character(length(x$tree$tip.label) + 1)]

timing <- system.time(
  smaps.nesosisyphus.gr7 <- make.simmap.classe.td(
    tree = x$tree,
    pars = mle.gr7$par.full,
    tip.weights = x$tip.states,
    k = 16,
    state.labels = state.labels,
    backbone.method = "tree-conditioned",
    root.prior = root.p,
    side.lineage.mode = "extinct_clades",
    side.lineage.tips = tips,
    nsim = NSIM,
    seed = NULL
  )
)

timing
timing["elapsed"]

## combine and save ----
smaps.sisyphini <- c(smaps.nesosisyphus.gr6, smaps.nesosisyphus.gr7)
saveRDS(smaps.sisyphini, file.path(results.dir, "smaps-sisyphini-full.rds"))

## plot ----
s <- smaps.sisyphini[[15]]
phytools::plotSimmap(s$tree, colors = colors, split.vertical = TRUE, ftype = "i", fsize = 0.5)
legend(
  "topleft",
  legend = c("A", "S", "R", "A.S", "A.R", "S.R", "Other"),
  fill = c(colors[c("A", "S", "R", "A.S", "A.R", "S.R")], "grey80"),
  cex = 0.7, bty = "n"
)
table(s$tips$fate)




# Nesovinsonia ----

x <- nesovin.gr6
plot(x$tree)
NSIM <- N.GR6
state.labels <- rownames(x$node.probs)
tips <- x$tree$tip.label[startsWith(x$tree$tip.label, "Nesovin")]
root.p = x$node.probs[, as.character(length(x$tree$tip.label) + 1)]

timing <- system.time(
  smaps.nesovin.gr6 <- make.simmap.classe.td(
    tree = x$tree,
    pars = mle.gr6$par.full,
    tip.weights = x$tip.states,
    k = 16,
    state.labels = state.labels,
    backbone.method = "tree-conditioned",
    root.prior = root.p,
    side.lineage.mode = "extinct_clades",
    side.lineage.tips = tips,
    nsim = NSIM,
    seed = NULL
  )
)

timing
timing["elapsed"]


x <- nesovin.gr7
plot(x$tree)
NSIM <- N.GR7
state.labels <- rownames(x$node.probs)
tips <- x$tree$tip.label[startsWith(x$tree$tip.label, "Nesovin")]
root.p = x$node.probs[, as.character(length(x$tree$tip.label) + 1)]

timing <- system.time(
  smaps.nesovin.gr7 <- make.simmap.classe.td(
    tree = x$tree,
    pars = mle.gr7$par.full,
    tip.weights = x$tip.states,
    k = 16,
    state.labels = state.labels,
    backbone.method = "tree-conditioned",
    root.prior = root.p,
    side.lineage.mode = "extinct_clades",
    side.lineage.tips = tips,
    nsim = NSIM,
    seed = NULL
  )
)

timing
timing["elapsed"]

## combine and save ----
smaps.nesovin <- c(smaps.nesovin.gr6, smaps.nesovin.gr7)
saveRDS(smaps.nesovin, file.path(results.dir, "smaps-nesovin-full.rds"))

## plot ----
s <- smaps.nesovin[[3]]
phytools::plotSimmap(s$tree, colors = colors, split.vertical = TRUE, ftype = "i", fsize = 0.5)
legend(
  "topleft",
  legend = c("A", "S", "R", "A.S", "A.R", "S.R", "Other"),
  fill = c(colors[c("A", "S", "R", "A.S", "A.R", "S.R")], "grey80"),
  cex = 0.7, bty = "n"
)
axisPhylo()
table(s$tips$fate)


# tt <- lineages.through.time.simmap(
#   smaps.nesovin[[3]],
#   times = seq(0, 70, by = 1),
#   states = c("S", "A.S", "S.R"),
#   tips = "Nesovinsionia_vinsoni _STL4"
# )
#
#
# tt$summary



