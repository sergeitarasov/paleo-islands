# Setup and input data ----
# Run from the divextra-sunken project root.
library(divextra)
library(ape)
library(phytools)

results.dir <- "asr/results/"
dir.create(results.dir, recursive = TRUE, showWarnings = FALSE)
phy <- readRDS("data/phy-271.rds")
geo_regions_abc <- readRDS("data/geo_regions-271.rds")
geo_regions <- setNames(
  unname(c(Afr = 1L, OP = 2L, Mada = 3L, Aus = 4L, Maur = 6L)[
    as.character(geo_regions_abc)]), names(geo_regions_abc)
)

# State labels and plot colours ----
state.labels <- c("A", "E", "M", "U", "S", "R", "A.E", "A.M", "A.U",
                  "A.S", "A.R", "E.M", "E.U", "E.S", "E.R", "S.R")
colors <- setNames(rep("grey80", length(state.labels)), state.labels)
colors[c("A", "S", "R", "A.S", "A.R", "S.R")] <- c(
  "#0072B2", "#E69F00", "#009E73", "#CC79A7", "#D55E00", "#56B4E9"
)

# Fitted models and AIC weights ----
mle.gr6 <- readRDS("runs/rds/AfrSunken16_gr-6_r-1.rds")
mle.gr7 <- readRDS("runs/rds/AfrSunken16_gr-7_r-1.rds")
model.weights <- data.frame(
  model = c("gr6", "gr7"),
  logLik = c(mle.gr6$lnLik, mle.gr7$lnLik),
  npar = c(mle.gr6$NPAR, mle.gr7$NPAR),
  AIC = c(mle.gr6$AIC, mle.gr7$AIC)
)
stopifnot(all(is.finite(model.weights$AIC)),
          all(abs(model.weights$AIC -
                    (2 * model.weights$npar - 2 * model.weights$logLik)) < 1e-6))
model.weights$delta.AIC <- model.weights$AIC - min(model.weights$AIC)
model.weights$weight <- exp(-model.weights$delta.AIC / 2)
model.weights$weight <- model.weights$weight / sum(model.weights$weight)
print(model.weights)
# Weights are normalized across these TWO fitted models only.
saveRDS(model.weights, file.path(results.dir, "model-weights.rds"))

# Whole-tree ASR: gr6 ----
lik.td <- make.classe.td(
  phy, geo_regions, k = 16, n.epoch = 2, strict = FALSE,
  sampling.f = c(.039, .027, .446, .008, 1e-6, .833, rep(1e-6, 10)),
  control = list(backend = "gslode")
)
asr.gr6 <- asr.marginal.classe(
  lik.td, mle.gr6$par.full,
  root = diversitree::ROOT.GIVEN, root.p = c(1, rep(0, 15))
)
dimnames(asr.gr6) <- list(state.labels, as.character(Ntip(phy) + seq_len(phy$Nnode)))

# Whole-tree ASR: gr7 ----
asr.gr7 <- asr.marginal.classe(
  lik.td, mle.gr7$par.full,
  root = diversitree::ROOT.GIVEN, root.p = c(1, rep(0, 15))
)
dimnames(asr.gr7) <- dimnames(asr.gr6)

# Whole-tree AIC-weighted average ----
# Average reconstructed probabilities, NOT model parameters or sampled states.
stopifnot(identical(dimnames(asr.gr6), dimnames(asr.gr7)))
asr.avg <- model.weights$weight[model.weights$model == "gr6"] * asr.gr6 +
  model.weights$weight[model.weights$model == "gr7"] * asr.gr7
stopifnot(all(is.finite(asr.avg)), all(asr.avg >= 0),
          max(abs(colSums(asr.avg) - 1)) < 1e-8)
saveRDS(asr.gr6, file.path(results.dir, "asr.gr6.rds"))
saveRDS(asr.gr7, file.path(results.dir, "asr.gr7.rds"))
saveRDS(asr.avg, file.path(results.dir, "asr.avg.rds"))

# Reusable ASR plot and export functions ----
plot.asr.result <- function(x, title = "") {
  phytools::plotTree(
    x$tree, fsize = if (Ntip(x$tree) > 100) .65 else .7,
    ftype = "i", lwd = 2, offset = .6, mar = c(2, 1, 3, 1)
  )
  ape::nodelabels(
    node = as.integer(colnames(x$node.probs)),
    pie = t(x$node.probs[state.labels, , drop = FALSE]),
    piecol = colors[state.labels],
    cex = if (Ntip(x$tree) > 100) .15 else .35
  )
  ape::tiplabels(
    pie = outer(as.integer(x$tip.states[x$tree$tip.label]),
                seq_along(state.labels), "==") * 1,
    piecol = colors[state.labels], cex = if (Ntip(x$tree) > 100) .15 else .35
  )
  graphics::title(main = title, cex.main = 1)
  ape::axisPhylo()
}

save.asr.result <- function(x, name, title, model) {
  # Each RDS contains the tree, node probabilities, tip states and node mapping.
  # Also retain state colours and averaging provenance for later reuse.
  x$analysis <- list(model = model, model.weights = model.weights,
                     state.labels = state.labels, colors = colors)
  saveRDS(x, file.path(results.dir, paste0(name, ".rds")))
  height <- max(8, Ntip(x$tree) * .16 + 2)
  grDevices::pdf(file.path(results.dir, paste0(name, ".pdf")),
                width = 16, height = height, useDingbats = FALSE)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::layout(matrix(1:2, ncol = 1), heights = c(height - .65, .65))
  plot.asr.result(x, title)
  graphics::par(mar = c(0, 0, 0, 0))
  graphics::plot.new()
  graphics::legend(
    "center", legend = c("A", "S", "R", "A.S", "A.R", "S.R", "Other"),
    fill = c(colors[c("A", "S", "R", "A.S", "A.R", "S.R")], "grey80"),
    cex = .85, bty = "n", horiz = TRUE
  )
  invisible(x)
}

# Export whole-tree results ----
root.node <- setdiff(phy$edge[, 1], phy$edge[, 2])
full.gr6 <- extract.clade.asr(phy, root.node, asr.gr6, geo_regions, state.labels)
full.gr7 <- extract.clade.asr(phy, root.node, asr.gr7, geo_regions, state.labels)
full.avg <- extract.clade.asr(phy, root.node, asr.avg, geo_regions, state.labels)
average.title <- sprintf("AIC-weighted average (gr6 %.3f; gr7 %.3f)",
                         model.weights$weight[1], model.weights$weight[2])
save.asr.result(full.gr6, "full.gr6", "Whole-tree ASR - gr6", "gr6")
save.asr.result(full.gr7, "full.gr7", "Whole-tree ASR - gr7", "gr7")
save.asr.result(full.avg, "full.avg", paste("Whole-tree ASR -", average.title), "AIC average")

# Sisyphini: extract and export ----
sisyphini.gr6 <- extract.clade.asr(phy, 495, asr.gr6, geo_regions, state.labels)
sisyphini.gr7 <- extract.clade.asr(phy, 495, asr.gr7, geo_regions, state.labels)
sisyphini.avg <- extract.clade.asr(phy, 495, asr.avg, geo_regions, state.labels)
save.asr.result(sisyphini.gr6, "sisyphini.gr6", "Sisyphini ASR - gr6", "gr6")
save.asr.result(sisyphini.gr7, "sisyphini.gr7", "Sisyphini ASR - gr7", "gr7")
save.asr.result(sisyphini.avg, "sisyphini.avg", paste("Sisyphini ASR -", average.title), "AIC average")

# Nesovinsonia: extract and export ----
nesovin.gr6 <- extract.clade.asr(phy, 489, asr.gr6, geo_regions, state.labels)
nesovin.gr7 <- extract.clade.asr(phy, 489, asr.gr7, geo_regions, state.labels)
nesovin.avg <- extract.clade.asr(phy, 489, asr.avg, geo_regions, state.labels)
save.asr.result(nesovin.gr6, "nesovin.gr6", "Nesovinsonia ASR - gr6", "gr6")
save.asr.result(nesovin.gr7, "nesovin.gr7", "Nesovinsonia ASR - gr7", "gr7")
save.asr.result(nesovin.avg, "nesovin.avg", paste("Nesovinsonia ASR -", average.title), "AIC average")

# Interactive plotting ----
# Example: plot.asr.result(sisyphini.avg, "Sisyphini - AIC-weighted ASR")
# Saved plots and data are in asr/results; no simulations are run by this script.
