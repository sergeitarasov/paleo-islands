# Rectangular ASR plots for the two focal clades ----
# Run from the paleo-islands project root. No inference is rerun.
library(ape)
library(divextra)

# Saved full-tree AIC-weighted reconstruction ----
x <- readRDS("asr/results/full.avg.rds")
state.labels <- rownames(x$node.probs)
colors <- c(
  A = "#23A4B8", E = "#E3C64A", M = "#8B5AA5", U = "#365C9A",
  S = "#E68613", R = "#C44E52",
  A.E = "#702D53", A.M = "#E58C72", A.U = "#D77DB0",
  A.S = "#A75D19", A.R = "#EFBB90", E.M = "#244E40",
  E.U = "#6F85CA", E.S = "#9A9C38", E.R = "#9DC8E5", S.R = "#76B99A"
)

# Extract matching nodes; retain their full-tree ASR probabilities ----
sisyphini.tips <- grep("^Sisyphus|^Sisyphini|^Nesosisyphus", x$tree$tip.label, value = TRUE)
nesovin.tips <- grep("^Nesovinsionia|^Nesovinsonia|^Tanzanolus", x$tree$tip.label, value = TRUE)
sisyphini <- extract.clade.asr(x$tree, getMRCA(x$tree, sisyphini.tips),
                              x$node.probs, x$tip.states, state.labels)
nesovin <- extract.clade.asr(x$tree, getMRCA(x$tree, nesovin.tips),
                            x$node.probs, x$tip.states, state.labels)
stopifnot(setequal(sisyphini$tree$tip.label, sisyphini.tips),
          setequal(nesovin$tree$tip.label, nesovin.tips))
clades <- list(sisyphini = sisyphini, "nesovinsonia-tanzanolus" = nesovin)

# Daughter-state ASR on the full tree ----
# Numeric branch time is measured BACKWARDS from the younger branch end.
# Its full length is therefore immediately after the parent's speciation.
# Use the same fitted runs and conditioning as asr-marginal.R.
mle.gr6 <- readRDS("ml-geosse/runs/rds/AfrSunken16_gr-6_r-1.rds")
mle.gr7 <- readRDS("ml-geosse/runs/rds/AfrSunken16_gr-7_r-1.rds")
weights <- x$analysis$model.weights
lik <- make.classe.td(x$tree, x$tip.states, k = 16, n.epoch = 2,
                     strict = FALSE,
                     sampling.f = c(.039, .027, .446, .008, 1e-6, .833, rep(1e-6, 10)),
                     control = list(backend = "gslode"))
fits <- list(gr6 = mle.gr6, gr7 = mle.gr7)
node.asr <- list(gr6 = readRDS("asr/results/full.gr6.rds"),
                 gr7 = readRDS("asr/results/full.gr7.rds"))
checks <- data.frame()
for (clade in names(clades)) {
  focal <- clades[[clade]]
  # Match subtree nodes to their original full-tree IDs, including tips.
  original <- integer(Ntip(focal$tree) + focal$tree$Nnode)
  original[seq_len(Ntip(focal$tree))] <- match(focal$tree$tip.label, x$tree$tip.label)
  original[focal$node.map$subtree.node] <- focal$node.map$original.node
  daughters <- focal$tree$edge[, 2]
  daughter.probs <- list()
  for (model in names(fits)) {
    probabilities <- matrix(NA_real_, 16, length(daughters),
                            dimnames = list(state.labels, as.character(daughters)))
    for (j in seq_along(daughters)) {
      node <- original[daughters[j]]
      length <- x$tree$edge.length[match(node, x$tree$edge[, 2])]
      probabilities[, j] <- asr.marginal.classe_branch(
        lik, fits[[model]]$par.full, node.id = node, relative.t.branch = length,
        root = diversitree::ROOT.GIVEN, root.p = c(1, rep(0, 15)))
      # At the younger end of internal branches, recover the saved node ASR.
      if (node > Ntip(x$tree)) {
        end <- asr.marginal.classe_branch(
          lik, fits[[model]]$par.full, node.id = node, relative.t.branch = 0,
          root = diversitree::ROOT.GIVEN, root.p = c(1, rep(0, 15)))
        error <- max(abs(end - node.asr[[model]]$node.probs[, as.character(node)]))
        checks <- rbind(checks, data.frame(clade, model, node, error))
      }
    }
    stopifnot(all(is.finite(probabilities)), min(probabilities) >= 0,
              max(abs(colSums(probabilities) - 1)) < 1e-8)
    daughter.probs[[model]] <- probabilities
  }
  focal$daughter.probs <- weights$weight[weights$model == "gr6"] * daughter.probs$gr6 +
    weights$weight[weights$model == "gr7"] * daughter.probs$gr7
  focal$daughter.probs.by.model <- daughter.probs
  focal$analysis <- x$analysis
  clades[[clade]] <- focal
}
stopifnot(max(checks$error) < 1e-5)
saveRDS(clades, "asr/results/focal-daughter-asr.rds")
write.csv(checks, "asr/results/focal-daughter-asr-checks.csv", row.names = FALSE)

# Display-only collapse of the sister clade to Nesosisyphus ----
# Keep its stem and crown age. No likelihoods or ASR are recalculated.
focal <- clades$sisyphini
other.tips <- grep("^Nesosisyphus", focal$tree$tip.label, invert = TRUE)
collapsed.node <- getMRCA(focal$tree, other.tips)
removed <- collapsed.node
frontier <- collapsed.node
while (length(frontier) > 0) {
  frontier <- focal$tree$edge[focal$tree$edge[, 1] %in% frontier, 2]
  removed <- c(removed, frontier)
}
keep.edges <- !focal$tree$edge[, 1] %in% removed
edges <- focal$tree$edge[keep.edges, , drop = FALSE]
kept.tips <- c(setdiff(seq_len(Ntip(focal$tree)), other.tips), collapsed.node)
kept.nodes <- sort(unique(edges[, 1]))
old.ids <- c(kept.tips, kept.nodes)
collapsed.tree <- list(edge = matrix(match(edges, old.ids), ncol = 2),
  edge.length = focal$tree$edge.length[keep.edges],
  tip.label = c(focal$tree$tip.label[setdiff(kept.tips, collapsed.node)],
                paste0("Other Sisyphini (", length(other.tips), " tips)")),
  Nnode = length(kept.nodes))
class(collapsed.tree) <- "phylo"
collapsed <- list(tree = reorder(collapsed.tree),
  node.probs = focal$node.probs[, as.character(kept.nodes), drop = FALSE],
  daughter.probs = focal$daughter.probs[, as.character(edges[, 2]), drop = FALSE],
  tip.states = focal$tip.states[collapsed.tree$tip.label],
  collapsed.tip = length(kept.tips),
  collapsed.probs = focal$node.probs[, as.character(collapsed.node)])
colnames(collapsed$node.probs) <- as.character(match(kept.nodes, old.ids))
colnames(collapsed$daughter.probs) <- as.character(match(edges[, 2], old.ids))
stopifnot(length(other.tips) == length(removed[removed <= Ntip(focal$tree)]))
clades[["nesosisyphus-collapsed"]] <- collapsed

# Save separate publication figures ----
dir.create("asr/plots", recursive = TRUE, showWarnings = FALSE)
for (clade in names(clades)) {
  focal <- clades[[clade]]
  tree <- ladderize(focal$tree)
  height <- max(node.depth.edgelength(tree))
  tip.states <- focal$tip.states[tree$tip.label]
  tip.col <- rep("#252525", Ntip(tree))
  tip.col[grep("^Nesovinsionia|^Nesovinsonia|^Nesosisyphus", tree$tip.label)] <- colors["R"]
  # Correct the genus spelling for display only; preserve saved identifiers.
  tree$tip.label <- sub("Nesovinsionia", "Nesovinsonia", tree$tip.label)
  tree$tip.label <- gsub("_", " ", tree$tip.label)
  tree$tip.label <- gsub(" +", " ", tree$tip.label)
  figure.height <- if (clade == "sisyphini") 9 else if (clade == "nesosisyphus-collapsed") 5.5 else 4.5
  for (format in c("pdf", "png")) {
    filename <- paste0("asr/plots/asr-", clade, "-daughters.", format)
    if (format == "pdf") pdf(filename, width = 12, height = figure.height, useDingbats = FALSE)
    if (format == "png") png(filename, width = 12, height = figure.height,
                             units = "in", res = 220)
    layout(matrix(c(1, 2), ncol = 1), heights = c(figure.height - 1.2, 1.2))
    par(mar = c(3.5, 1.5, 1, 1), family = "sans")
    display.labels <- tree$tip.label
    if (!is.null(focal$collapsed.tip)) display.labels[focal$collapsed.tip] <- ""
    display.tree <- tree
    display.tree$tip.label <- display.labels
    plot(display.tree, type = "phylogram", direction = "rightwards",
         edge.color = "#555555", edge.width = 1.2,
         font = 3, cex = .85, tip.color = tip.col,
         label.offset = height * .025, x.lim = c(-.04, 1.95) * height,
         y.lim = c(.5, Ntip(tree) + .5), no.margin = FALSE)
    nodelabels(node = as.integer(colnames(focal$node.probs)),
               pie = t(focal$node.probs[state.labels, , drop = FALSE]),
               piecol = colors[state.labels], cex = .45)
    # Pies show exact branch-start ASR at each horizontal/vertical corner.
    coordinates <- get("last_plot.phylo", envir = ape::.PlotPhyloEnv)
    parents <- tree$edge[, 1]
    daughters <- tree$edge[, 2]
    pie.x <- coordinates$xx[parents]
    pie.y <- coordinates$yy[daughters]
    for (j in seq_along(daughters)) {
      ape:::floating.pie.asp(pie.x[j], pie.y[j],
        focal$daughter.probs[, as.character(daughters[j])],
        radius = height * .013, col = colors[state.labels])
    }
    observed <- which(!is.na(tip.states))
    tiplabels(tip = observed, pch = 19,
              col = colors[state.labels[tip.states[observed]]], cex = 1.1)
    if (!is.null(focal$collapsed.tip)) {
      i <- focal$collapsed.tip
      # Triangle spans the collapsed crown age to the present.
      polygon(c(coordinates$xx[i], height, height),
              c(coordinates$yy[i], coordinates$yy[i] - .3, coordinates$yy[i] + .3),
              col = "#EEEEEE", border = "#777777")
      ape:::floating.pie.asp(coordinates$xx[i], coordinates$yy[i],
                            focal$collapsed.probs, radius = height * .018,
                            col = colors[state.labels])
      text(height * 1.025, coordinates$yy[i], tree$tip.label[i],
           adj = 0, cex = .85)
    }
    axisPhylo(backward = TRUE, cex.axis = .85)
    mtext("Age (Ma)", side = 1, line = 2.3, cex = .9)
    par(mar = c(0, 0, 0, 0))
    plot.new()
    legend("center", legend = gsub(".", "-", state.labels, fixed = TRUE), fill = colors[state.labels],
           border = NA, ncol = 8, bty = "n", cex = .9)
    mtext("Large pies: before speciation; small pies: daughter states immediately after speciation",
          side = 1, line = -1, cex = .7)
    dev.off()
  }
}
