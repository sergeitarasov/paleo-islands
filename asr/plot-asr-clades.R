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
  figure.height <- if (clade == "sisyphini") 9 else 4.5
  for (format in c("pdf", "png")) {
    filename <- paste0("asr/plots/asr-", clade, ".", format)
    if (format == "pdf") pdf(filename, width = 12, height = figure.height, useDingbats = FALSE)
    if (format == "png") png(filename, width = 12, height = figure.height,
                             units = "in", res = 220)
    layout(matrix(c(1, 2), ncol = 1), heights = c(figure.height - 1.2, 1.2))
    par(mar = c(3.5, 1.5, 1, 1), family = "sans")
    plot(tree, type = "phylogram", direction = "rightwards",
         edge.color = "#555555", edge.width = 1.2,
         font = 3, cex = .85, tip.color = tip.col,
         label.offset = height * .025, x.lim = c(-.04, 1.95) * height,
         y.lim = c(.5, Ntip(tree) + .5), no.margin = FALSE)
    nodelabels(node = as.integer(colnames(focal$node.probs)),
               pie = t(focal$node.probs[state.labels, , drop = FALSE]),
               piecol = colors[state.labels], cex = .55)
    tiplabels(pch = 19, col = colors[state.labels[tip.states]], cex = 1.1)
    axisPhylo(backward = TRUE, cex.axis = .85)
    mtext("Age (Ma)", side = 1, line = 2.3, cex = .9)
    par(mar = c(0, 0, 0, 0))
    plot.new()
    legend("center", legend = state.labels, fill = colors[state.labels],
           border = NA, ncol = 8, bty = "n", cex = .9)
    dev.off()
  }
}
