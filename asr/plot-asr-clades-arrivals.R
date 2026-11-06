# Focal ASR trees on one common age axis ----
# Run from the paleo-islands project root after plot-asr-clades-daughters.R.
# Uses saved full-tree ASR; no inference is rerun and the trees are not joined.
library(ape)
clades <- readRDS("asr/results/focal-daughter-asr.rds")
state.labels <- rownames(clades$sisyphini$node.probs)
colors <- c(
  A = "#23A4B8", E = "#E3C64A", M = "#8B5AA5", U = "#365C9A",
  S = "#E68613", R = "#C44E52",
  A.E = "#702D53", A.M = "#E58C72", A.U = "#D77DB0",
  A.S = "#A75D19", A.R = "#EFBB90", E.M = "#244E40",
  E.U = "#6F85CA", E.S = "#9A9C38", E.R = "#9DC8E5", S.R = "#76B99A"
)

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


# Saved arrivals from posterior-parameter maps, AIC-weighted across gr6/gr7.
arrival.ages <- read.csv("smaps/summaries/posterior-aic-saved-asr/first-arrival-ages.csv")
arrival.summary <- read.csv("smaps/summaries/posterior-aic-saved-asr/first-arrival-summary.csv")
# Bounds reflect the clade roots and the model's 9 Ma Mauritius epoch.
upper.ages <- c(max(node.depth.edgelength(clades[["nesovinsonia-tanzanolus"]]$tree)), 9,
                max(node.depth.edgelength(clades$sisyphini$tree)), 9)

selected <- c("nesovinsonia-tanzanolus", "nesosisyphus-collapsed")
for (format in c("pdf", "png")) {
  filename <- paste0("asr/plots/asr-focal-clades-common-timeline-arrivals.", format)
  if (format == "pdf") pdf(filename, width = 13, height = 10.5, useDingbats = FALSE)
  if (format == "png") png(filename, width = 13, height = 10.5, units = "in", res = 220)
  layout(matrix(c(1, 2), ncol = 1), heights = c(9, 1.5))
  par(mar = c(3.5, 1.2, 1, 1), family = "sans")
  plot.new()
  plot.window(xlim = c(65, -41), ylim = c(-6.5, 10.5), xaxs = "i", yaxs = "i")
  # Light guides and a common present-day line.
  for (age in seq(0, 60, by = 10)) {
    segments(age, -6.1, age, 10.2, col = "#E4E7EA", lty = 3, lwd = .7)
  }
  axis(1, at = seq(0, 60, by = 10), labels = seq(0, 60, by = 10), cex.axis = .85)
  mtext("Age (Ma)", side = 1, line = 2.3, at = 30, cex = .95)

  for (clade in selected) {
    focal <- clades[[clade]]
    tree <- reorder(ladderize(focal$tree), "cladewise")
    n <- Ntip(tree)
    node.age <- max(node.depth.edgelength(tree)) - node.depth.edgelength(tree)
    # Set tip positions; internal nodes sit midway between their two children.
    y <- numeric(n + tree$Nnode)
    tip.order <- tree$edge[tree$edge[, 2] <= n, 2]
    y[tip.order] <- rev(seq_len(n))
    if (clade == "nesovinsonia-tanzanolus") y[seq_len(n)] <- y[seq_len(n)] + 7
    postorder <- reorder(tree, "postorder")
    for (parent in unique(postorder$edge[, 1])) {
      children <- tree$edge[tree$edge[, 1] == parent, 2]
      y[parent] <- mean(y[children])
    }
    for (parent in unique(tree$edge[, 1])) {
      children <- tree$edge[tree$edge[, 1] == parent, 2]
      segments(node.age[parent], min(y[children]), node.age[parent], max(y[children]),
               col = "#555555", lwd = 1.1)
    }
    for (j in seq_len(nrow(tree$edge))) {
      parent <- tree$edge[j, 1]
      child <- tree$edge[j, 2]
      segments(node.age[parent], y[child], node.age[child], y[child],
               col = "#555555", lwd = 1.1)
    }
    if (!is.null(focal$collapsed.tip)) {
      i <- focal$collapsed.tip
      polygon(c(node.age[i], 0, 0), c(y[i], y[i] - .25, y[i] + .25),
              col = "#EEEEEE", border = "#777777")
      ape:::floating.pie.asp(node.age[i], y[i], focal$collapsed.probs,
                            radius = .7, col = colors[state.labels])
    }
    for (node in as.integer(colnames(focal$node.probs))) {
      ape:::floating.pie.asp(node.age[node], y[node], focal$node.probs[, as.character(node)],
                            radius = .7, col = colors[state.labels])
    }
    for (j in seq_len(nrow(tree$edge))) {
      parent <- tree$edge[j, 1]
      child <- tree$edge[j, 2]
      ape:::floating.pie.asp(node.age[parent], y[child],
        focal$daughter.probs[, as.character(child)], radius = .46,
        col = colors[state.labels])
    }
    for (i in seq_len(n)) {
      label <- gsub("_", " ", tree$tip.label[i])
      label <- gsub(" +", " ", label)
      label <- sub("Nesovinsionia", "Nesovinsonia", label)
      colour <- if (grepl("^Nesovinsonia|^Nesosisyphus", label)) colors["R"] else "#252525"
      text(-1.2, y[i], label, adj = 0, col = colour, font = 3, cex = .85)
      if (is.null(focal$collapsed.tip) || i != focal$collapsed.tip) {
        state <- focal$tip.states[tree$tip.label[i]]
        points(0, y[i], pch = 19, col = colors[state.labels[state]], cex = .85)
      }
    }
  }
  # Arrival-time densities below the trees, on exactly the same age scale.
  text(63, -.6, "First recorded arrivals from posterior-parameter maps", adj = 0,
       cex = .9, font = 2, col = "#444444")
  for (i in seq_len(nrow(arrival.summary))) {
    row <- arrival.summary[i, ]
    ages <- arrival.ages$age[arrival.ages$clade == row$clade & arrival.ages$route == row$route]
    ages <- ages[is.finite(ages)]
    stopifnot(length(ages) == row$n.arrivals, all(ages >= 0), all(ages <= upper.ages[i]))
    curve <- density(c(ages, -ages, 2 * upper.ages[i] - ages), bw = bw.nrd0(ages),
                     from = 0, to = upper.ages[i], n = 2048)
    # Each density has the same plotted height; height is not arrival probability.
    curve$y <- .7 * curve$y / max(curve$y)
    baseline <- -1.7 - (i - 1) * 1.35
    region <- if (grepl("Sunken", row$route)) "S" else "R"
    colour <- colors[region]
    polygon(c(0, curve$x, upper.ages[i]), c(baseline, baseline + curve$y, baseline),
            col = adjustcolor(colour, alpha.f = .25), border = NA)
    lines(curve$x, baseline + curve$y, col = colour, lwd = 1.1)
    segments(row$lower, baseline - .10, row$upper, baseline - .10, col = colour, lwd = 2.5)
    points(row$median, baseline - .10, pch = 19, col = colour, cex = .75)
    destination <- if (region == "S") "Sunken Islands" else "Mauritius"
    text(-1.2, baseline + .35, paste(row$clade, destination, sep = " - "),
         adj = 0, cex = .8, col = colour)
    text(-1.2, baseline, sprintf("%.2f Ma [%.2f-%.2f]", row$median, row$lower, row$upper),
         adj = 0, cex = .75, col = "#555555")
  }
  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend("center", legend = gsub(".", "-", state.labels, fixed = TRUE),
         fill = colors[state.labels], border = NA, ncol = 8, bty = "n", cex = .9)
  mtext("Large pies: before speciation; small pies: daughter states immediately after speciation",
        side = 1, line = -3.2, cex = .7)
  mtext("Arrivals: dot = median; bar = 95% interval; densities conditional on a matching transition",
        side = 1, line = -1.9, cex = .7)
  dev.off()
}
