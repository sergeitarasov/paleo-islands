# Whole-tree ancestral ranges: semicircular fan ----
# Run from the paleo-islands project root. Uses saved ASR, not new inference.
library(ape)

# Reconstruction ----
# To plot one model, change full.avg.rds to full.gr6.rds or full.gr7.rds.
x <- readRDS("asr/results/full.avg.rds")
output.name <- "full-asr-fan-average"
phy <- ladderize(x$tree)
state.labels <- c("A", "E", "M", "U", "S", "R", "A.E", "A.M", "A.U",
                  "A.S", "A.R", "E.M", "E.U", "E.S", "E.R", "S.R")
node.probs <- x$node.probs[state.labels, , drop = FALSE]
tip.states <- x$tip.states[phy$tip.label]
stopifnot(!anyNA(tip.states), max(abs(colSums(node.probs) - 1)) < 1e-8)

# Colours: all 16 geographic states are shown ----
# Single regions: Africa turquoise-blue, Eurasia-Oriental golden yellow, Madagascar purple,
# Australia deep blue, Sunken Islands orange, Mauritius brick red.
# Joint ranges have their own colours; none are pooled as Other.
colors <- c(
  A = "#23A4B8", E = "#E3C64A", M = "#8B5AA5", U = "#365C9A",
  S = "#E68613", R = "#C44E52",
  A.E = "#702D53", A.M = "#E58C72", A.U = "#D77DB0",
  A.S = "#A75D19", A.R = "#EFBB90", E.M = "#244E40",
  E.U = "#6F85CA", E.S = "#9A9C38", E.R = "#9DC8E5", S.R = "#76B99A"
)
region.names <- c("Africa", "Eurasia-Oriental", "Madagascar", "Australia",
                  "Sunken Islands", "Mauritius")

# Display settings ----
# 271 labels need a large vector figure; keep specimen identifiers for traceability.
tip.cex <- 0.48
node.cex <- 0.12
height <- max(node.depth.edgelength(phy))
phy$tip.label <- gsub("_", " ", phy$tip.label)
plots.dir <- "asr/plots"
dir.create(plots.dir, recursive = TRUE, showWarnings = FALSE)

# Save PDF and PNG ----
for (format in c("pdf", "png")) {
  filename <- file.path(plots.dir, paste0(output.name, ".", format))
  if (format == "pdf") pdf(filename, width = 22, height = 14, useDingbats = FALSE)
  if (format == "png") png(filename, width = 22, height = 14, units = "in", res = 220)
  layout(matrix(c(1, 2), ncol = 1), heights = c(11.8, 2.2))
  par(mar = c(0.3, 0.3, 0.3, 0.3), family = "sans")
  # ape leaves one tip-spacing gap in fan plots; compensate for an exact semicircle.
  plot(phy, type = "fan", open.angle = 180 - 360 / Ntip(phy),
       rotate.tree = 0, use.edge.length = TRUE, edge.color = "#606060",
       edge.width = 0.65, font = 3, cex = tip.cex, label.offset = height * .015,
       tip.color = "#252525", underscore = TRUE,
       x.lim = c(-1.36, 1.36) * height, y.lim = c(-.10, 1.35) * height)
  # Radial time guides: root-to-tip distance converted to age before present.
  # This age scale requires an ultrametric tree.
  stopifnot(is.ultrametric(phy))
  ages <- seq(0, height, by = 20)
  angle <- seq(0, pi, length.out = 500)
  for (age in ages) {
    radius <- height - age
    lines(radius * cos(angle), radius * sin(angle),
          col = "#BCC6CF", lty = 3, lwd = 1.2)
  }
  nodelabels(node = as.integer(colnames(node.probs)), pie = t(node.probs),
             piecol = colors[state.labels], cex = node.cex)
  tiplabels(pch = 21, bg = colors[state.labels[tip.states]], col = "white",
            cex = 0.75, lwd = 0.35)
  # Timeline along the right-hand base; 0 Ma is the present-day tip radius.
  timeline.y <- -.04 * height
  segments(0, timeline.y, height, timeline.y, col = "#6B7480", lwd = .7)
  axis.ages <- c(ages, height)
  for (age in axis.ages) {
    radius <- height - age
    segments(radius, timeline.y, radius, timeline.y - .008 * height,
             col = "#6B7480", lwd = .7)
    text(radius, timeline.y - .012 * height, format(round(age, 1), trim = TRUE),
         adj = c(.5, 1), cex = .8, col = "#46515C")
  }
  text(height / 2, -.095 * height, "Age (Ma)", cex = .9, col = "#46515C")

  # Complete geographic-state key ----
  par(mar = c(0, 1, 0, 1))
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  legend("top", legend = paste(state.labels[1:6], region.names, sep = "  "),
         fill = colors[1:6], border = NA, ncol = 3, bty = "n", cex = 1.0,
         x.intersp = .7, y.intersp = 1.1)
  legend(x = .5, y = .46, xjust = .5, legend = state.labels[7:16],
         fill = colors[7:16], border = NA, ncol = 5, bty = "n", cex = 1.0,
         x.intersp = .7, y.intersp = 1.1)
  text(.5, .08, "Node pies: marginal ancestral-state probabilities; tip circles: observed ranges",
       cex = .9, col = "#555555")
  dev.off()
}
# The PDF is vector-based: labels and node pies remain sharp when enlarged.
