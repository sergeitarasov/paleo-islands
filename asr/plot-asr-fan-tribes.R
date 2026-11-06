# Whole-tree ancestral ranges: semicircular fan ----
# Run from the paleo-islands project root. Uses saved ASR, not new inference.
library(ape)
library(sf)
# Install once if needed: install.packages("rnaturalearthdata")
library(rnaturalearthdata)

# Reconstruction ----
# To plot one model, change full.avg.rds to full.gr6.rds or full.gr7.rds.
x <- readRDS("asr/results/full.avg.rds")
output.name <- "full-asr-fan-tribes"
phy <- ladderize(x$tree)
# Revised tribal assignments; Nanosini and Epilissini follow the author's update.
# Source spellings are preserved to match the saved tree exactly.
classification <- read.csv("asr/tribe-assignments.csv", stringsAsFactors = FALSE)
genus <- sub("[_ ].*", "", phy$tip.label)
tribe <- classification$tribe[match(genus, classification$tree.genus)]
# Tip-specific assignment supplied by the author; preserve the source tip name.
# Other Arachnodes remain Epilissini.
tribe[phy$tip.label == "Arachnodes_aff metallicus _STL10196"] <- "Nanosini"
stopifnot(!anyNA(tribe))
state.labels <- c("A", "E", "M", "U", "S", "R", "A.E", "A.M", "A.U",
                  "A.S", "A.R", "E.M", "E.U", "E.S", "E.R", "S.R")
# Identify the focal genera using saved tip names (including the source typo).
focal.nesovin <- grep("^Nesovinsionia_|^Nesovinsonia_", phy$tip.label)
focal.nesosisyphus <- grep("^Nesosisyphus[_ ]", phy$tip.label)
focal.tips <- c(focal.nesovin, focal.nesosisyphus)
stopifnot(length(focal.nesovin) == 1, length(focal.nesosisyphus) == 4)
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

# Simplified biogeographic colour key ----
# Schematic realm boundaries, not country borders or a palaeogeographic map.
# E combines Palaearctic + Oriental; U includes Wallacea and Australasia.
# North Africa/Sahara belongs to E. The Sahara-Sahel transition is approximate.
# Realm context: https://doi.org/10.1071/SB22002
# Coastlines: Natural Earth 1:110m, distributed by rnaturalearthdata.
world <- rnaturalearthdata::countries110
world <- world[world$continent %in% c("Africa", "Europe", "Asia", "Oceania"), ]
# Dissolve countries BEFORE splitting islands, removing all political boundaries.
land <- st_make_valid(st_geometry(st_transform(world, 6933)))
land <- st_union(land)
land <- st_cast(land, "POLYGON")
land.area <- as.numeric(st_area(land)) / 1e6
land.centres <- st_coordinates(st_transform(st_centroid(land), 4326))
# Keep major landmasses, but remove small islands and northern archipelagos.
land <- land[land.area >= 15000 & land.centres[, 2] < 58 & land.centres[, 1] > -25]
land <- st_transform(land, 4326)
# Planar geometry in degrees is sufficient for this schematic inset.
land <- st_set_crs(land, NA)
land <- st_simplify(land, dTolerance = .15, preserveTopology = TRUE)

# Approximate southern Sahara boundary; do not use national boundaries.
africa.mask <- matrix(c(
  -20, -40, 52, -40, 52, 12, 44, 12, 39, 18, 33, 22,
  25, 21, 15, 18, 5, 18, -5, 18, -17, 20, -20, 20, -20, -40
), ncol = 2, byrow = TRUE)
africa.mask <- st_sfc(st_polygon(list(africa.mask)))
# Wallace's line runs east of Borneo and Bali, west of Sulawesi and Lombok.
# The Philippines remain Oriental. This explicitly avoids splitting Indonesia
# according to its national border or colouring all of Indonesia as Oriental.
australasia.mask <- matrix(c(
  100, -50, 180, -50, 180, 0, 130, 0, 127, 5,
  120, 5, 119, 2, 118, 0, 117, -3, 115.7, -8, 115.7, -11, 100, -12, 100, -50
), ncol = 2, byrow = TRUE)
australasia.mask <- st_sfc(st_polygon(list(australasia.mask)))
madagascar.mask <- st_as_sfc(st_bbox(c(xmin = 43, ymin = -26,
                                      xmax = 51, ymax = -11)))
map.africa <- st_intersection(land, africa.mask)
map.australasia <- st_intersection(land, australasia.mask)
map.madagascar <- st_intersection(land, madagascar.mask)

# Display settings ----
# 271 labels need a large vector figure; keep specimen identifiers for traceability.
tip.cex <- 0.48
node.cex <- 0.18
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
       edge.width = 0.65, show.tip.label = FALSE, font = 3, cex = tip.cex, label.offset = height * .015,
       tip.color = "#252525", underscore = TRUE,
       x.lim = c(-1.48, 1.48) * height, y.lim = c(-.10, 1.48) * height)
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
  # Solid red circles identify the five focal Mauritian tips.
  tiplabels(tip = focal.tips, pch = 19, col = colors["R"], cex = 1.2)
  # Brackets mark contiguous tribal membership in THIS tree.
  # Repeated labels are intentional when a tribe is split in this topology.
  coordinates <- get("last_plot.phylo", envir = ape::.PlotPhyloEnv)
  angles <- atan2(coordinates$yy[seq_len(Ntip(phy))],
                  coordinates$xx[seq_len(Ntip(phy))])
  order.tips <- order(angles)
  blocks <- rle(tribe[order.tips])
  ends <- cumsum(blocks$lengths)
  starts <- c(1, head(ends, -1) + 1)
  centres <- numeric(length(starts))
  for (j in seq_along(starts)) {
    centres[j] <- mean(range(angles[order.tips[starts[j]:ends[j]]]))
  }
  # Spread nearby names slightly while retaining their exact bracket endpoints.
  label.angles <- centres
  gap <- 2.5 * pi / 180
  for (j in 2:length(label.angles)) {
    label.angles[j] <- max(label.angles[j], label.angles[j - 1] + gap)
  }
  if (tail(label.angles, 1) > pi) {
    label.angles[length(label.angles)] <- pi
    for (j in (length(label.angles) - 1):1) {
      label.angles[j] <- min(label.angles[j], label.angles[j + 1] - gap)
    }
  }
  for (j in seq_along(starts)) {
    a <- range(angles[order.tips[starts[j]:ends[j]]])
    if (diff(a) == 0) a <- a + c(-.002, .002)
    arc <- seq(a[1], a[2], length.out = 80)
    lines(1.035 * height * cos(arc), 1.035 * height * sin(arc),
          col = "#4E5964", lwd = 1.2)
    for (end in a) {
      segments(1.025 * height * cos(end), 1.025 * height * sin(end),
               1.045 * height * cos(end), 1.045 * height * sin(end),
               col = "#4E5964", lwd = .7)
    }
    # Nesovinsonia has a separate horizontal callout below.
    if (blocks$values[j] == "Nesovinsoniini") next
    segments(1.045 * height * cos(centres[j]), 1.045 * height * sin(centres[j]),
             1.10 * height * cos(label.angles[j]), 1.10 * height * sin(label.angles[j]),
             col = "#A0A7AD", lwd = .5)
    degrees <- label.angles[j] * 180 / pi
    rotation <- if (degrees > 90) degrees - 180 else degrees
    alignment <- if (degrees > 90) 1 else 0
    label <- blocks$values[j]
    label.col <- "#26323D"
    label.font <- 1
    label.size <- 1.05
    text(1.11 * height * cos(label.angles[j]), 1.11 * height * sin(label.angles[j]),
         label, srt = rotation, adj = c(alignment, .5),
         cex = label.size, col = label.col, font = label.font)
  }

  # Horizontal genus label with a red leader to the Mauritian tip.
  segments(coordinates$xx[focal.nesovin], coordinates$yy[focal.nesovin],
           -1.18 * height, .74 * height, col = colors["R"], lwd = 1.0)
  text(-1.20 * height, .74 * height, "Nesovinsonia", adj = c(1, .5),
       cex = 1.15, font = 4, col = colors["R"])

  # Separate genus callout within Sisyphini, pointing to the four focal tips.
  focal.angle <- mean(angles[focal.nesosisyphus])
  segments(1.01 * height * cos(focal.angle), 1.01 * height * sin(focal.angle),
           -1.18 * height, .42 * height, col = colors["R"], lwd = 1.0)
  text(-1.20 * height, .42 * height, "Nesosisyphus", adj = c(1, .5),
       cex = 1.15, font = 4, col = colors["R"])

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
  legend("top", legend = paste0(region.names, " (", state.labels[1:6], ")"),
         fill = colors[1:6], border = NA, ncol = 3, bty = "n", cex = 1.0,
         x.intersp = .7, y.intersp = 1.1)
  legend(x = .5, y = .46, xjust = .5, legend = gsub(".", "-", state.labels[7:16], fixed = TRUE),
         fill = colors[7:16], border = NA, ncol = 5, bty = "n", cex = 1.0,
         x.intersp = .7, y.intersp = 1.1)
  text(.5, .08, "Node pies: marginal ancestral-state probabilities; tip circles: observed ranges",
       cex = .9, col = "#555555")
  # World inset occupies the unused upper-right corner of the fan panel.
  par(fig = c(.73, .98, .72, .95), mar = c(1.8, .2, 1.7, .2), new = TRUE)
  plot(land, col = colors["E"], border = NA,
       xlim = c(-22, 180), ylim = c(-50, 78), axes = FALSE, asp = 1,
       xaxs = "i", yaxs = "i")
  plot(map.africa, col = colors["A"], border = NA, add = TRUE)
  plot(map.australasia, col = colors["U"], border = NA, add = TRUE)
  plot(map.madagascar, col = colors["M"], border = NA, add = TRUE)
  # Mauritius is too small at this scale: enlarge its symbol and use a leader.
  points(57.55, -20.2, pch = 19, col = colors["R"], cex = 1.1)
  segments(60, -23, 85, -46, col = colors["R"], lwd = .8)
  text(88, -47, "Mauritius", adj = 0, cex = .75, col = colors["R"])
  # Schematic position on the Mascarene Plateau, not reconstructed coastlines.
  points(60, -11, pch = 19, col = colors["S"], cex = 1.3)
  segments(62, -11, 77, -29, col = colors["S"], lwd = .8)
  text(79, -29, "Sunken Islands", adj = 0, cex = .75, col = colors["S"])
  mtext("Island markers schematic; not to scale", side = 1, line = .3,
        cex = .65, col = "#666666")
  dev.off()
}
# The PDF is vector-based: labels and node pies remain sharp when enlarged.
