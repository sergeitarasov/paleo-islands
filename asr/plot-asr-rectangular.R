# Whole ingroup: rectangular biogeographic tree ----
# Run from the repository root: Rscript asr/plot-asr-rectangular.R
# Optional arguments: project directory, output directory (for a separate preview).
# Uses saved reconstruction probabilities; does not rerun inference.
library(ape)
library(readxl)

script_argument <- grep("^--file=", commandArgs(), value = TRUE)
arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) >= 1) {
  setwd(normalizePath(arguments[1]))
} else if (length(script_argument) == 1 && script_argument != "--file=-") {
  script_file <- sub("^--file=", "", script_argument)
  setwd(dirname(dirname(normalizePath(script_file))))
}
output_dir <- "asr/plots/full-asr-rectangular"
if (length(arguments) >= 2) output_dir <- arguments[2]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Read the existing 271-tip ASR before changing any display names.
x <- readRDS("asr/results/full.avg.rds")
phy <- ladderize(x$tree)
source_labels <- phy$tip.label
names_data <- as.data.frame(read_excel("data/biogeo.xlsx", sheet = "BioGeo"))
name_rows <- match(source_labels, names_data$to)
stopifnot(!anyNA(name_rows), !anyDuplicated(name_rows))
tip_data <- names_data[name_rows, c("from", "to", "tree_id", "geo")]
tip_data$source_label <- source_labels
tip_data$translated_label <- trimws(source_labels)
tip_data$translated_label <- gsub("[[:space:]]+", "_", tip_data$translated_label)
tip_data$translated_label <- gsub("_+", "_", tip_data$translated_label)

# Author-approved identifications and spelling changes (2026-11-03).
# Specimen IDs are retained; unnamed taxa in other genera are not relabelled.
approved_names <- data.frame(
  tree_id = c("STL18", "STL23", "STL40", "STL42", "STL9", "STL10223", "STL10224"),
  taxon = c("Sisyphus_sp1", "Sisyphus_sp6", "Sisyphus_sp8", "Sisyphus_sp10",
            "Sisyphus_gazanus", "Tanzanolus_leleupi", "Tanzanolus_leleupi")
)
approved_match <- match(tip_data$tree_id, approved_names$tree_id)
changed_tips <- which(!is.na(approved_match))
tip_data$translated_label[changed_tips] <- paste(
  approved_names$taxon[approved_match[changed_tips]],
  tip_data$tree_id[changed_tips], sep = "_")
tip_data$display_label <- gsub("_", " ", tip_data$translated_label, fixed = TRUE)
stopifnot(length(changed_tips) == 6L, !anyDuplicated(tip_data$translated_label))
for (tip_index in seq_len(nrow(tip_data))) {
  compact_id <- gsub("[^[:alnum:]]", "", tip_data$tree_id[tip_index])
  compact_label <- gsub("[^[:alnum:]]", "", tip_data$translated_label[tip_index])
  stopifnot(endsWith(compact_label, compact_id))
}

# Resolve state indices using ORIGINAL names. Node numbers survive ladderizing.
state_labels <- c("A", "E", "M", "U", "S", "R", "A.E", "A.M", "A.U",
                  "A.S", "A.R", "E.M", "E.U", "E.S", "E.R", "S.R")
stopifnot(identical(x$analysis$state.labels, state_labels))
node_probs <- x$node.probs[state_labels, , drop = FALSE]
node_ids <- as.integer(colnames(node_probs))
tip_states <- x$tip.states[source_labels]
n_tip <- Ntip(phy)
n_node <- Nnode(phy)
stopifnot(n_tip == 271L, n_node == 270L, is.ultrametric(phy))
stopifnot(setequal(node_ids, n_tip + seq_len(n_node)))
stopifnot(!anyNA(tip_states), all(tip_states %in% seq_along(state_labels)))
stopifnot(all(is.finite(node_probs)), all(node_probs >= 0), all(node_probs <= 1))
stopifnot(max(abs(colSums(node_probs) - 1)) < 1e-8)
stopifnot(identical(node_probs, x$node.probs))
tip_data$state <- state_labels[tip_states]
geo_lookup <- c(Afr = "A", OP = "E", Mada = "M", Aus = "U", Maur = "R")
stopifnot(identical(unname(geo_lookup[tip_data$geo]), unname(tip_data$state)))

# Use the current geographic palette in plot-asr-fan.R / plot-asr-clades-combined.R.
colours <- c(
  A = "#23A4B8", E = "#E3C64A", M = "#8B5AA5", U = "#365C9A",
  S = "#E68613", R = "#C44E52",
  A.E = "#702D53", A.M = "#E58C72", A.U = "#D77DB0",
  A.S = "#A75D19", A.R = "#EFBB90", E.M = "#244E40",
  E.U = "#6F85CA", E.S = "#9A9C38", E.R = "#9DC8E5", S.R = "#76B99A"
)
region_names <- c("Africa", "Eurasia-Oriental", "Madagascar", "Australia",
                  "Paleoislands", "Mauritius")

# Conventional rectangular layout: each internal node is midway between its
# children. Depths and ancestral-range probabilities are unchanged.
phy <- reorder(phy, "cladewise")
tip_order <- phy$edge[phy$edge[, 2] <= n_tip, 2]
y <- numeric(n_tip + n_node)
y[tip_order] <- rev(seq_len(n_tip))
postorder <- reorder(phy, "postorder")
for (parent in unique(postorder$edge[, 1])) {
  children <- phy$edge[phy$edge[, 1] == parent, 2]
  y[parent] <- mean(y[children])
}
depth <- node.depth.edgelength(phy)
tree_height <- max(depth[seq_len(n_tip)])
tip_data$plot_row <- y[seq_len(n_tip)]
phy$tip.label <- tip_data$translated_label
write.csv(tip_data, file.path(output_dir, "tip-labels.csv"), row.names = FALSE)
write.csv(data.frame(node = node_ids, node_tag = paste0("ASR_node_", node_ids),
                     node_age_Ma = tree_height - depth[node_ids],
                     t(node_probs), check.names = FALSE),
          file.path(output_dir, "node-probabilities.csv"), row.names = FALSE)
# Stable node tags keep probability-table joins valid after a tree reader
# renumbers nodes. These tags are identifiers, not support values.
export_tree <- phy
export_tree$node.label <- paste0("ASR_node_", n_tip + seq_len(n_node))
# Quote full labels to preserve parenthesized subgenera in plain Newick.
newick_tree <- export_tree
newick_tree$tip.label <- sprintf("tip_token_%03d", seq_len(n_tip))
newick_text <- write.tree(newick_tree, digits = 12)
for (tip_index in seq_len(n_tip)) {
  quoted_name <- gsub("'", "''", phy$tip.label[tip_index], fixed = TRUE)
  quoted_name <- paste0("'", quoted_name, "':")
  token <- paste0(newick_tree$tip.label[tip_index], ":")
  newick_text <- sub(token, quoted_name, newick_text, fixed = TRUE)
}
writeLines(newick_text, file.path(output_dir, "full-asr-rectangular-Ma.tre"))
write.nexus(export_tree, file = file.path(output_dir, "full-asr-rectangular-Ma.nex"),
            translate = TRUE, digits = 12)
newick_check <- read.tree(file.path(output_dir, "full-asr-rectangular-Ma.tre"))
newick_check$tip.label <- sub("^'", "", newick_check$tip.label)
newick_check$tip.label <- sub("'$", "", newick_check$tip.label)
newick_check$tip.label <- gsub("''", "'", newick_check$tip.label, fixed = TRUE)
stopifnot(isTRUE(all.equal(export_tree, newick_check, use.edge.length = TRUE, tolerance = 1e-9)))
nexus_check <- read.nexus(file.path(output_dir, "full-asr-rectangular-Ma.nex"))
stopifnot(identical(nexus_check$tip.label, phy$tip.label))
stopifnot(isTRUE(all.equal(export_tree, nexus_check, use.edge.length = TRUE, tolerance = 1e-9)))

# Draw pies with a physical radius, independent of the page aspect ratio.
draw_node_pie <- function(x, y, probabilities, radius_points, colours) {
  radius_x <- xinch(radius_points / 72)
  radius_y <- yinch(radius_points / 72)
  boundaries <- 2 * pi * c(0, cumsum(probabilities))
  for (state_index in seq_along(probabilities)) {
    if (probabilities[state_index] <= 0) next
    n_segments <- max(3, ceiling(120 * probabilities[state_index]))
    angle <- seq(boundaries[state_index], boundaries[state_index + 1], length.out = n_segments)
    polygon(c(x, x + radius_x * cos(angle)),
            c(y, y + radius_y * sin(angle)),
            col = colours[state_index], border = NA)
  }
  angle <- seq(0, 2 * pi, length.out = 100)
  lines(x + radius_x * cos(angle), y + radius_y * sin(angle),
        col = "#53606B", lwd = 0.35)
}

pdf_width <- 17
pdf_height <- 40
tip_points <- 9
node_radius_points <- 2.7
plot_file <- file.path(output_dir, "full-asr-rectangular.pdf")
cairo_pdf(plot_file, width = pdf_width, height = pdf_height,
          pointsize = 10, family = "Arial", onefile = TRUE)
layout(matrix(c(1, 2), ncol = 1), heights = c(pdf_height - 1.3, 1.3))
par(mar = c(3.5, 2, 3.7, 2), xaxs = "i", yaxs = "i", lend = "butt")
plot.new()
plot.window(xlim = c(-4, tree_height + 59), ylim = c(0.2, n_tip + 0.8))
ages <- seq(0, 100, 10)
axis_depths <- tree_height - ages
abline(v = axis_depths, col = "#E9EDF0", lwd = 0.5)
abline(v = tree_height, col = "#BCC7CF", lwd = 0.7)

for (parent in unique(phy$edge[, 1])) {
  children <- phy$edge[phy$edge[, 1] == parent, 2]
  segments(depth[parent], min(y[children]), depth[parent], max(y[children]),
           col = "#4C5663", lwd = 0.65)
}
segments(depth[phy$edge[, 1]], y[phy$edge[, 2]],
         depth[phy$edge[, 2]], y[phy$edge[, 2]], col = "#4C5663", lwd = 0.65)

# Bound pie radii by nearest symbol distance to avoid neighbouring pies
# obscuring one another. All state probabilities are retained, without cutoff.
symbol_x <- depth
symbol_y <- y
distance_inches <- as.matrix(dist(cbind(symbol_x / xinch(1), symbol_y / yinch(1))))
diag(distance_inches) <- Inf
node_radii <- rep(node_radius_points, length(node_ids))
for (node_index in seq_along(node_ids)) {
  closest_points <- min(distance_inches[node_ids[node_index], ]) * 72
  node_radii[node_index] <- min(node_radius_points, 0.44 * closest_points)
  draw_node_pie(depth[node_ids[node_index]], y[node_ids[node_index]],
                node_probs[, node_index], node_radii[node_index], colours[state_labels])
}
for (tip_index in seq_len(n_tip)) {
  closest_points <- min(distance_inches[tip_index, ]) * 72
  tip_radius <- min(2.4, 0.44 * closest_points)
  draw_node_pie(depth[tip_index], y[tip_index],
                as.numeric(seq_along(state_labels) == tip_states[tip_index]),
                tip_radius, colours[state_labels])
}

# Full names include every original specimen identifier. Keep voucher text
# together with the name so each displayed tip is unambiguous.
label_x <- tree_height + 1.6
text(label_x, y[seq_len(n_tip)], tip_data$display_label,
     adj = c(0, 0.5), cex = tip_points / 10, col = "#253744")
rightmost_label <- label_x + strwidth(tip_data$display_label, cex = tip_points / 10)
stopifnot(max(rightmost_label) < par("usr")[2])
axis(1, at = axis_depths, labels = ages, lwd = 0, lwd.ticks = 0.6,
     cex.axis = 0.95, col.axis = "#253744", col.ticks = "#78848F")
axis(3, at = axis_depths, labels = ages, lwd = 0, lwd.ticks = 0.6,
     cex.axis = 0.95, col.axis = "#253744", col.ticks = "#78848F")
mtext("Age (Ma)", side = 1, line = 2.2, at = tree_height / 2, col = "#253744")
mtext("Dung beetle biogeography", side = 3, line = 2.2, adj = 0,
      cex = 1.5, font = 2, col = "#253744")

# Compact key: single-region codes and all ten permitted two-region ranges.
par(mar = c(0.1, 2, 0.2, 2), xaxs = "i", yaxs = "i")
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1))
single_x <- rep(c(0.02, 0.35, 0.68), times = 2)
single_y <- rep(c(0.87, 0.63), each = 3)
for (state_index in seq_len(6)) {
  points(single_x[state_index], single_y[state_index], pch = 21,
         bg = colours[state_index], col = "#53606B", cex = 1.0, lwd = 0.4)
  text(single_x[state_index] + 0.015, single_y[state_index],
       paste(state_labels[state_index], region_names[state_index], sep = "  "),
       adj = c(0, 0.5), cex = 0.95, col = "#253744")
}
joint_x <- rep(seq(0.02, 0.82, length.out = 5), times = 2)
joint_y <- rep(c(0.35, 0.11), each = 5)
for (joint_index in seq_len(10)) {
  state_index <- joint_index + 6
  points(joint_x[joint_index], joint_y[joint_index], pch = 21,
         bg = colours[state_index], col = "#53606B", cex = 1.0, lwd = 0.4)
  text(joint_x[joint_index] + 0.015, joint_y[joint_index], state_labels[state_index],
       adj = c(0, 0.5), cex = 0.95, col = "#253744")
}
dev.off()

write.csv(data.frame(node = node_ids, radius_points = node_radii),
          file.path(output_dir, "node-pie-sizes.csv"), row.names = FALSE)
validation <- c(
  paste("Tips:", n_tip, "; node pies:", n_node),
  paste("Tip identifiers checked:", nrow(tip_data)),
  paste("Approved names changed:", length(changed_tips)),
  paste("Root age (Ma):", tree_height),
  paste("Maximum probability-sum error:", max(abs(colSums(node_probs) - 1))),
  paste("Pie radius range (points):", paste(round(range(node_radii), 3), collapse = " to ")),
  paste("PDF dimensions (inches):", pdf_width, "x", pdf_height),
  "All 16 states retained; no categories pooled or probabilities recalculated.",
  "Original saved RDS tree, probabilities, and spreadsheet unchanged."
)
writeLines(validation, file.path(output_dir, "validation.txt"))
write.csv(data.frame(file = c("asr/results/full.avg.rds", "data/biogeo.xlsx"),
                     md5 = unname(tools::md5sum(c("asr/results/full.avg.rds", "data/biogeo.xlsx")))),
          file.path(output_dir, "input-checksums.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(output_dir, "session-info.txt"))
cat(paste(validation, collapse = "\n"), "\n")
