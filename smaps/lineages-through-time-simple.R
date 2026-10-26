# Setup ----
# Open the paleo-islands RStudio project and run this script section by section.
library(divextra)
library(ggplot2)

# Read maps ----
smaps.nesovin <- readRDS("smaps/maps/posterior-aic-saved-asr/smaps-nesovin-full.rds")
smaps.nesosisyphus <- readRDS("smaps/maps/posterior-aic-saved-asr/smaps-nesosisyphus-full.rds")

# Settings ----
# Count observed backbone AND extinct side lineages in any of these states.
states <- c("A.S", "S", "S.R")
time.step <- 0.1  # Ma; change to 0.5 or 1 for a coarser grid.
height.nesovin <- max(smaps.nesovin[[1]]$lineages$start.age)
height.nesosisyphus <- max(smaps.nesosisyphus[[1]]$lineages$start.age)
times <- seq(0, max(height.nesovin, height.nesosisyphus), by = time.step)
times <- sort(unique(c(times, height.nesovin, height.nesosisyphus, 9)))

# Nesovinsonia ----
ltt.nesovin <- lineages.through.time.simmap(
  smaps.nesovin,
  times = times[times <= height.nesovin],
  states = states, quantile.type = 7
)
ltt.nesovin$summary

# Nesosisyphus ----
ltt.nesosisyphus <- lineages.through.time.simmap(
  smaps.nesosisyphus,
  times = times[times <= height.nesosisyphus],
  states = states, quantile.type = 7
)
ltt.nesosisyphus$summary

# Combined richness ----
# Pair maps by their shared model and posterior parameter draw.
stopifnot(length(smaps.nesovin) == length(smaps.nesosisyphus))
for (i in seq_along(smaps.nesovin)) {
  stopifnot(identical(smaps.nesovin[[i]]$posterior.draw$model,
                      smaps.nesosisyphus[[i]]$posterior.draw$model),
            identical(smaps.nesovin[[i]]$posterior.draw$posterior.row,
                      smaps.nesosisyphus[[i]]$posterior.draw$posterior.row))
}
# Only combine ages covered by BOTH clades; missing older history is not zero.
common.times <- times[times <= min(height.nesovin, height.nesosisyphus)]
counts.combined <- ltt.nesovin$counts[match(common.times, ltt.nesovin$times), , drop = FALSE] +
  ltt.nesosisyphus$counts[match(common.times, ltt.nesosisyphus$times), , drop = FALSE]

richness.nesovin <- ltt.nesovin$summary
richness.nesosisyphus <- ltt.nesosisyphus$summary
richness.combined <- data.frame(
  time = common.times, mean = rowMeans(counts.combined),
  median = NA_real_, lower = NA_real_, upper = NA_real_,
  n.maps = ncol(counts.combined)
)
# Sum counts first; do not add separate medians or interval endpoints.
for (i in seq_along(common.times)) {
  richness.combined$median[i] <- median(counts.combined[i, ])
  richness.combined$lower[i] <- quantile(counts.combined[i, ], .025, type = 7)
  richness.combined$upper[i] <- quantile(counts.combined[i, ], .975, type = 7)
}
richness.combined

# Full time range for the combined display ----
# Before the Nesosisyphus root, only Nesovinsonia is reconstructed.
# This is not an assumption that the older Nesosisyphus richness was zero.
boundary <- max(common.times)
older <- ltt.nesovin$times > boundary
counts.display <- rbind(counts.combined, ltt.nesovin$counts[older, , drop = FALSE])
times.display <- c(common.times, ltt.nesovin$times[older])

# Exact numerical summaries ----
counts <- list(Nesovinsonia = ltt.nesovin$counts,
               Nesosisyphus = ltt.nesosisyphus$counts, Combined = counts.display)
ages <- list(Nesovinsonia = ltt.nesovin$times,
             Nesosisyphus = ltt.nesosisyphus$times, Combined = times.display)
results <- list()
for (name in names(counts)) {
  n <- counts[[name]]
  x <- data.frame(time = ages[[name]], mean = rowMeans(n), median = 0,
                  lower = 0, upper = 0, minimum = 0, maximum = 0,
                  coverage = name)
  for (i in seq_len(nrow(n))) {
    x$median[i] <- median(n[i, ])
    x$lower[i] <- quantile(n[i, ], .025, type = 7)
    x$upper[i] <- quantile(n[i, ], .975, type = 7)
    x$minimum[i] <- min(n[i, ])
    x$maximum[i] <- max(n[i, ])
  }
  if (name == "Combined") {
    x$coverage <- ifelse(x$time > boundary, "Nesovinsonia only", "Both clades")
  }
  results[[name]] <- x
}

# Smooth the combined display only ----
# The exact summaries in results remain unchanged.
# Smoothed bounds are illustrative, not exact pointwise quantiles or extrema.
# Do not smooth across the 9 Ma epoch or the change in clade coverage.
combined.display <- results$Combined
segment <- cut(combined.display$time, breaks = c(-Inf, 9, boundary, Inf))
for (i in seq_len(nrow(combined.display))) {
  distance <- results$Combined$time - results$Combined$time[i]
  neighbours <- which(segment == segment[i] & abs(distance) <= 0.9)
  weights <- dnorm(distance[neighbours] / 0.3)
  for (column in c("mean", "lower", "upper", "minimum", "maximum")) {
    combined.display[i, column] <- weighted.mean(
      results$Combined[neighbours, column], weights
    )
  }
}
endpoints <- c(1, nrow(combined.display))
combined.display[endpoints, ] <- results$Combined[endpoints, ]

# Plot and save all three figures ----
# Grey: full map range. Blue: 95% map-distribution interval, not a CI of the mean.
# The y-axis uses a square-root scale; tick labels remain lineage counts.
plots.dir <- "smaps/plots/posterior-aic-saved-asr"
dir.create(plots.dir, recursive = TRUE, showWarnings = FALSE)
plots <- list()
for (name in names(results)) {
  x <- results[[name]]
  subtitle <- name
  caption <- "Exact pointwise 95% interval and full range; no smoothing."
  ticks <- seq(0, max(x$time), by = 10)
  if (name == "Nesosisyphus") ticks <- seq(0, max(x$time), by = 5)
  if (name == "Combined") {
    x <- combined.display
    # Join adjacent grid slices visually; smoothing remains separate by coverage.
    x$coverage <- "Combined display"
    subtitle <- sprintf("Older than %.2f Ma: Nesovinsonia only; younger: both clades", boundary)
    caption <- "Display-smoothed bands and mean (bandwidth 0.3 Ma); numerical summaries remain exact."
  }
  ymax <- max(x$maximum) * 1.10
  figure <- ggplot(x, aes(x = time, group = coverage)) +
    geom_ribbon(aes(ymin = minimum, ymax = maximum, fill = "Full range")) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = "95% interval")) +
    geom_line(aes(y = mean, colour = "Mean"), linewidth = 0.8) +
    geom_vline(xintercept = 9, colour = "#B85C38", linetype = "dashed", linewidth = 0.7) +
    annotate("text", x = 9.8, y = ymax * 0.95,
             label = "Mauritius emergence: 9 Ma", hjust = 1, size = 3.4, colour = "#8B442B") +
    scale_x_reverse(limits = c(max(x$time), 0), breaks = ticks,
                    expand = expansion(mult = c(.01, .01))) +
    scale_y_sqrt(limits = c(0, ymax), breaks = c(0, 1, 2, 5, 10, 20, 50, 100, 135),
                 expand = c(0, 0)) +
    scale_fill_manual(values = c("Full range" = "#EEEEEE", "95% interval" = "#ACCDDF"),
                      breaks = c("Full range", "95% interval")) +
    scale_colour_manual(values = c(Mean = "#174A6E")) +
    labs(x = "Age (Ma)", y = "Sunken-island lineages (square-root scale)",
         subtitle = subtitle, caption = caption) +
    theme_classic(base_size = 12) +
    theme(legend.position = "bottom", legend.title = element_blank(),
          plot.subtitle = element_text(size = 10), plot.caption = element_text(size = 8, hjust = 0))
  plots[[name]] <- figure
  print(figure)
  ggsave(file.path(plots.dir, paste0("fan-", tolower(name), ".pdf")), figure, width = 9, height = 5.7)
  ggsave(file.path(plots.dir, paste0("fan-", tolower(name), ".png")), figure, width = 9, height = 5.7, dpi = 200)
}
# Redisplay a figure in RStudio, for example: plots$Combined
# Exact statistics, including medians: results$Combined

# Combined plot with geological epochs and periods ----
# Install once if needed: install.packages("deeptime")
# Uses the ICS interval boundaries bundled with deeptime.
plot.geological <- plots$Combined +
  deeptime::coord_geo(
    xlim = c(max(results$Combined$time), 0),
    dat = list("epochs", "periods"), pos = list("bottom", "bottom"),
    abbrv = TRUE, size = 3, height = grid::unit(1.3, "lines"),
    alpha = 0.65, skip = "Holocene", center_end_labels = TRUE
  )
plots$Combined.geological <- plot.geological
print(plot.geological)
ggsave(file.path(plots.dir, "fan-combined-geological.pdf"),
       plot.geological, width = 9, height = 6.3)
ggsave(file.path(plots.dir, "fan-combined-geological.png"),
       plot.geological, width = 9, height = 6.3, dpi = 200)
