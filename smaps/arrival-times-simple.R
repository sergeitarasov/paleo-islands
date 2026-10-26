# First recorded arrivals on extant-lineage histories ----
# Run from the paleo-islands project root. Side lineages were pruned previously.
library(divextra)
library(ggplot2)

# Read observed-only maps ----
nesovin <- readRDS("smaps/maps/posterior-aic-saved-asr/nesovin-observed.rds")
nesosisyphus <- readRDS("smaps/maps/posterior-aic-saved-asr/nesosisyphus-observed.rds")
# To recreate these objects from full maps, prune.unobserved.simmap removes
# generated side lineages while preserving transitions on the retained paths.

# Exact directed transitions; both event types count equally ----
# First arrival = oldest matching event per map, not all events pooled.
# A map with no matching event has NA age, not age zero.
arrival.nesovin.sunken <- transitions.through.time.simmap(
  nesovin, breaks = seq(0, 65, by = 0.5), from = "A", to = "A.S",
  event.types = c("anagenetic", "cladogenetic"), level = .95, quantile.type = 7
)
arrival.nesovin.mauritius <- transitions.through.time.simmap(
  nesovin, breaks = seq(0, 65, by = 0.5), from = "S", to = "S.R",
  event.types = c("anagenetic", "cladogenetic"), level = .95, quantile.type = 7
)
arrival.nesosisyphus.sunken <- transitions.through.time.simmap(
  nesosisyphus, breaks = seq(0, 65, by = 0.5), from = "A", to = "A.S",
  event.types = c("anagenetic", "cladogenetic"), level = .95, quantile.type = 7
)
arrival.nesosisyphus.mauritius <- transitions.through.time.simmap(
  nesosisyphus, breaks = seq(0, 65, by = 0.5), from = "S", to = "S.R",
  event.types = c("anagenetic", "cladogenetic"), level = .95, quantile.type = 7
)

# Summaries: conditional on observing the specified transition ----
arrivals <- list(arrival.nesovin.sunken, arrival.nesovin.mauritius,
                 arrival.nesosisyphus.sunken, arrival.nesosisyphus.mauritius)
clades <- c("Nesovinsonia", "Nesovinsonia", "Nesosisyphus", "Nesosisyphus")
routes <- rep(c("Sunken islands: A -> A.S", "Mauritius: S -> S.R"), 2)
age.tables <- list()
summary.tables <- list()
for (i in seq_along(arrivals)) {
  a <- arrivals[[i]]
  age <- a$per.map$first.age
  finite <- age[is.finite(age)]
  age.tables[[i]] <- data.frame(clade = clades[i], route = routes[i],
                                map = seq_along(age), age = age)
  summary.tables[[i]] <- data.frame(
    clade = clades[i], route = routes[i], n.maps = length(age),
    n.arrivals = length(finite), n.without = sum(!is.finite(age)),
    p.arrival = length(finite) / length(age),
    mean = if (length(finite)) mean(finite) else NA_real_,
    median = if (length(finite)) median(finite) else NA_real_,
    lower = if (length(finite)) quantile(finite, .025, type = 7) else NA_real_,
    upper = if (length(finite)) quantile(finite, .975, type = 7) else NA_real_,
    n.anagenetic = sum(a$per.map$n.anagenetic),
    n.cladogenetic = sum(a$per.map$n.cladogenetic),
    n.root.in.destination = a$overall$n.root.in.destination,
    n.root.unknown = sum(is.na(a$per.map$root.in.destination))
  )
  stopifnot(all(a$per.map$n.outside == 0))
}
arrival.ages <- do.call(rbind, age.tables)
arrival.summary <- do.call(rbind, summary.tables)
rownames(arrival.summary) <- NULL
print(arrival.summary)

# Histograms preserve the empirical distributions without KDE smoothing ----
# Age intervals are 2.5th-97.5th map quantiles, not CIs for the mean.
# Initial occupancy of a destination is not counted as an observed transition;
# colonization before the root is outside the reconstructed history.
# These exact routes do not include other ways of arriving in S or R.
arrival.ages$clade <- factor(arrival.ages$clade, levels = c("Nesovinsonia", "Nesosisyphus"))
arrival.ages$route <- factor(arrival.ages$route, levels = unique(routes))
arrival.summary$clade <- factor(arrival.summary$clade, levels = levels(arrival.ages$clade))
arrival.summary$route <- factor(arrival.summary$route, levels = levels(arrival.ages$route))
arrival.summary$label <- sprintf("Mean: %.2f Ma\n95%% CI: %.2f-%.2f Ma",
                                 arrival.summary$mean, arrival.summary$lower, arrival.summary$upper)
plot.arrivals <- ggplot(arrival.ages[is.finite(arrival.ages$age), ], aes(x = age)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.5, boundary = 0,
                 fill = "#ACCDDF", colour = "white", linewidth = .15) +
  geom_vline(data = arrival.summary, aes(xintercept = lower, colour = "95% interval"), linetype = "dotted", linewidth = .7) +
  geom_vline(data = arrival.summary, aes(xintercept = upper, colour = "95% interval"), linetype = "dotted", linewidth = .7) +
  geom_vline(data = arrival.summary, aes(xintercept = mean, colour = "Mean"), linetype = "dotted", linewidth = .7) +
  geom_label(data = arrival.summary, aes(x = -Inf, y = Inf, label = label),
            inherit.aes = FALSE, hjust = 1.05, vjust = 1.1, size = 3, linewidth = 0, fill = "white") +
  facet_wrap(vars(clade, route), ncol = 2, scales = "free") +
  scale_x_reverse() +
  scale_colour_manual(values = c(Mean = "#174A6E", "95% interval" = "#B85C38")) +
  labs(x = "First recorded arrival (Ma before present)", y = "Conditional density",
       caption = "CI denotes the central 95% distribution of arrival ages, not a CI of the mean. Extant-lineage paths only.") +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        strip.background = element_blank(), strip.text = element_text(size = 10),
        plot.caption = element_text(size = 8, hjust = 0))
print(plot.arrivals)

# Save plots, ages and exact summaries ----
plots.dir <- "smaps/plots/posterior-aic-saved-asr"
summaries.dir <- "smaps/summaries/posterior-aic-saved-asr"
dir.create(plots.dir, recursive = TRUE, showWarnings = FALSE)
dir.create(summaries.dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(arrivals, file.path(summaries.dir, "first-arrivals-extant.rds"))
write.csv(arrival.ages, file.path(summaries.dir, "first-arrival-ages.csv"), row.names = FALSE)
write.csv(arrival.summary, file.path(summaries.dir, "first-arrival-summary.csv"), row.names = FALSE)
ggsave(file.path(plots.dir, "first-arrivals-extant.pdf"), plot.arrivals, width = 10, height = 7)
ggsave(file.path(plots.dir, "first-arrivals-extant.png"), plot.arrivals, width = 10, height = 7, dpi = 250)

# KDE distributions with the empirical 95% age interval filled ----
# Reflect at the temporal boundaries to avoid density beyond the root or 9 Ma.
# S -> S.R is allowed only in the recent epoch (0-9 Ma) in these fitted models.
# Bandwidths are selected from original ages, before reflection.
upper.ages <- c(max(nesovin[[1]]$lineages$start.age), 9,
                max(nesosisyphus[[1]]$lineages$start.age), 9)
density.tables <- list()
boundary.tables <- list()
interval.tables <- list()
for (i in seq_along(arrivals)) {
  ages <- arrivals[[i]]$per.map$first.age
  ages <- ages[is.finite(ages)]
  stopifnot(length(ages) > 1, all(ages >= 0), all(ages <= upper.ages[i]))
  bandwidth <- bw.nrd0(ages)
  reflected <- c(ages, -ages, 2 * upper.ages[i] - ages)
  curve <- density(reflected, bw = bandwidth, from = 0, to = upper.ages[i], n = 2048)
  curve$y <- 3 * curve$y
  area <- sum(diff(curve$x) * (head(curve$y, -1) + tail(curve$y, -1)) / 2)
  curve$y <- curve$y / area
  density.tables[[i]] <- data.frame(age = curve$x, density = curve$y,
                                    clade = clades[i], route = routes[i])
  boundary.tables[[i]] <- data.frame(
    age = c(0, upper.ages[i]), density = curve$y[c(1, length(curve$y))],
    clade = clades[i], route = routes[i]
  )
  # Insert exact interval boundaries, so the fill starts/ends at the quantiles.
  lower <- arrival.summary$lower[i]
  upper <- arrival.summary$upper[i]
  inside <- curve$x > lower & curve$x < upper
  interval.ages <- c(lower, curve$x[inside], upper)
  interval.tables[[i]] <- data.frame(
    age = interval.ages, density = approx(curve$x, curve$y, xout = interval.ages)$y,
    clade = clades[i], route = routes[i]
  )
}
arrival.density <- do.call(rbind, density.tables)
arrival.interval <- do.call(rbind, interval.tables)
arrival.boundaries <- do.call(rbind, boundary.tables)
arrival.density$clade <- factor(arrival.density$clade, levels = levels(arrival.ages$clade))
arrival.density$route <- factor(arrival.density$route, levels = levels(arrival.ages$route))
arrival.interval$clade <- factor(arrival.interval$clade, levels = levels(arrival.ages$clade))
arrival.interval$route <- factor(arrival.interval$route, levels = levels(arrival.ages$route))
arrival.boundaries$clade <- factor(arrival.boundaries$clade, levels = levels(arrival.ages$clade))
arrival.boundaries$route <- factor(arrival.boundaries$route, levels = levels(arrival.ages$route))
plot.arrivals.kde <- ggplot(arrival.density, aes(x = age, y = density)) +
  geom_area(fill = "#EEEEEE") +
  geom_area(data = arrival.interval, aes(fill = "95% age interval")) +
  geom_line(colour = "#174A6E", linewidth = .65) +
  geom_vline(data = arrival.summary, aes(xintercept = mean, colour = "Mean"), linetype = "dotted", linewidth = .65) +
  geom_label(data = arrival.summary, aes(x = -Inf, y = Inf, label = label),
            inherit.aes = FALSE, hjust = 1.05, vjust = 1.1, size = 3, linewidth = 0, fill = "white") +
  # Keep the Sunken-island boundary outline; Mauritius ends at 9 Ma without an extra line.
  geom_segment(data = arrival.boundaries[arrival.boundaries$route == "Sunken islands: A -> A.S", ],
               aes(x = age, xend = age, y = 0, yend = density),
               inherit.aes = FALSE, colour = "#174A6E", linewidth = .65) +
  facet_wrap(vars(clade, route), ncol = 2, scales = "free") +
  scale_x_reverse(expand = expansion(mult = c(.02, .02))) +
  scale_y_continuous(expand = expansion(mult = c(0, .12))) +
  scale_fill_manual(values = c("95% age interval" = "#ACCDDF")) +
  scale_colour_manual(values = c(Mean = "#174A6E")) +
  labs(x = "First recorded arrival (Ma before present)", y = "Conditional density",
       caption = "KDE bounded by the tree and model temporal limits. Blue: empirical 2.5th-97.5th age percentiles; summaries use unsmoothed ages.\nCI describes arrival-age uncertainty, not uncertainty in the mean; conditional on a matching event.") +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        strip.background = element_blank(), strip.text = element_text(size = 10),
        plot.caption = element_text(size = 8, hjust = 0))
print(plot.arrivals.kde)
ggsave(file.path(plots.dir, "first-arrivals-extant-kde.pdf"), plot.arrivals.kde, width = 10, height = 7)
ggsave(file.path(plots.dir, "first-arrivals-extant-kde.png"), plot.arrivals.kde, width = 10, height = 7, dpi = 250)

# Shared full timespan for direct comparison ----
# Empty space before a clade's root is outside its reconstructed history.
oldest.age <- max(nesovin[[1]]$lineages$start.age,
                  nesosisyphus[[1]]$lineages$start.age)
plot.arrivals.full <- plot.arrivals +
  facet_wrap(vars(clade, route), ncol = 2, scales = "free_y") +
  scale_x_reverse(limits = c(ceiling(oldest.age), 0), breaks = seq(0, 60, by = 10),
                  expand = expansion(mult = c(.01, .01)))
plot.arrivals.kde.full <- plot.arrivals.kde +
  facet_wrap(vars(clade, route), ncol = 2, scales = "free_y") +
  scale_x_reverse(limits = c(ceiling(oldest.age), 0), breaks = seq(0, 60, by = 10),
                  expand = expansion(mult = c(.01, .01)))
print(plot.arrivals.full)
print(plot.arrivals.kde.full)
ggsave(file.path(plots.dir, "first-arrivals-extant-full-span.pdf"), plot.arrivals.full, width = 11, height = 7)
ggsave(file.path(plots.dir, "first-arrivals-extant-full-span.png"), plot.arrivals.full, width = 11, height = 7, dpi = 250)
ggsave(file.path(plots.dir, "first-arrivals-extant-kde-full-span.pdf"), plot.arrivals.kde.full, width = 11, height = 7)
ggsave(file.path(plots.dir, "first-arrivals-extant-kde-full-span.png"), plot.arrivals.kde.full, width = 11, height = 7, dpi = 250)
