# Open the paleo-islands project and run from the project root.
library(readxl)
library(ggplot2)
library(ggrepel)

# All islands; one pooled power relationship ----
islands <- as.data.frame(read_excel("data/richness/island-richness-area.xlsx"))
islands$log_area <- log(islands$area_km2)
islands$log_richness <- log(islands$n_sp)
fit <- lm(log_richness ~ log_area, data = islands)

# Back-transformed fitted median and its 95% confidence band ----
curve <- data.frame(area_km2 = exp(seq(min(islands$log_area), max(islands$log_area), length.out = 300)))
curve$log_area <- log(curve$area_km2)
interval <- predict(fit, newdata = curve, interval = "confidence", level = .95)
curve$richness <- exp(interval[, "fit"])
curve$lower <- exp(interval[, "lwr"])
curve$upper <- exp(interval[, "upr"])

# Names are retained even in crowded regions of the plot ----
figure <- ggplot(islands, aes(area_km2, n_sp)) +
  geom_ribbon(data = curve, aes(x = area_km2, ymin = lower, ymax = upper),
              inherit.aes = FALSE, fill = "grey85", alpha = .65) +
  geom_line(data = curve, aes(area_km2, richness), inherit.aes = FALSE,
            colour = "grey25", linewidth = .85) +
  geom_point(aes(colour = region), size = 2.5) +
  geom_text_repel(aes(label = island, colour = region), seed = 26,
    size = 3.1, max.overlaps = Inf, max.iter = 100000, max.time = 10,
    box.padding = .45, point.padding = .2, min.segment.length = 0,
    segment.colour = "grey65", segment.size = .25, show.legend = FALSE) +
  scale_x_log10(breaks = c(100, 1000, 10000, 100000, 1000000),
    labels = c("100", "1,000", "10,000", "100,000", "1,000,000"),
    expand = expansion(mult = c(.12, .16))) +
  scale_y_log10(breaks = c(1, 2, 5, 10, 20, 50, 100, 200, 500),
    expand = expansion(mult = c(.22, .15))) +
  scale_colour_manual(values = c(Australasian = "#0072B2", Caribbean = "#D55E00", "Indian Ocean" = "#009E73")) +
  labs(x = expression("Island area (km"^2*")"), y = "Species richness", colour = NULL,
    caption = "Pooled log-log regression; grey band: 95% confidence interval for fitted median richness.") +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom", plot.caption = element_text(size = 9, hjust = 0),
        plot.margin = margin(12, 20, 12, 12))
if (interactive()) print(figure)

dir.create("area-estimation/plots/pilot-regression", recursive = TRUE, showWarnings = FALSE)
ggsave("area-estimation/plots/pilot-regression/all-islands-labelled.pdf", figure, width = 13, height = 10)
ggsave("area-estimation/plots/pilot-regression/all-islands-labelled.png", figure, width = 13, height = 10, dpi = 180)
write.csv(cbind(coef(summary(fit)), confint(fit)),
  "area-estimation/results/pilot-regression/pooled-regression-coefficients.csv")
