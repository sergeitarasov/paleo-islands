# Open the paleo-islands project and run from the project root.
library(readxl)
library(ggplot2)
library(ggrepel)

# Sensitivity: exclude Australasian richness < 6 and Reunion ----
islands <- as.data.frame(read_excel("data/richness/island-richness-area.xlsx"))
islands$log_area <- log(islands$area_km2)
islands$log_richness <- log(islands$n_sp)
# Compare the pooled power regressions ----
full <- islands
remove <- (full$region == "Australasian" & full$n_sp < 6) | full$island == "Réunion"
islands <- full[!remove, ]
fit.full <- lm(log_richness ~ log_area, data = full)
fit <- lm(log_richness ~ log_area, data = islands)

results.dir <- "area-estimation/results/filtered-regression"
plots.dir <- "area-estimation/plots/filtered-regression"
dir.create(results.dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots.dir, recursive = TRUE, showWarnings = FALSE)
write.csv(full[remove, ], file.path(results.dir, "excluded-islands.csv"), row.names = FALSE)
write.csv(islands, file.path(results.dir, "retained-islands.csv"), row.names = FALSE)
metrics <- data.frame()
fits <- list(All = fit.full, Filtered = fit)
for (name in names(fits)) {
  f <- fits[[name]]
  ci <- confint(f, "log_area")
  metrics <- rbind(metrics, data.frame(dataset = name, n = nobs(f),
    R2 = summary(f)$r.squared, adjusted_R2 = summary(f)$adj.r.squared,
    residual_SD_log = sigma(f), area_exponent = coef(f)["log_area"],
    exponent_lower = ci[1], exponent_upper = ci[2]))
}

# Fair prediction comparison: exactly the same 44 held-out islands.
# Full training includes the excluded islands, but never the held-out island.
predictions <- data.frame()
for (i in seq_len(nrow(islands))) {
  test <- islands[i, , drop = FALSE]
  train.full <- full[full$island != test$island, ]
  train.filtered <- islands[-i, ]
  f.full <- lm(log_richness ~ log_area, data = train.full)
  f.filtered <- lm(log_richness ~ log_area, data = train.filtered)
  pr.full <- predict(f.full, test, interval = "prediction")
  pr.filtered <- predict(f.filtered, test, interval = "prediction")
  predictions <- rbind(predictions, data.frame(island = test$island, observed = test$n_sp,
    log_observed = test$log_richness,
    full_log = pr.full[1], filtered_log = pr.filtered[1],
    full_mean = exp(pr.full[1] + sigma(f.full)^2/2),
    filtered_mean = exp(pr.filtered[1] + sigma(f.filtered)^2/2),
    full_covered = test$log_richness >= pr.full[2] & test$log_richness <= pr.full[3],
    filtered_covered = test$log_richness >= pr.filtered[2] & test$log_richness <= pr.filtered[3]))
}
cv <- data.frame(training = c("All", "Filtered"), n_test = nrow(predictions),
  log_RMSE = c(sqrt(mean((predictions$log_observed-predictions$full_log)^2)),
               sqrt(mean((predictions$log_observed-predictions$filtered_log)^2))),
  count_RMSE = c(sqrt(mean((predictions$observed-predictions$full_mean)^2)),
                 sqrt(mean((predictions$observed-predictions$filtered_mean)^2))),
  count_MAE = c(mean(abs(predictions$observed-predictions$full_mean)),
                mean(abs(predictions$observed-predictions$filtered_mean))),
  predictive_coverage = c(mean(predictions$full_covered), mean(predictions$filtered_covered)))
write.csv(metrics, file.path(results.dir, "fit-comparison.csv"), row.names = FALSE)
write.csv(cv, file.path(results.dir, "common-test-prediction.csv"), row.names = FALSE)
write.csv(predictions, file.path(results.dir, "held-out-predictions.csv"), row.names = FALSE)
saveRDS(fits, file.path(results.dir, "models.rds"))
print(metrics)
print(cv)


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
    caption = "Filtered pooled regression (44 islands); grey band: 95% confidence interval for fitted median richness.") +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom", plot.caption = element_text(size = 9, hjust = 0),
        plot.margin = margin(12, 20, 12, 12))
if (interactive()) print(figure)

dir.create("area-estimation/plots/filtered-regression", recursive = TRUE, showWarnings = FALSE)
ggsave("area-estimation/plots/filtered-regression/filtered-islands-labelled.pdf", figure, width = 13, height = 10)
ggsave("area-estimation/plots/filtered-regression/filtered-islands-labelled.png", figure, width = 13, height = 10, dpi = 180)
write.csv(cbind(coef(summary(fit)), confint(fit)),
  "area-estimation/results/filtered-regression/pooled-regression-coefficients.csv")

# Overlay full and filtered fits for direct comparison ----
comparison.curves <- data.frame()
for (name in names(fits)) {
  grid <- data.frame(log_area = seq(min(full$log_area), max(full$log_area), length.out = 300))
  grid$area_km2 <- exp(grid$log_area)
  pr <- predict(fits[[name]], grid, interval = "confidence")
  grid$richness <- exp(pr[,1]); grid$lower <- exp(pr[,2]); grid$upper <- exp(pr[,3])
  grid$fit <- name
  comparison.curves <- rbind(comparison.curves, grid)
}
full$status <- ifelse(remove, "Excluded", "Retained")
p <- ggplot(full, aes(area_km2, n_sp)) +
  geom_point(aes(shape = status), colour = "grey45", size = 2) +
  geom_ribbon(data = comparison.curves, aes(ymin = lower, ymax = upper, fill = fit, x = area_km2),
    inherit.aes = FALSE, alpha = .12) +
  geom_line(data = comparison.curves, aes(area_km2, richness, colour = fit), inherit.aes = FALSE, linewidth = .8) +
  scale_x_log10() + scale_y_log10() +
  scale_colour_manual(values = c(All = "#0072B2", Filtered = "#D55E00")) +
  scale_fill_manual(values = c(All = "#0072B2", Filtered = "#D55E00")) +
  scale_shape_manual(values = c(Excluded = 4, Retained = 16)) +
  labs(x = expression("Island area (km"^2*")"), y = "Species richness", colour = NULL, fill = NULL, shape = NULL) +
  theme_classic(base_size = 12) + theme(legend.position = "bottom")
ggsave(file.path(plots.dir, "full-versus-filtered.pdf"), p, width = 9, height = 7)
ggsave(file.path(plots.dir, "full-versus-filtered.png"), p, width = 9, height = 7, dpi = 180)

# Short, reproducible report ----
report <- c("# Richness-based exclusion sensitivity", "",
  "Excluded Australasian records with fewer than six species and Reunion: 16 + 1 records; 44 of 61 remain. The source workbook and previous analyses are unchanged.", "",
  sprintf("Adjusted R²: %.3f to %.3f. Log residual SD: %.3f to %.3f.", metrics$adjusted_R2[1], metrics$adjusted_R2[2], metrics$residual_SD_log[1], metrics$residual_SD_log[2]),
  sprintf("Area exponent: %.3f (95%% CI %.3f–%.3f) before; %.3f (%.3f–%.3f) after filtering.",metrics$area_exponent[1],metrics$exponent_lower[1],metrics$exponent_upper[1],metrics$area_exponent[2],metrics$exponent_lower[2],metrics$exponent_upper[2]), "",
  "Predictive comparison uses the SAME 44 held-out retained islands, with each target removed from both training sets. This avoids comparing errors on different test islands. Selection is still based on observed richness, so the comparison is conditional on passing the filter.", "",
  sprintf("Log-scale leave-one-out RMSE: %.3f (full training) versus %.3f (filtered training), a %.1f%% reduction.",cv$log_RMSE[1],cv$log_RMSE[2],100*(1-cv$log_RMSE[2]/cv$log_RMSE[1])),
  sprintf("Count-scale RMSE: %.2f versus %.2f; MAE: %.2f versus %.2f. Held-out 95%% prediction-interval coverage: %.1f%% versus %.1f%%.",cv$count_RMSE[1],cv$count_RMSE[2],cv$count_MAE[1],cv$count_MAE[2],100*cv$predictive_coverage[1],100*cv$predictive_coverage[2]), "",
  "A better fit after this filter does not demonstrate that the excluded records are erroneous or understudied. Removing low response values selectively from one region changes the sample and can bias its richness–area relationship. Confirm inventory completeness independently before using the filtered fit as the primary calibration. AIC values are not compared across different datasets. These confidence intervals do not account for data-driven selection of the exclusion rule.", "",
  "This run refits the pooled power regression requested for the single-panel plot. It does not refit the grouped Bayesian area-estimation models; their previous results remain unchanged. The reconstructed Sunken Islands richness is concentrated near the exclusion threshold, making extrapolation of this filtered calibration particularly sensitive to the inventory assumption.", "",
  'Run: source("area-estimation/scripts/06-filtered-regression.R") from the repository root. Results and figures are in results/filtered-regression and plots/filtered-regression.')
writeLines(report, "area-estimation/reports/filtered-regression.md")
