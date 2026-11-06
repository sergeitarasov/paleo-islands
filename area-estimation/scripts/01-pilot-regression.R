# Run from the paleo-islands project root in RStudio or with Rscript.
library(readxl)
library(MASS)
library(ggplot2)
library(flexmix)

# Data ----
d <- as.data.frame(read_excel("data/richness/island-richness-area.xlsx"))
d$region <- factor(d$region)
d$type <- factor(d$type)
d$log_area <- log(d$area_km2)
d$log_richness <- log(d$n_sp)
stopifnot(all(d$n_sp > 0), all(d$area_km2 > 0), !anyNA(d))
out <- "area-estimation/results/pilot-regression"
plots.dir <- "area-estimation/plots/pilot-regression"
logs.dir <- "area-estimation/logs/pilot-regression"
dir.create(plots.dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs.dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
write.csv(d, file.path(out, "input-data.csv"), row.names = FALSE)
write.csv(as.data.frame(table(d$region, d$type)), file.path(out, "group-counts.csv"), row.names = FALSE)

# Candidate models: power relationship, with grouped intercepts or slopes ----
formulas <- list(
  area = n_sp ~ log_area,
  region = n_sp ~ log_area + region,
  type = n_sp ~ log_area + type,
  region_type = n_sp ~ log_area + region + type,
  region_slopes = n_sp ~ log_area * region + type,
  type_slopes = n_sp ~ log_area * type + region
)
fits <- list()
comparison <- data.frame()
predictions <- data.frame()
set.seed(261026)

# Leave one island out: the same holdouts for every model.
# AICc is compared ONLY within response families (count vs log richness).
for (family in c("lognormal", "negative_binomial", "poisson")) {
  for (model in names(formulas)) {
    f <- formulas[[model]]
    if (family == "lognormal") f <- update(f, log_richness ~ .)
    if (family == "lognormal") fit <- lm(f, data = d)
    if (family == "negative_binomial") fit <- glm.nb(f, data = d, control = glm.control(maxit = 100))
    if (family == "poisson") fit <- glm(f, data = d, family = poisson())
    id <- paste(family, model, sep = ": ")
    fits[[id]] <- fit
    k <- attr(logLik(fit), "df")
    row <- data.frame(model = id, family = family, k = k,
                      AICc = AIC(fit) + 2*k*(k+1)/(nrow(d)-k-1), BIC = BIC(fit))
    for (i in seq_len(nrow(d))) {
      train <- d[-i, ]
      test <- d[i, , drop = FALSE]
      if (family == "lognormal") {
        cv <- lm(f, data = train)
        log_prediction <- as.numeric(predict(cv, test))
        prediction <- exp(log_prediction + sigma(cv)^2/2)
      }
      if (family == "negative_binomial") {
        cv <- glm.nb(f, data = train, control = glm.control(maxit = 100))
        prediction <- as.numeric(predict(cv, test, type = "response"))
        log_prediction <- log(prediction)
      }
      if (family == "poisson") {
        cv <- glm(f, data = train, family = poisson())
        prediction <- as.numeric(predict(cv, test, type = "response"))
        log_prediction <- log(prediction)
      }
      predictions <- rbind(predictions, data.frame(model = id, island = test$island,
                              observed = test$n_sp, predicted = prediction,
                              predicted_log = log_prediction))
    }
    comparison <- rbind(comparison, row)
  }
}

# Two latent groups: common residual variance avoids tiny-variance components.
# Labels are inferred from the training data, never from the held-out richness.
mix_driver <- FLXMRglmfix(family = "gaussian", varFix = TRUE)
mix <- stepFlexmix(log_richness ~ log_area, data = d, k = 2, nrep = 50,
                   model = mix_driver, control = list(iter.max = 1000, minprior = 0.1), verbose = FALSE)
if (inherits(mix, "stepFlexmix")) mix <- getModel(mix, "BIC")
fits[["mixture: two powers"]] <- mix
k <- attr(logLik(mix), "df")
comparison <- rbind(comparison, data.frame(model = "mixture: two powers", family = "lognormal",
  k = k, AICc = AIC(mix) + 2*k*(k+1)/(nrow(d)-k-1), BIC = BIC(mix)))
for (i in seq_len(nrow(d))) {
  cv <- stepFlexmix(log_richness ~ log_area, data = d[-i, ], k = 2, nrep = 15,
    model = mix_driver, control = list(iter.max = 1000, minprior = 0.1), verbose = FALSE)
  if (inherits(cv, "stepFlexmix")) cv <- getModel(cv, "BIC")
  pars <- parameters(cv)
  means <- pars[1, ] + pars[2, ] * d$log_area[i]
  weights <- prior(cv)
  prediction <- sum(weights * exp(means + pars[3, ]^2/2))
  predictions <- rbind(predictions, data.frame(model = "mixture: two powers", island = d$island[i],
    observed = d$n_sp[i], predicted = prediction, predicted_log = sum(weights * means)))
}
write.csv(data.frame(island = d$island, region = d$region, type = d$type,
                    posterior(mix)), file.path(out, "mixture-membership.csv"), row.names = FALSE)

# Scores: log-scale prediction of a typical island, and raw mean-count error ----
for (i in seq_len(nrow(comparison))) {
  p <- predictions[predictions$model == comparison$model[i], ]
  comparison$log_RMSE[i] <- sqrt(mean((log(p$observed)-p$predicted_log)^2))
  comparison$log_mean_RMSE[i] <- sqrt(mean((log(p$observed)-log(p$predicted))^2))
  comparison$count_RMSE[i] <- sqrt(mean((p$observed-p$predicted)^2))
  comparison$count_MAE[i] <- mean(abs(p$observed-p$predicted))
}
comparison <- comparison[order(comparison$log_mean_RMSE), ]
write.csv(comparison, file.path(out, "model-comparison.csv"), row.names = FALSE)
write.csv(predictions, file.path(out, "leave-one-out-predictions.csv"), row.names = FALSE)
saveRDS(fits, file.path(out, "models.rds"))

# Inspect the grouped power model, including influence and sensitivity ----
fit <- fits[["lognormal: region_type"]]
coef_table <- cbind(coef(summary(fit)), confint(fit))
write.csv(coef_table, file.path(out, "grouped-power-coefficients.csv"))
write.csv(data.frame(island = d$island, cooks_distance = cooks.distance(fit),
  leverage = hatvalues(fit)), file.path(out, "influence.csv"), row.names = FALSE)
sets <- list(all = d,
  individual_islands = d[!d$island %in% c("Aru Islands", "Kei Islands", "Guadeloupe", "New Caledonia"), ],
  without_mauritius_reunion = d[!d$island %in% c("Mauritius", "Réunion"), ])
sensitivity <- data.frame()
for (name in names(sets)) {
  fit_s <- lm(log_richness ~ log_area + region + type, data = sets[[name]])
  ci <- confint(fit_s, "log_area")
  sensitivity <- rbind(sensitivity, data.frame(subset = name, n = nrow(sets[[name]]),
    area_slope = coef(fit_s)["log_area"], lower = ci[1], upper = ci[2]))
}
write.csv(sensitivity, file.path(out, "sensitivity.csv"), row.names = FALSE)

# Publication-style pilot plots ----
theme_set(theme_classic(base_size = 11))
palette <- c(Australasian = "#0072B2", Caribbean = "#D55E00", "Indian Ocean" = "#009E73")
grid <- expand.grid(area_km2 = exp(seq(min(d$log_area), max(d$log_area), length.out = 180)),
                    region = levels(d$region), type = levels(d$type))
grid$log_area <- log(grid$area_km2)
ci <- predict(fit, grid, interval = "confidence")
grid$fit <- exp(ci[, 1]); grid$lower <- exp(ci[, 2]); grid$upper <- exp(ci[, 3])
# Restrict lines to each observed region/type combination's area range.
grid$keep <- FALSE
for (i in seq_len(nrow(grid))) {
  z <- d[d$region == grid$region[i] & d$type == grid$type[i], ]
  if (nrow(z) >= 2) grid$keep[i] <- grid$area_km2[i] >= min(z$area_km2) & grid$area_km2[i] <= max(z$area_km2)
}
p1 <- ggplot(d, aes(area_km2, n_sp, colour = region)) +
  geom_ribbon(data = grid[grid$keep, ], aes(y = fit, ymin = lower, ymax = upper, fill = region), alpha = .12, colour = NA) +
  geom_line(data = grid[grid$keep, ], aes(y = fit), linewidth = .65) + geom_point(size = 2.1) +
  facet_wrap(~type) + scale_x_log10(breaks = c(100, 10000, 1000000), labels = c("100", "10k", "1m")) + scale_y_log10() +
  scale_colour_manual(values = palette) + scale_fill_manual(values = palette) +
  labs(x = expression("Island area (km"^2*")"), y = "Species richness", colour = NULL, fill = NULL) +
  theme(legend.position = "bottom")
p2 <- ggplot(comparison, aes(log_mean_RMSE, reorder(model, -log_mean_RMSE), colour = family)) +
  geom_point(size = 2.5) + labs(x = "Prediction error: log-scale RMSE (lower is better)", y = NULL, colour = NULL) +
  theme(legend.position = "bottom")
p3 <- ggplot(predictions[predictions$model %in% c("lognormal: type", "negative_binomial: type_slopes", "mixture: two powers"), ],
             aes(observed, predicted)) + geom_abline(slope = 1, intercept = 0, colour = "grey60") +
  geom_point(size = 2, colour = "#0072B2") + scale_x_log10() + scale_y_log10() + facet_wrap(~model) +
  labs(x = "Observed richness", y = "Leave-one-out predicted mean richness")
plots <- list(species_area = p1, model_comparison = p2, predictive_check = p3)
for (name in names(plots)) {
  ggsave(file.path(plots.dir, paste0(name, ".pdf")), plots[[name]], width = 10, height = 5.5)
  ggsave(file.path(plots.dir, paste0(name, ".png")), plots[[name]], width = 10, height = 5.5, dpi = 180)
}
capture.output(sessionInfo(), file = file.path(logs.dir, "session-info.txt"))
print(comparison)
print(coef_table)
print(sensitivity)
print(summary(mix))
# Additional diagnostics and the count-model figure ----
nb <- fits[["negative_binomial: type_slopes"]]
write.csv(cbind(coef(summary(nb)), confint.default(nb)), file.path(out, "negative-binomial-coefficients.csv"))
ci <- predict(nb, grid, type = "link", se.fit = TRUE)
grid$fit <- exp(ci$fit)
grid$lower <- exp(ci$fit - 1.96*ci$se.fit)
grid$upper <- exp(ci$fit + 1.96*ci$se.fit)
p4 <- ggplot(d, aes(area_km2, n_sp, colour = region)) +
  geom_ribbon(data = grid[grid$keep, ], aes(y = fit, ymin = lower, ymax = upper, fill = region), alpha = .12, colour = NA) +
  geom_line(data = grid[grid$keep, ], aes(y = fit), linewidth = .65) + geom_point(size = 2.1) +
  facet_wrap(~type) + scale_x_log10(breaks = c(100, 10000, 1000000), labels = c("100", "10k", "1m")) +
  scale_y_log10() + scale_colour_manual(values = palette) + scale_fill_manual(values = palette) +
  labs(x = expression("Island area (km"^2*")"), y = "Species richness", colour = NULL, fill = NULL) +
  theme(legend.position = "bottom")
ggsave(file.path(plots.dir, "negative_binomial.pdf"), p4, width = 10, height = 5.5)
ggsave(file.path(plots.dir, "negative_binomial.png"), p4, width = 10, height = 5.5, dpi = 180)

# A sensitivity check, not a reason to discard unusual observations.
sets$without_aru <- d[d$island != "Aru Islands", ]
sets$without_three_influential <- d[!d$island %in% c("Aru Islands", "Madagascar", "Réunion"), ]
checks <- data.frame()
for (name in names(sets)) {
  for (model in names(formulas)) {
    f <- update(formulas[[model]], log_richness ~ .)
    z <- lm(f, data = sets[[name]])
    k <- attr(logLik(z), "df")
    n <- nrow(sets[[name]])
    checks <- rbind(checks, data.frame(subset = name, model = model, n = n,
      AICc = AIC(z) + 2*k*(k+1)/(n-k-1)))
  }
}
write.csv(checks, file.path(out, "sensitivity-model-selection.csv"), row.names = FALSE)
poisson_fit <- fits[["poisson: type_slopes"]]
cat("Poisson Pearson dispersion:", sum(residuals(poisson_fit, type = "pearson")^2)/df.residual(poisson_fit), "\n")
cat("Type-only power slope and confidence interval:\n")
print(cbind(coef(fits[["lognormal: type"]]), confint(fits[["lognormal: type"]])))
a <- predictions[predictions$model == "negative_binomial: type", ]
b <- predictions[predictions$model == "negative_binomial: type_slopes", ]
difference <- (log(a$observed)-log(a$predicted))^2 - (log(b$observed)-log(b$predicted))^2
cat("Mean paired log squared-error improvement and descriptive SE:", mean(difference), sd(difference)/sqrt(length(difference)), "\n")
