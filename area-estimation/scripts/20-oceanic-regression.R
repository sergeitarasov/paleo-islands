# Run from paleo-islands. Retain all records classified as Oceanic.
library(readxl)
library(ggplot2)

# Data ----
d <- as.data.frame(read_excel('data/richness/island-richness-area.xlsx', sheet = 'Island isolation'))
d$log_richness <- log(d$n_sp)
d$log_area <- log(d$area_km2)
d$log_mainland <- log1p(d$distance_mainland_km)
d$group <- factor(ifelse(d$region == 'Caribbean', 'Caribbean', 'Other islands'))
d$type <- factor(d$type)
ocean <- d[d$type == 'Oceanic', ]
stopifnot(all(complete.cases(ocean[, c('log_richness','log_area','log_mainland')])))
out <- 'area-estimation/results/oceanic'
plots.dir <- 'area-estimation/plots/oceanic'
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(plots.dir, recursive = TRUE, showWarnings = FALSE)
write.csv(ocean, file.path(out, 'regression-input.csv'), row.names = FALSE)

# Candidates ----
forms <- list('Area' = log_richness ~ log_area,
  'Area + isolation' = log_richness ~ log_area + log_mainland,
  'Area + group' = log_richness ~ log_area + group,
  'Area + group + isolation' = log_richness ~ log_area + group + log_mainland,
  'Area x group' = log_richness ~ log_area * group,
  'Area x group + isolation' = log_richness ~ log_area * group + log_mainland)
scores <- data.frame()
coefficients <- data.frame()
predictions <- data.frame()
models <- list()
for (model in names(forms)) {
  fit <- lm(forms[[model]], data = ocean)
  models[[model]] <- fit
  k <- attr(logLik(fit), 'df')
  scores <- rbind(scores, data.frame(model, n = nrow(ocean),
    AICc = AIC(fit) + 2*k*(k+1)/(nrow(ocean)-k-1),
    adjusted_R2 = summary(fit)$adj.r.squared))
  ct <- summary(fit)$coefficients
  ci <- confint(fit)
  coefficients <- rbind(coefficients, data.frame(model, term = rownames(ct),
    estimate = ct[,1], lower = ci[,1], upper = ci[,2], p = ct[,4]))
}

# Same held-out oceanic islands for both training datasets ----
# Each held-out island is removed from ALL training data.
cv.forms <- forms
cv.forms[['Area + island type']] <- log_richness ~ log_area + type
for (model in names(cv.forms)) {
  for (training_set in c('Oceanic only', 'All island types')) {
    if (model == 'Area + island type' && training_set == 'Oceanic only') next
    train_data <- ocean
    if (training_set == 'All island types') train_data <- d
    for (i in seq_len(nrow(ocean))) {
      train <- train_data[train_data$island != ocean$island[i], ]
      fit <- lm(cv.forms[[model]], data = train)
      mu <- as.numeric(predict(fit, newdata = ocean[i, ]))
      sigma <- sqrt(sum(residuals(fit)^2)/nobs(fit))
      predictions <- rbind(predictions, data.frame(model, training_set,
        island = ocean$island[i], observed_log = ocean$log_richness[i], predicted_log = mu,
        log_density = dnorm(ocean$log_richness[i], mu, sigma, log = TRUE)))
    }
  }
}
cv <- data.frame()
for (model in names(cv.forms)) {
  for (training_set in c('Oceanic only', 'All island types')) {
    z <- predictions[predictions$model == model & predictions$training_set == training_set, ]
    if (nrow(z) == 0) next
    cv <- rbind(cv, data.frame(model, training_set,
      log_RMSE = sqrt(mean((z$observed_log-z$predicted_log)^2)),
      negative_log_density = -mean(z$log_density)))
  }
}
for (i in seq_len(nrow(scores))) {
  z <- cv[cv$model == scores$model[i] & cv$training_set == 'Oceanic only', ]
  scores$log_RMSE[i] <- z$log_RMSE
}
scores$delta_AICc <- scores$AICc - min(scores$AICc)
scores <- scores[order(scores$AICc), ]
write.csv(scores, file.path(out, 'model-comparison.csv'), row.names = FALSE)
write.csv(cv, file.path(out, 'matched-oceanic-validation.csv'), row.names = FALSE)
write.csv(coefficients, file.path(out, 'coefficients.csv'), row.names = FALSE)
write.csv(predictions, file.path(out, 'held-out-predictions.csv'), row.names = FALSE)
saveRDS(models, file.path(out, 'models.rds'))

# Figures ----
theme_set(theme_classic(base_size = 12))
p <- ggplot(cv, aes(log_RMSE, model, colour = training_set)) +
  geom_point(size = 3, position = position_dodge(width = .5)) +
  scale_colour_manual(values = c('All island types' = '#D55E00', 'Oceanic only' = '#0072B2')) +
  labs(x = 'Held-out oceanic islands: log RMSE (lower is better)', y = NULL, colour = 'Training data') +
  theme(legend.position = 'bottom')
ggsave(file.path(plots.dir, 'predictive-comparison.pdf'), p, width = 10, height = 5)
ggsave(file.path(plots.dir, 'predictive-comparison.png'), p, width = 10, height = 5, dpi = 180)

# Best-AICc regression at a common isolation distance ----
best <- scores$model[1]
fit <- models[[best]]
grid <- expand.grid(log_area = seq(min(ocean$log_area), max(ocean$log_area), length.out = 150),
                    group = levels(ocean$group))
grid$log_mainland <- median(ocean$log_mainland)
if (best %in% c('Area', 'Area + isolation')) grid <- grid[grid$group == levels(ocean$group)[1], ]
pr <- predict(fit, newdata = grid, interval = 'confidence')
grid$median <- exp(pr[, 'fit'])
grid$lower <- exp(pr[, 'lwr'])
grid$upper <- exp(pr[, 'upr'])
q <- ggplot(ocean, aes(area_km2, n_sp, colour = group)) +
  geom_ribbon(data = grid, aes(x = exp(log_area), ymin = lower, ymax = upper),
              inherit.aes = FALSE, fill = 'grey60', alpha = .20) +
  geom_line(data = grid, aes(exp(log_area), median), inherit.aes = FALSE, colour = 'grey20', linewidth = .7) +
  geom_point(size = 2.5, alpha = .85) +
  scale_x_log10() + scale_y_log10() +
  scale_colour_manual(values = c('#D55E00', '#0072B2')) +
  scale_fill_manual(values = c('#D55E00', '#0072B2')) +
  labs(x = expression('Island area (km'^2*')'), y = 'Species richness', colour = NULL, fill = NULL,
       caption = paste('Model:', best, '| Lines: median richness; bands: 95% confidence intervals.',
                       'Isolation held at the oceanic-sample median where applicable.')) +
  theme(legend.position = 'top')
ggsave(file.path(plots.dir, 'regression.pdf'), q, width = 10, height = 6)
ggsave(file.path(plots.dir, 'regression.png'), q, width = 10, height = 6, dpi = 180)

# Report ----
report <- c('# Oceanic-only richness regression', '',
 paste('All', nrow(ocean), 'oceanic islands retained; no richness-based exclusions. Region counts:'),
 capture.output(table(ocean$region)), '',
 'Response: log richness; Gaussian errors. Area uses natural log; mainland distance uses log(1 + km). Group is Caribbean versus other islands. Age is omitted because of sparse, heterogeneous proxies. Geographic classifications were not changed.', '',
 '| Model | AICc | Adjusted R2 | LOO log RMSE |', '|---|---:|---:|---:|')
for (i in seq_len(nrow(scores))) report <- c(report,
 sprintf('| %s | %.2f | %.3f | %.3f |', scores$model[i], scores$AICc[i], scores$adjusted_R2[i], scores$log_RMSE[i]))
report <- c(report, '', '## Matched predictive comparison', '',
 'Both approaches predict the same 31 held-out oceanic islands. Oceanic-only training contains the other 30 oceanic islands; full-data training contains the other 60 islands. No held-out observation is used in fitting. AICc is compared only within the oceanic subset, never across datasets.', '',
 '| Model | Training data | Oceanic LOO log RMSE | Negative log predictive density |', '|---|---|---:|---:|')
for (i in seq_len(nrow(cv))) report <- c(report,
 sprintf('| %s | %s | %.3f | %.3f |', cv$model[i], cv$training_set[i], cv$log_RMSE[i], cv$negative_log_density[i]))
report <- c(report, '',
 'These are exploratory candidate comparisons, not nested validation of the model-selection procedure. Island-level validation may benefit from shared regional/archipelago history and does not establish transfer to an extinct archipelago. The prediction score uses plug-in Gaussian residual variance, matching preceding analyses. Curves show conditional median richness, not inverse area estimates; no Sunken Islands area estimates were updated. Age and the large-landmass alternative were not searched in this focused comparison.')
report <- c(report, '', '## Interpretation', '',
 'Area alone has the lowest AICc (54.76), adjusted R2 = 0.763 and LOO log RMSE = 0.551. Adding isolation is competitive (delta AICc = 1.53), but gives worse LOO RMSE (0.563). The most complex model reduces RMSE to 0.545, only about 1.1%, with delta AICc = 3.69. Prefer area alone as the parsimonious oceanic-only candidate.',
 'Restricting training to oceanic islands does not improve point prediction on the same held-out oceanic records: area-only RMSE is 0.551 versus 0.539 using all island types, and the previously selected full-data area x group + isolation model yields 0.518. However oceanic area-only predicts uncertainty better than full-data area-only by the plug-in log-density score (0.909 versus 1.027; lower is better). Thus different metrics need not agree; there is no clear general predictive benefit from discarding the other islands. Oceanic-only remains a useful biological sensitivity analysis.')
writeLines(report, 'area-estimation/reports/oceanic-regression.md')
print(scores)
print(cv)
