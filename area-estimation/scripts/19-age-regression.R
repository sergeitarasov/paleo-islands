# Run from paleo-islands. No richness-based or geographic exclusions.
library(readxl)
library(ggplot2)

# Data ----
d <- as.data.frame(read_excel('data/richness/island-richness-area.xlsx', sheet = 'Island isolation'))
d$log_richness <- log(d$n_sp)
d$log_area <- log(d$area_km2)
d$log_mainland <- log1p(d$distance_mainland_km)
d$log_age <- log(d$island_age_ma)
age <- d[!is.na(d$island_age_ma), ]
# Bounds are not exact ages. Non-oceanic ages are not zero.
samples <- list('All point-age proxies' = age,
                'Without conflicting Sumbawa' = age[age$island != 'Sumbawa', ])
forms <- list('Area' = log_richness ~ log_area,
              'Area + age' = log_richness ~ log_area + log_age,
              'Area + mainland' = log_richness ~ log_area + log_mainland,
              'Area + mainland + age' = log_richness ~ log_area + log_mainland + log_age)
out <- 'area-estimation/results/island-age'
plots.dir <- 'area-estimation/plots/island-age'
dir.create(plots.dir, recursive = TRUE, showWarnings = FALSE)
scores <- data.frame()
coefficients <- data.frame()
predictions <- data.frame()
for (sample in names(samples)) {
  a <- samples[[sample]]
  for (model in names(forms)) {
    fit <- lm(forms[[model]], data = a)
    predicted <- numeric(nrow(a))
    for (i in seq_len(nrow(a))) {
      training <- lm(forms[[model]], data = a[-i, ])
      predicted[i] <- predict(training, newdata = a[i, ])
    }
    k <- attr(logLik(fit), 'df')
    aicc <- NA_real_
    if (nrow(a) > k + 1) aicc <- AIC(fit) + 2*k*(k+1)/(nrow(a)-k-1)
    scores <- rbind(scores, data.frame(sample, model, n = nrow(a), AICc = aicc,
      adjusted_R2 = summary(fit)$adj.r.squared,
      LOO_log_RMSE = sqrt(mean((a$log_richness - predicted)^2))))
    ct <- summary(fit)$coefficients
    ci <- confint(fit)
    coefficients <- rbind(coefficients, data.frame(sample, model, term = rownames(ct),
      estimate = ct[,1], lower = ci[,1], upper = ci[,2], p = ct[,4]))
    predictions <- rbind(predictions, data.frame(sample, model, island = a$island,
      observed_log = a$log_richness, predicted_log = predicted))
  }
}
write.csv(age, file.path(out, 'regression-input.csv'), row.names = FALSE)
write.csv(scores, file.path(out, 'model-comparison.csv'), row.names = FALSE)
write.csv(coefficients, file.path(out, 'coefficients.csv'), row.names = FALSE)
write.csv(predictions, file.path(out, 'held-out-predictions.csv'), row.names = FALSE)

# Comparison plot ----
p <- ggplot(scores, aes(LOO_log_RMSE, model, colour = sample)) +
  geom_point(size = 3, position = position_dodge(width = .4)) +
  scale_colour_manual(values = c('#0072B2', '#D55E00')) +
  labs(x = 'Leave-one-island-out log RMSE (lower is better)', y = NULL, colour = NULL) +
  theme_classic(base_size = 12) + theme(legend.position = 'bottom')
ggsave(file.path(plots.dir, 'age-predictive-comparison.pdf'), p, width = 9, height = 4)
ggsave(file.path(plots.dir, 'age-predictive-comparison.png'), p, width = 9, height = 4, dpi = 180)
print(scores)
print(coefficients[coefficients$term == 'log_age', ])

# Short report ----
report <- c('# Age and isolation: unfiltered richness analysis', '',
'Isolation was refitted on all 61 islands, without the previous richness/geographic exclusions (script 18). Mainland distance improves the area-only model: AICc 171.81 to 164.35; leave-one-out log RMSE 0.969 to 0.915. For Caribbean/rest-specific slopes, mainland distance improves AICc 172.47 to 160.45 and log RMSE 0.955 to 0.865 (9.4% reduction). This is the best of the 12 isolation candidates by AICc and log RMSE. Distance to a landmass of at least 10,000 square km alone does not improve these baselines. Isolation does not improve held-out performance after accounting for island type. See isolation-regression.md for all candidates and uncertainty.', '',
'Age cannot be evaluated across all 61 islands: only seven have point-age proxies. The following comparisons use identical complete-case islands within each sample. The three one-sided bounds were not treated as exact ages, and missing ages were not imputed as zero. This is exploratory: volcanic construction, rock ages and reconstructed land appearance are heterogeneous covariates.', '',
'| Sample | Model | n | AICc | Adjusted R2 | LOO log RMSE |',
'|---|---|---:|---:|---:|---:|')
for (i in seq_len(nrow(scores))) {
 report <- c(report, sprintf('| %s | %s | %d | %.2f | %.3f | %.3f |',
   scores$sample[i], scores$model[i], scores$n[i], scores$AICc[i],
   scores$adjusted_R2[i], scores$LOO_log_RMSE[i]))
}
report <- c(report, '',
'Adding age worsens AICc and leave-one-out prediction on the seven-island subset. The age coefficient in area + age is 0.496 (ordinary model-based 95% CI -1.561 to 2.554; p = 0.540). Excluding conflicting Sumbawa gives the same qualitative result. AICc is undefined for the most complex six-island model because n-k-1 is zero. Small samples and high leverage make these predictive comparisons unstable. Do not interpret failure to improve fit as proof that geological age is biologically unimportant.', '',
'For now, mainland distance is the useful new candidate predictor; age is not supported by the available data. The nearest-landmass name and database ID identify the distance target and were not fitted as high-cardinality categorical predictors. The 10,000 km2 threshold is constant and cannot explain variation. Contemporary distances do not establish palaeodistances for Sunken Islands, and no Bayesian inverse area analysis was changed. No AIC comparison is made between the seven-island age models and the 61-island isolation models.')
writeLines(report, 'area-estimation/reports/age-isolation-comparison.md')
