# Bayesian inverse prediction, preserving the upstream richness distribution.
library(ggplot2)
models <- readRDS("area-estimation/results/bayesian-selected21-10ma/regression-posteriors.rds")
frequency <- read.csv("area-estimation/results/inputs/richness-10ma-frequency.csv")
d <- read.csv("area-estimation/results/bayesian-selected21-10ma/calibration-data.csv")
d$x <- log(d$area_km2/10000)
# z is log10(area in km2). Grid integration includes endpoint weights.
z <- seq(0, 7, by = .01)
area <- 10^z
area_priors <- data.frame(z = z,
  log_uniform_10_1million = as.numeric(z >= 1 & z <= 6),
  log_uniform_100_100k = as.numeric(z >= 2 & z <= 5),
  log_uniform_1_10million = rep(1, length(z)),
  lognormal_10k = dnorm(z, 4, 1) * (z >= 1 & z <= 6),
  uniform_area_10_1million = 10^(z-6) * (z >= 1 & z <= 6))
for (name in names(area_priors)[-1]) {
  positive <- which(area_priors[[name]] > 0)
  area_priors[[name]][range(positive)] <- area_priors[[name]][range(positive)]/2
  area_priors[[name]] <- area_priors[[name]]/sum(area_priors[[name]])
}
# Summaries are calculated by interpolation of the grid CDF.
area_summary <- function(mass) {
  c(mean = sum(area*mass), median = 10^approx(cumsum(mass),z,xout=.5,ties="ordered")$y,
    lower = 10^approx(cumsum(mass),z,xout=.025,ties="ordered")$y,
    upper = 10^approx(cumsum(mass),z,xout=.975,ties="ordered")$y,
    p_below_100 = sum(mass[area<100]), p_above_100k = sum(mass[area>1e5]))
}
curves <- data.frame(); summaries <- data.frame(); numerical <- data.frame()
predictive <- data.frame()
started <- proc.time()[3]
for (name in names(models)) {
  model <- models[[name]]
  # Stratified across all four retained chains; full draws remain saved.
  idx <- round(seq(1, nrow(model$draws), length.out = 16000))
  draws <- model$draws[idx, ]
  p <- ncol(draws)-1
  # Conditional posterior predictive checks on the modern islands.
  X <- model.matrix(model$formula, d)
  mu <- X %*% t(draws[,seq_len(p)])
  set.seed(261026)
  replicated <- matrix(NA_real_, nrow(d), nrow(draws))
  for (j in seq_len(nrow(draws))) {
    if (startsWith(name,"simple")) replicated[,j] <- exp(rnorm(nrow(d), mu[,j], exp(draws[j,p+1])))
    else replicated[,j] <- rnbinom(nrow(d), mu = exp(mu[,j]), size = exp(draws[j,p+1]))
  }
  for (i in seq_len(nrow(d))) {
    predictive <- rbind(predictive, data.frame(model=name,island=d$island[i],observed=d$n_sp[i],
      median=median(replicated[i,]),lower=quantile(replicated[i,],.025),upper=quantile(replicated[i,],.975)))
  }
  for (type in "Selected 21 islands") {
    new <- data.frame(x = log(area/10000))
    X <- model.matrix(model$formula, new)
    eta <- X %*% t(draws[,seq_len(p)])
    scale <- matrix(exp(draws[,p+1]), nrow = length(z), ncol = nrow(draws), byrow = TRUE)
    likelihood <- matrix(0, length(z), nrow(frequency))
    likelihood_half <- likelihood
    for (j in seq_len(nrow(frequency))) {
      n <- frequency$n_species[j]
      if (startsWith(name, "simple")) value <- dnorm(log(n), eta, scale)
      else value <- dnbinom(n, mu = exp(eta), size = scale)
      likelihood[,j] <- rowMeans(value)
      likelihood_half[,j] <- rowMeans(value[,seq(1,16000,by=2)])
    }
    # Marginalize regression uncertainty FIRST. Normalize each richness value
    # separately, THEN mix with the unchanged map-derived richness frequencies.
    for (prior_name in names(area_priors)[-1]) {
      posterior_n <- likelihood * area_priors[[prior_name]]
      posterior_n <- sweep(posterior_n, 2, colSums(posterior_n), "/")
      mass <- as.vector(posterior_n %*% frequency$weight)
      half <- likelihood_half * area_priors[[prior_name]]
      half <- sweep(half, 2, colSums(half), "/")
      mass_half <- as.vector(half %*% frequency$weight)
      stats <- area_summary(mass)
      stats_half <- area_summary(mass_half)
      summaries <- rbind(summaries, data.frame(model=name,type=type,prior=prior_name,t(stats)))
      numerical <- rbind(numerical, data.frame(model=name,type=type,prior=prior_name,
        median_change_percent=100*(stats_half["median"]/stats["median"]-1),
        upper_change_percent=100*(stats_half["upper"]/stats["upper"]-1)))
      curves <- rbind(curves, data.frame(model=name,type=type,prior=prior_name,area_km2=area,
        density_log10_area=mass/.01, prior_density_log10_area=area_priors[[prior_name]]/.01))
    }
    cat(name,type,"done\n")
  }
}
write.csv(summaries,"area-estimation/results/bayesian-selected21-10ma/area-summary.csv",row.names=FALSE)
write.csv(curves,"area-estimation/results/bayesian-selected21-10ma/area-distributions.csv",row.names=FALSE)
write.csv(numerical,"area-estimation/results/bayesian-selected21-10ma/integration-check.csv",row.names=FALSE)
write.csv(predictive,"area-estimation/results/bayesian-selected21-10ma/posterior-predictive.csv",row.names=FALSE)


# Plots and numerical report ----
theme_set(theme_classic(base_size = 12))
base <- summaries[summaries$prior == "log_uniform_10_1million", ]
show <- curves[curves$model == "simple_prior1" & curves$prior == "log_uniform_10_1million" &
                 curves$area_km2 >= 10 & curves$area_km2 <= 1e6, ]
show$density_log10_area[show$area_km2 %in% c(10,1e6)] <- 2*show$density_log10_area[show$area_km2 %in% c(10,1e6)]
p <- ggplot(show, aes(area_km2, density_log10_area)) +
  geom_area(fill = "#0072B2", alpha = .2) + geom_line(colour = "#0072B2", linewidth = .8) +
  geom_vline(xintercept = base$median[1], linetype = 2, colour = "#0072B2") +
  scale_x_log10() + labs(x = expression("Effective area (km"^2*")"), y = "Density per log10 area")
ggsave("area-estimation/plots/bayesian-selected21-10ma/area-distribution.pdf", p, width = 8, height = 5)
ggsave("area-estimation/plots/bayesian-selected21-10ma/area-distribution.png", p, width = 8, height = 5, dpi = 180)
p <- ggplot(summaries, aes(median, prior, colour = model)) +
  geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y", width = .2,
                position = position_dodge(width = .5)) +
  geom_point(position = position_dodge(width = .5), size = 2) + scale_x_log10() +
  labs(x = expression("Median and 95% interval (km"^2*")"), y = NULL, colour = "Regression prior") +
  theme(legend.position = "bottom")
ggsave("area-estimation/plots/bayesian-selected21-10ma/prior-sensitivity.pdf", p, width = 10, height = 5)
ggsave("area-estimation/plots/bayesian-selected21-10ma/prior-sensitivity.png", p, width = 10, height = 5, dpi = 180)
for (m in unique(curves$model)) {
  for (pr in unique(curves$prior)) {
    mass <- curves$density_log10_area[curves$model == m & curves$prior == pr] * .01
    stopifnot(abs(sum(mass)-1) < 1e-8, all(is.finite(mass)))
  }
}
diag <- read.csv("area-estimation/results/bayesian-selected21-10ma/convergence.csv")
report <- c("# Sunken Islands effective area: selected 21-island calibration", "",
"Calibrated on 16 Caribbean records plus Réunion, Mauritius, Madagascar, New Caledonia and New Zealand (18 species). NZ area is provisionally 263,310 km2; see selected-21.md for source. No mainland-distance predictor is included. Model: ln(N) = alpha + beta ln(A/10,000) + epsilon. This is the pooled area-only model preferred in the requested subset comparison. The sampling distribution is lognormal; inverse prediction includes residual variation among islands, not only uncertainty in the fitted mean.", "",
"Used the full saved combined 2,000-map richness distribution at 10 Ma, including the tails. For each richness value, integrated the likelihood over regression draws, multiplied by the area prior and normalized. Then mixed these conditional area distributions with the unchanged map frequencies. This modular calculation preserves upstream richness probabilities rather than updating them with the island calibration data.", "",
"Regression priors: alpha ~ Normal(ln10, SD 2), beta ~ Normal(0.4, SD 0.5), sigma ~ half-Normal(SD 1.5). The wider-prior analysis doubles all these SDs. Four chains per prior, 5,000 warmup plus 25,000 retained iterations each; fixed-covariance random-walk Metropolis with scale adaptation during warmup only. Log-sigma target includes its Jacobian. Baseline area prior: log-uniform from 10 to 1,000,000 km2. No geology was used as a calibration target.", "",
"| Regression prior | Area prior | Mean km2 | Median km2 | 95% interval km2 |", "|---|---|---:|---:|---:|")
for (i in seq_len(nrow(summaries))) report <- c(report, sprintf("| %s | %s | %.0f | %.0f | %.0f–%.0f |", summaries$model[i], summaries$prior[i], summaries$mean[i], summaries$median[i], summaries$lower[i], summaries$upper[i]))
report <- c(report, "", sprintf("Maximum Rhat %.5f; minimum bulk/tail ESS %.0f / %.0f. Half-draw integration changes the median/upper limit by at most %.2f%% / %.2f%%. All area distributions normalize to one.", max(diag$rhat), min(diag$ess_bulk), min(diag$ess_tail), max(abs(numerical$median_change_percent)), max(abs(numerical$upper_change_percent))), "",
"Interpretation: effective richness-equivalent area conditional on the selected calibration subset, regression, priors and upstream maps. This is not an independent measurement of emerged geological area. Model-selection uncertainty and systematic richness/inventory error are not propagated. Earlier analyses and inputs were left unchanged.")
writeLines(report, "area-estimation/reports/bayesian-selected21-10ma.md")
print(base)
print(summaries)
