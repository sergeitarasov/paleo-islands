# Reuse the fitted modern-island regression; change only the target age and richness.
dir.create("area-estimation/results/bayesian-isolation-scenarios-10ma",recursive=TRUE,showWarnings=FALSE)
dir.create("area-estimation/plots/bayesian-isolation-scenarios-10ma",recursive=TRUE,showWarnings=FALSE)
# Bayesian inverse prediction, preserving the upstream richness distribution.
library(ggplot2)
models <- readRDS("area-estimation/results/bayesian-isolation-scenarios/regression-posteriors.rds")
frequency <- read.csv("area-estimation/results/inputs/richness-10ma-frequency.csv")
d <- read.csv("area-estimation/results/bayesian-isolation-scenarios/calibration-data.csv")
d$x <- log(d$area_km2/10000)
d$group <- factor(d$group, levels = c("Caribbean", "Other islands"))
d$log_distance <- log1p(d$distance_mainland_km) - log1p(1000)
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
  idx <- round(seq(1, nrow(model$draws), length.out = 8000))
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
  for (distance in c(500, 1000, 1500, 2000, 2500)) {
    type <- paste(distance, "km from mainland")
    new <- data.frame(x = log(area/10000),
      group = factor("Other islands", levels = c("Caribbean", "Other islands")),
      log_distance = log1p(distance) - log1p(1000))
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
      likelihood_half[,j] <- rowMeans(value[,seq(1,8000,by=2)])
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
write.csv(summaries,"area-estimation/results/bayesian-isolation-scenarios-10ma/area-summary.csv",row.names=FALSE)
write.csv(curves,"area-estimation/results/bayesian-isolation-scenarios-10ma/area-distributions.csv",row.names=FALSE)
write.csv(numerical,"area-estimation/results/bayesian-isolation-scenarios-10ma/integration-check.csv",row.names=FALSE)
write.csv(predictive,"area-estimation/results/bayesian-isolation-scenarios-10ma/posterior-predictive.csv",row.names=FALSE)


# Scenario table and figure ----
base <- summaries[summaries$model == "simple_prior1" & summaries$prior == "log_uniform_10_1million", ]
base$distance_km <- c(500,1000,1500,2000,2500)
write.csv(base,"area-estimation/results/bayesian-isolation-scenarios-10ma/distance-summary.csv",row.names=FALSE)
show <- curves[curves$model=="simple_prior1" & curves$prior=="log_uniform_10_1million" & curves$area_km2>=10 & curves$area_km2<=1e6, ]
show$type <- factor(show$type, levels=base$type)
show$density_log10_area[show$area_km2 %in% c(10,1e6)] <- 2 * show$density_log10_area[show$area_km2 %in% c(10,1e6)]
p <- ggplot(show,aes(area_km2,density_log10_area,colour=type)) + geom_line(linewidth=.8) +
 scale_x_log10() + scale_colour_discrete(labels=c("500 km","1,000 km","1,500 km","2,000 km","2,500 km")) + theme_classic(base_size=12) +
 labs(x=expression("Effective area (km"^2*")"),y="Density per log10 area",colour="Distance to Africa") +
 theme(legend.position="bottom")
ggsave("area-estimation/plots/bayesian-isolation-scenarios-10ma/area-distributions.png",p,width=10,height=6,dpi=180)
ggsave("area-estimation/plots/bayesian-isolation-scenarios-10ma/area-distributions.pdf",p,width=10,height=6)
diag <- read.csv("area-estimation/results/bayesian-isolation-scenarios/convergence.csv")
report <- c("# Sunken Islands area under mainland-distance scenarios", "",
"All 61 unfiltered modern islands calibrate ln(N) = alpha + beta ln(A/10,000) + Caribbean/rest intercept and slope differences + gamma [ln(1+distance_km)-ln(1001)]. Sunken Islands use the non-Caribbean group. These are illustrative coast-to-Africa distance scenarios at 10 Ma, NOT a reconstructed palaeodistance distribution. No single scenario is selected or averaged.", "",
"The full saved 2,000-map combined richness distribution at 10 Ma is retained. Regression uncertainty is marginalized first within each richness value; each conditional area distribution is normalized before mixing with the unchanged map weights. Residual variation among islands is included. The baseline area prior is log-uniform on 10–1,000,000 km2. The lognormal observation model and continuous log-richness inverse likelihood follow previous Bayesian area analyses.", "",
"| Distance to Africa (km) | Median km2 | Mean km2 | 95% credible interval km2 |", "|---:|---:|---:|---:|")
for(i in seq_len(nrow(base))) report <- c(report,sprintf("| %d | %.0f | %.0f | %.0f–%.0f |",base$distance_km[i],base$median[i],base$mean[i],base$lower[i],base$upper[i]))
report <- c(report,"",sprintf("Maximum Rhat %.5f; minimum bulk/tail ESS %.0f / %.0f. Four chains per prior width; 5,000 warmup and 25,000 retained draws each. Integration uses 8,000 draws; the half-draw check changes median/upper limits by at most %.2f%% / %.2f%%.",max(diag$rhat),min(diag$ess_bulk),min(diag$ess_tail),max(abs(numerical$median_change_percent)),max(abs(numerical$upper_change_percent))),"",
"Priors: reference-group intercept Normal(ln10, SD2); reference area slope Normal(0.4, SD0.5); group intercept difference Normal(0, SD1.5); slope difference Normal(0, SD0.5); distance coefficient Normal(0, SD1.5); residual SD half-Normal(SD1.5). Wider-prior fit doubles all SDs. Rate-free, fixed-covariance random-walk Metropolis uses adaptation only in warmup; the log-sigma Jacobian is included. Area-prior and wider-regression-prior sensitivity are saved in area-summary.csv.","",
"These are effective richness-equivalent areas, conditional on geographic distance, calibration model, priors and upstream maps. They are not physical land-area measurements. Geographic distance could itself depend on island geometry/area; scenarios treat distance as externally fixed. Model-selection uncertainty and palaeodistance uncertainty are not propagated. The current model wins forward richness prediction, but the simpler area-only model had a smaller median inverse-recovery error in exploratory validation. Existing estimates were not overwritten.")
writeLines(report,"area-estimation/reports/bayesian-isolation-scenarios-10ma.md")
for(m in unique(curves$model)) for(t in unique(curves$type)) for(pr in unique(curves$prior)) {
 mass <- curves$density_log10_area[curves$model==m & curves$type==t & curves$prior==pr]*.01
 stopifnot(abs(sum(mass)-1)<1e-8,all(is.finite(mass)))
}
print(base)

old <- read.csv("area-estimation/results/bayesian-isolation-scenarios/distance-summary.csv")
comparison <- data.frame(distance_km=base$distance_km,median_9ma=old$median,median_10ma=base$median,
 mean_9ma=old$mean,mean_10ma=base$mean,lower_9ma=old$lower,lower_10ma=base$lower,upper_9ma=old$upper,upper_10ma=base$upper)
write.csv(comparison,"area-estimation/results/bayesian-isolation-scenarios-10ma/comparison-9-vs-10ma.csv",row.names=FALSE)
sources <- c("smaps/summaries/posterior-aic-saved-asr/richness-through-time.rds",
 "area-estimation/results/inputs/richness-10ma-frequency.csv",
 "area-estimation/results/bayesian-isolation-scenarios/regression-posteriors.rds")
write.csv(data.frame(source=sources,md5=unname(tools::md5sum(sources))),"area-estimation/results/bayesian-isolation-scenarios-10ma/source-files.csv",row.names=FALSE)
