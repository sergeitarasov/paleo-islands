# Bayesian inverse prediction, preserving the upstream richness distribution.
library(ggplot2)
models <- readRDS("area-estimation/results/bayesian-filtered/regression-posteriors.rds")
frequency <- read.csv("area-estimation/results/inputs/richness-9ma-frequency.csv")
d <- read.csv("area-estimation/results/bayesian-filtered/calibration-data.csv")
d$group <- factor(d$group, levels = c("Caribbean", "Other islands"))
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
  idx <- round(seq(1, nrow(model$draws), length.out = 4000))
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
  for (type in "Non-Caribbean") {
    new <- data.frame(x = log(area/10000),
      group = factor("Other islands", levels = levels(d$group)))
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
      likelihood_half[,j] <- rowMeans(value[,seq(1,4000,by=2)])
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
write.csv(summaries,"area-estimation/results/bayesian-filtered/area-summary.csv",row.names=FALSE)
write.csv(curves,"area-estimation/results/bayesian-filtered/area-distributions.csv",row.names=FALSE)
write.csv(numerical,"area-estimation/results/bayesian-filtered/integration-check.csv",row.names=FALSE)
write.csv(predictive,"area-estimation/results/bayesian-filtered/posterior-predictive.csv",row.names=FALSE)


# Compare with previous estimates under identical area priors ----
old_summary <- read.csv("area-estimation/results/bayesian/area-summary.csv")
old_curves <- read.csv("area-estimation/results/bayesian/area-distributions.csv")
old_summary <- old_summary[old_summary$type == "Oceanic" & grepl("prior1", old_summary$model), ]
old_curves <- old_curves[old_curves$type == "Oceanic" & grepl("prior1", old_curves$model), ]
old_summary$analysis <- ifelse(startsWith(old_summary$model,"simple"),"Previous simple: 61 islands","Previous flexible: 61 islands")
old_curves$analysis <- ifelse(startsWith(old_curves$model,"simple"),"Previous simple: 61 islands","Previous flexible: 61 islands")
new_summary <- summaries[summaries$model == "simple_prior1", ]
new_curves <- curves[curves$model == "simple_prior1", ]
new_summary$analysis <- "Filtered grouped: 39 islands"
new_curves$analysis <- "Filtered grouped: 39 islands"
comparison <- rbind(old_summary,new_summary)
comparison_curves <- rbind(old_curves,new_curves)
write.csv(comparison,"area-estimation/results/bayesian-filtered/comparison-with-previous.csv",row.names=FALSE)
base <- comparison[comparison$prior == "log_uniform_10_1million", ]
base$interval_width <- base$upper - base$lower
base$interval_fold <- base$upper / base$lower
write.csv(base,"area-estimation/results/bayesian-filtered/baseline-comparison.csv",row.names=FALSE)

# All distributions are normalized; no richness draws were removed.
for (name in unique(curves$model)) {
  for (pr in unique(curves$prior)) {
    zc <- curves[curves$model==name & curves$prior==pr, ]
    stopifnot(abs(sum(zc$density_log10_area)*.01-1)<1e-8)
  }
}
stopifnot(sum(frequency$n_maps)==2000,abs(sum(frequency$weight)-1)<1e-10)

# Plots ----
cols <- c("Previous simple: 61 islands"="#0072B2", "Previous flexible: 61 islands"="#D55E00",
          "Filtered grouped: 39 islands"="#009E73")
theme_set(theme_classic(base_size=12))
zplot <- comparison_curves[comparison_curves$prior == "log_uniform_10_1million" &
  comparison_curves$area_km2 >= 10 & comparison_curves$area_km2 <= 1e6, ]
p1 <- ggplot(zplot,aes(area_km2,density_log10_area,colour=analysis,fill=analysis)) +
  geom_area(alpha=.10,position="identity") + geom_line(linewidth=.85) +
  scale_x_log10(breaks=c(10,100,1000,1e4,1e5,1e6),labels=c("10","100","1k","10k","100k","1m")) +
  scale_colour_manual(values=cols) + scale_fill_manual(values=cols) +
  labs(x=expression("Effective area (km"^2*")"),y="Density per log10 area",colour=NULL,fill=NULL,
    caption="Same 9 Ma richness distribution and log-uniform area prior (10-1,000,000 km²).\nPrevious models: oceanic assumption; new model: non-Caribbean group, no geological-type term.") +
  theme(legend.position="bottom",legend.direction="vertical",plot.caption=element_text(hjust=0,size=9))
p2 <- ggplot(base,aes(median,analysis,colour=analysis)) +
  geom_segment(aes(x=lower,xend=upper,yend=analysis),linewidth=.8) + geom_point(size=3) +
  scale_x_log10() + scale_colour_manual(values=cols) +
  labs(x=expression("Median and 95% credible interval (km"^2*")"),y=NULL) + theme(legend.position="none")
prior_names <- c("log_uniform_10_1million","log_uniform_100_100k","log_uniform_1_10million","lognormal_10k","uniform_area_10_1million")
comparison$y <- match(comparison$prior,prior_names) +
  c(-.2,0,.2)[match(comparison$analysis,names(cols))]
p3 <- ggplot(comparison,aes(median,y,colour=analysis)) +
  geom_segment(aes(x=lower,xend=upper,yend=y),linewidth=.6) + geom_point(size=2) +
  scale_x_log10() + scale_colour_manual(values=cols) +
  scale_y_continuous(breaks=1:5,labels=c("Log-uniform: 10-1m","Log-uniform: 100-100k","Log-uniform: 1-10m","Lognormal: median 10k","Area-uniform: 10-1m")) +
  labs(x=expression("Median and 95% credible interval (km"^2*")"),y=NULL,colour=NULL) +
  theme(legend.position="bottom",legend.direction="vertical")
figures <- list(area_comparison=p1,area_intervals=p2,prior_comparison=p3)
for (name in names(figures)) {
  ggsave(paste0("area-estimation/plots/bayesian-filtered/",name,".pdf"),figures[[name]],width=10,height=6.5)
  ggsave(paste0("area-estimation/plots/bayesian-filtered/",name,".png"),figures[[name]],width=10,height=6.5,dpi=180)
}

# Short report ----
diag <- read.csv("area-estimation/results/bayesian-filtered/convergence.csv")
min_area <- min(d$area_km2[d$group=="Other islands"])
min_richness <- min(d$n_sp[d$group=="Other islands"])
p_low_richness <- sum(frequency$weight[frequency$n_species<min_richness])
znew <- new_curves[new_curves$prior=="log_uniform_10_1million", ]
p_low_area <- sum(znew$density_log10_area[znew$area_km2<min_area])*.01
report <- c("# Filtered non-Caribbean area reconstruction at 9 Ma", "",
"Used 39 islands (16 Caribbean, 23 other), retaining the agreed exclusions. The Bayesian lognormal regression has separate intercepts and slopes for Caribbean and non-Caribbean groups and a shared residual standard deviation. Inverse prediction uses only the non-Caribbean relationship. It does not include island-type or separate Indian Ocean effects.", "",
"The full 2,000-map combined richness distribution at exactly 9 Ma was retained: median 4, mean 5.5025, 95% interval 2–21, range 2–101. No map was discarded or reweighted. Regression coefficients and residual variation are integrated out before normalizing each richness-specific area distribution and combining by its map frequency. This is the same modular inverse-prediction procedure as before.", "",
"## Baseline comparison", "",
"Area prior: uniform on log10 area, 10–1,000,000 km². Previous estimates below use their oceanic-island setting; the new group model has no geological-type setting. Both data selection and model structure change, so their effects cannot be attributed separately.", "",
"| Analysis | Mean km² | Median km² | 95% interval km² |", "|---|---:|---:|---:|")
for(i in seq_len(nrow(base))) report <- c(report,sprintf("| %s | %.0f | %.0f | %.0f–%.0f |",base$analysis[i],base$mean[i],base$median[i],base$lower[i],base$upper[i]))
report <- c(report,"","## Area-prior sensitivity: new model", "",
"| Prior | Median km² | 95% interval km² |","|---|---:|---:|")
for(i in seq_len(nrow(new_summary))) report <- c(report,sprintf("| %s | %.0f | %.0f–%.0f |",new_summary$prior[i],new_summary$median[i],new_summary$lower[i],new_summary$upper[i]))
report <- c(report,"","## Checks and interpretation", "",
sprintf("Four independent chains for each of two regression-prior widths; 5,000 warmup and 25,000 retained iterations per chain. Max Rhat %.5f; minimum bulk/tail ESS %.0f / %.0f. All passed Rhat <1.01 and ESS >400. Area integration uses 4,000 posterior draws and a 0.01 log10-area grid; half-draw checks have maximum median/upper-limit change %.2f%% / %.2f%%. All distributions normalize to one.",max(diag$rhat),min(diag$ess_bulk),min(diag$ess_tail),max(abs(numerical$median_change_percent)),max(abs(numerical$upper_change_percent))),"",
"Centered predictor x=ln(area/10,000). Priors: Caribbean intercept Normal(ln10,2), Caribbean slope Normal(0.4,0.5), non-Caribbean intercept difference Normal(0,1.5), slope difference Normal(0,0.5), residual sigma half-Normal(0,1.5); values after commas are SDs. The wider alternative doubles every prior SD. Sampling log sigma includes the Jacobian. Proposal covariance is fixed from the posterior-mode Hessian and proposal scale is tuned only during warmup. Both chains and the wider-prior area estimates are saved.","",
sprintf("Important extrapolation: retained non-Caribbean islands have at least %d species and area at least %.0f km². %.1f%% of mapped richness values are below that minimum richness; %.1f%% of the baseline inferred area distribution lies below the minimum non-Caribbean calibration area. The Caribbean data help estimate shared residual scatter but do not replace missing low-richness non-Caribbean calibration islands.",min_richness,min_area,100*p_low_richness,100*p_low_area),"",
"Successive data-dependent exclusions reduce observed scatter and may narrow intervals without establishing better accuracy. Interpretation remains conditional on the excluded inventories being unsuitable calibration records. Coefficient intervals and area intervals do not include uncertainty about choosing exclusions or grouping. This is effective richness-equivalent area, not necessarily total exposed land area. Previous Bayesian results and original data remain unchanged.","",
"Reproduce from the project root: source('area-estimation/scripts/10-fit-filtered-bayesian.R'), then source('area-estimation/scripts/11-estimate-filtered-area.R'). Results, plots and logs are in the corresponding bayesian-filtered subdirectories.")
writeLines(report,"area-estimation/reports/bayesian-filtered-area.md")
print(base)
