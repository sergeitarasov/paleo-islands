# Bayesian inverse prediction, preserving the upstream richness distribution.
library(ggplot2)
models <- readRDS("area-estimation/results/bayesian-mixture-high/regression-posteriors.rds")
frequency <- read.csv("area-estimation/results/inputs/richness-9ma-frequency.csv")
d <- read.csv("area-estimation/results/bayesian-mixture-high/calibration-data.csv")
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
  for (type in "Higher-richness component") {
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
write.csv(summaries,"area-estimation/results/bayesian-mixture-high/area-summary.csv",row.names=FALSE)
write.csv(curves,"area-estimation/results/bayesian-mixture-high/area-distributions.csv",row.names=FALSE)
write.csv(numerical,"area-estimation/results/bayesian-mixture-high/integration-check.csv",row.names=FALSE)
write.csv(predictive,"area-estimation/results/bayesian-mixture-high/posterior-predictive.csv",row.names=FALSE)


# Comparison with the non-Caribbean analysis ----
previous <- read.csv("area-estimation/results/bayesian-filtered/area-summary.csv")
previous_curves <- read.csv("area-estimation/results/bayesian-filtered/area-distributions.csv")
previous <- previous[previous$model=="simple_prior1", ]
previous_curves <- previous_curves[previous_curves$model=="simple_prior1", ]
previous$analysis <- "Filtered non-Caribbean relationship"
previous_curves$analysis <- "Filtered non-Caribbean relationship"
current <- summaries[summaries$model=="simple_prior1", ]
current_curves <- curves[curves$model=="simple_prior1", ]
current$analysis <- "Selected higher-richness islands"
current_curves$analysis <- "Selected higher-richness islands"
comparison <- rbind(previous,current)
comparison_curves <- rbind(previous_curves,current_curves)
write.csv(comparison,"area-estimation/results/bayesian-mixture-high/comparison.csv",row.names=FALSE)
base <- comparison[comparison$prior=="log_uniform_10_1million", ]
write.csv(base,"area-estimation/results/bayesian-mixture-high/baseline-comparison.csv",row.names=FALSE)
for (m in unique(curves$model)) {
  for (pr in unique(curves$prior)) {
    zc <- curves[curves$model==m & curves$prior==pr, ]
    stopifnot(abs(sum(zc$density_log10_area)*.01-1)<1e-8)
  }
}

# Density and interval comparisons ----
theme_set(theme_classic(base_size=12))
colors <- c("Filtered non-Caribbean relationship"="#009E73", "Selected higher-richness islands"="#CC79A7")
zplot <- comparison_curves[comparison_curves$prior=="log_uniform_10_1million" &
  comparison_curves$area_km2>=10 & comparison_curves$area_km2<=1e6, ]
p1 <- ggplot(zplot,aes(area_km2,density_log10_area,colour=analysis,fill=analysis)) +
  geom_area(alpha=.12,position="identity") + geom_line(linewidth=.8) +
  scale_x_log10(breaks=c(10,100,1000,1e4,1e5,1e6),labels=c("10","100","1k","10k","100k","1m")) +
  scale_colour_manual(values=colors) + scale_fill_manual(values=colors) +
  labs(x=expression("Effective area (km"^2*")"),y="Density per log10 area",colour=NULL,fill=NULL,
    caption="Exploratory: mixture membership was selected using the same richness data.") +
  theme(legend.position="bottom",legend.direction="vertical",plot.caption=element_text(hjust=0))
prior_names <- c("log_uniform_10_1million","log_uniform_100_100k","log_uniform_1_10million","lognormal_10k","uniform_area_10_1million")
comparison$y <- match(comparison$prior,prior_names)+ifelse(comparison$analysis=="Selected higher-richness islands",.12,-.12)
p2 <- ggplot(comparison,aes(median,y,colour=analysis)) +
  geom_segment(aes(x=lower,xend=upper,yend=y),linewidth=.7) + geom_point(size=2.5) +
  scale_x_log10() + scale_colour_manual(values=colors) +
  scale_y_continuous(breaks=1:5,labels=c("Log-uniform: 10-1m","Log-uniform: 100-100k","Log-uniform: 1-10m","Lognormal: median 10k","Area-uniform: 10-1m")) +
  labs(x=expression("Median and 95% interval (km"^2*")"),y=NULL,colour=NULL) +
  theme(legend.position="bottom",legend.direction="vertical")
figures <- list(area_comparison=p1,prior_comparison=p2)
for (name in names(figures)) {
  ggsave(paste0("area-estimation/plots/bayesian-mixture-high/",name,".pdf"),figures[[name]],width=10,height=6.5)
  ggsave(paste0("area-estimation/plots/bayesian-mixture-high/",name,".png"),figures[[name]],width=10,height=6.5,dpi=180)
}

# Report ----
diag <- read.csv("area-estimation/results/bayesian-mixture-high/convergence.csv")
report <- c("# Exploratory area estimate from higher-richness mixture islands", "",
"Selected all 23 pink/red points assigned to the higher-richness component in the previous mixture fit: 7 Caribbean and 16 other islands. Refit a single Bayesian power regression to these observations, ignoring region because the requested subset defines this exploratory calibration. This is not a Bayesian mixture fit: hard assignments are treated as fixed, and membership/selection uncertainty is not propagated.", "",
"Used the unchanged complete 2,000-map combined richness distribution at 9 Ma (mean 5.5025, median 4, 95% interval 2–21, full range 2–101). All richness values and their frequencies were retained. The same modular inverse calibration marginalizes regression uncertainty, normalizes the area distribution conditional on each richness value, then combines by map frequencies.", "",
"## Baseline comparison", "",
"Area prior: uniform on log10 area between 10 and 1,000,000 km².", "",
"| Analysis | Mean km² | Median km² | 95% interval km² |", "|---|---:|---:|---:|")
for(i in seq_len(nrow(base))) report <- c(report,sprintf("| %s | %.0f | %.0f | %.0f–%.0f |",base$analysis[i],base$mean[i],base$median[i],base$lower[i],base$upper[i]))
report <- c(report,"","## Area-prior sensitivity", "", "| Prior | Median km² | 95% interval km² |", "|---|---:|---:|")
for(i in seq_len(nrow(current))) report <- c(report,sprintf("| %s | %.0f | %.0f–%.0f |",current$prior[i],current$median[i],current$lower[i],current$upper[i]))
report <- c(report,"",
sprintf("Four chains for each of two regression-prior widths, 5,000 warmup + 25,000 retained iterations each. Maximum Rhat %.5f; minimum bulk/tail ESS %.0f / %.0f. All passed Rhat <1.01 and ESS >400. Integration uses 4,000 draws and a 0.01 log10-area grid; half-draw changes in median/upper limit are at most %.2f%% / %.2f%%. All area distributions normalize to one.",max(diag$rhat),min(diag$ess_bulk),min(diag$ess_tail),max(abs(numerical$median_change_percent)),max(abs(numerical$upper_change_percent))),"",
"Model: ln(N)=alpha+beta*ln(A/10,000)+epsilon. Priors alpha ~ Normal(ln10,2), beta ~ Normal(0.4,0.5), sigma ~ half-Normal(0,1.5), with SDs after commas. The wider-prior fit doubles these SDs. Log-sigma sampling includes its Jacobian. MCMC proposal covariance comes from posterior-mode curvature; proposal scale is tuned only during warmup.","",
"Any increased precision is conditional on selecting the higher-richness pattern using the same response data. Hard selection removes low-richness residuals and ignores uncertain memberships, so it can understate population scatter and uncertainty. This experiment does not validate the selected subset as the correct analogue for Sunken Islands. The mixture itself predicted less well than the known-region model. Results are exploratory effective richness-equivalent area, not a primary geological estimate. No original data, prior analyses or stochastic maps were changed.","",
"Scripts: 13-fit-mixture-high-bayesian.R and 14-estimate-mixture-high-area.R. Results/plots/logs are in bayesian-mixture-high subdirectories. calibration-data.csv lists all included islands. Baseline and wider-regression-prior results remain separate.")
writeLines(report,"area-estimation/reports/bayesian-mixture-high.md")
print(base)
