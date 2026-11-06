# Bayesian inverse prediction, preserving the upstream richness distribution.
library(ggplot2)
models <- readRDS("area-estimation/results/bayesian/regression-posteriors.rds")
frequency <- read.csv("area-estimation/results/inputs/richness-9ma-frequency.csv")
d <- read.csv("area-estimation/results/inputs/island-calibration.csv")
d$type <- factor(d$type); d$region <- factor(d$region)
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
  for (type in levels(d$type)) {
    new <- data.frame(x = log(area/10000),
      type = factor(type, levels = levels(d$type)),
      region = factor("Indian Ocean", levels = levels(d$region)))
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
write.csv(summaries,"area-estimation/results/bayesian/area-summary.csv",row.names=FALSE)
write.csv(curves,"area-estimation/results/bayesian/area-distributions.csv",row.names=FALSE)
write.csv(numerical,"area-estimation/results/bayesian/integration-check.csv",row.names=FALSE)
write.csv(predictive,"area-estimation/results/bayesian/posterior-predictive.csv",row.names=FALSE)

# Main comparison: both models, all three geological interpretations.
theme_set(theme_classic(base_size=11))
base <- curves[curves$prior=="log_uniform_10_1million" & grepl("prior1",curves$model),]
base <- base[base$area_km2 >= 10 & base$area_km2 <= 1e6, ]
base$model <- ifelse(startsWith(base$model,"simple"),"Simple power model","Flexible count model")
p1 <- ggplot(base,aes(area_km2,density_log10_area,colour=model,fill=model)) +
  geom_area(position="identity",alpha=.12) + geom_line(linewidth=.7) + facet_wrap(~type,ncol=1) +
  scale_x_log10(breaks=c(10,100,1000,1e4,1e5,1e6),labels=c("10","100","1k","10k","100k","1m")) +
  scale_colour_manual(values=c("#D55E00","#0072B2")) + scale_fill_manual(values=c("#D55E00","#0072B2")) +
  labs(x=expression("Area (km"^2*")"),y="Density per log10 area",colour=NULL,fill=NULL) + theme(legend.position="bottom")
s <- summaries[grepl("prior1",summaries$model),]
prior_order <- c("log_uniform_10_1million", "log_uniform_100_100k", "log_uniform_1_10million", "lognormal_10k", "uniform_area_10_1million")
s$y <- match(s$prior, prior_order) + ifelse(startsWith(s$model,"simple"), -.12, .12)
s$model <- ifelse(startsWith(s$model,"simple"), "Simple", "Flexible")
p2 <- ggplot(s,aes(median,y,colour=model)) +
  geom_segment(aes(x=lower,xend=upper,yend=y),linewidth=.6) + geom_point(size=2) +
  facet_wrap(~type,ncol=1) + scale_x_log10() +
  scale_y_continuous(breaks=1:5, labels=c("Log-uniform: 10-1m", "Log-uniform: 100-100k", "Log-uniform: 1-10m", "Lognormal: median 10k", "Area-uniform: 10-1m")) +
  scale_colour_manual(values=c(Flexible="#D55E00",Simple="#0072B2")) +
  labs(x=expression("Median area and 95% credible interval (km"^2*")"),y=NULL,colour=NULL) + theme(legend.position="bottom")
s <- summaries[summaries$prior=="log_uniform_10_1million",]
p3 <- ggplot(s,aes(median,model,colour=model)) +
  geom_segment(aes(x=lower,xend=upper,yend=model)) + geom_point(size=2) + facet_wrap(~type,ncol=1) +
  scale_x_log10() + labs(x=expression("Area and 95% interval (km"^2*")"),y=NULL) + theme(legend.position="none")
plots <- list(area_model_comparison=p1,area_prior_sensitivity=p2,regression_prior_sensitivity=p3)
for (name in names(plots)) {
  ggsave(paste0("area-estimation/plots/bayesian/",name,".pdf"),plots[[name]],width=9,height=8)
  ggsave(paste0("area-estimation/plots/bayesian/",name,".png"),plots[[name]],width=9,height=8,dpi=160)
}
cat("Elapsed seconds:",proc.time()[3]-started,"\n")
print(summaries[summaries$prior=="log_uniform_10_1million",])
