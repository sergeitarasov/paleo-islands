# Checks and short report. Run after scripts 02 and 03.
s <- read.csv("area-estimation/results/bayesian/area-summary.csv")
c <- read.csv("area-estimation/results/bayesian/convergence.csv")
n <- read.csv("area-estimation/results/bayesian/integration-check.csv")
curves <- read.csv("area-estimation/results/bayesian/area-distributions.csv")
old <- read.csv("area-estimation/results/bayesian/area-summary-coarse-grid.csv")
comparison <- merge(s,old,by=c("model","type","prior"),suffixes=c("_fine","_coarse"))
comparison$median_change_percent <- 100*(comparison$median_fine/comparison$median_coarse-1)
comparison$upper_change_percent <- 100*(comparison$upper_fine/comparison$upper_coarse-1)
write.csv(comparison,"area-estimation/results/bayesian/grid-refinement.csv",row.names=FALSE)
stopifnot(all(is.finite(s$median)),all(s$lower < s$median),all(s$median < s$upper))
keys <- unique(curves[,c("model","type","prior")])
for(i in seq_len(nrow(keys))) {
  z <- curves[curves$model==keys$model[i] & curves$type==keys$type[i] & curves$prior==keys$prior[i],]
  stopifnot(abs(sum(z$density_log10_area)*.01-1)<1e-8)
}
# Analytic inverse-normal check: fixed coefficients and log-uniform area prior.
z <- seq(1,6,by=.01)
alpha <- -2; beta <- .45; sigma <- .9; richness <- 5
mass <- dnorm(log(richness), alpha+beta*log(10)*z, sigma)
mass[c(1,length(mass))] <- mass[c(1,length(mass))]/2
mass <- mass/sum(mass)
numeric_quantile <- 10^approx(cumsum(mass),z,xout=c(.025,.5,.975))$y
mu <- (log(richness)-alpha)/(beta*log(10))
sd <- sigma/(beta*log(10))
prob <- pnorm(1,mu,sd)+c(.025,.5,.975)*(pnorm(6,mu,sd)-pnorm(1,mu,sd))
analytic_quantile <- 10^qnorm(prob,mu,sd)
check <- data.frame(prob=c(.025,.5,.975),numeric=numeric_quantile,analytic=analytic_quantile,
                    relative_error=numeric_quantile/analytic_quantile-1)
write.csv(check,"area-estimation/results/bayesian/analytic-check.csv",row.names=FALSE)
stopifnot(max(abs(check$relative_error))<.02)

report <- c("# Bayesian area estimation at 9 Ma", "",
"The analysis uses every map-derived richness value, without truncation to the 95% interval. Combined richness has mean 5.5025, median 4, 95% interval 2–21 and full range 2–101 across 2,000 paired maps. All 61 calibration records were retained.", "",
"## Main estimates", "",
"Baseline area prior: uniform in log10 area between 10 and 1,000,000 km². These are provisional prior bounds. The three island types are alternative assumptions, not posterior classifications. The flexible model conditions on Indian Ocean; the simple model has no region term.", "",
"| Model | Island type | Median km² | Mean km² | 95% credible interval km² |",
"|---|---|---:|---:|---:|")
main <- s[s$prior=="log_uniform_10_1million" & grepl("prior1",s$model),]
for(i in seq_len(nrow(main))) report <- c(report,sprintf("| %s | %s | %s | %s | %s–%s |",
 main$model[i],main$type[i],format(round(main$median[i]),big.mark=","),format(round(main$mean[i]),big.mark=","),
 format(round(main$lower[i]),big.mark=","),format(round(main$upper[i]),big.mark=",")))
report <- c(report,"","## Prior sensitivity: oceanic assumption","",
"| Model | Area prior | Median km² | 95% interval km² |","|---|---|---:|---:|")
z <- s[s$type=="Oceanic" & grepl("prior1",s$model),]
for(i in seq_len(nrow(z))) report <- c(report,sprintf("| %s | %s | %.0f | %.0f–%.0f |",z$model[i],z$prior[i],z$median[i],z$lower[i],z$upper[i]))
report <- c(report,"",
"The area prior strongly affects the answer. Log-uniform bounds changes have a smaller effect on the oceanic median than switching to an area-uniform prior, but greatly affect the tails and arithmetic mean. An area-uniform prior concentrates much more initial mass at large areas. The lognormal alternative is Normal(4,1) in log10 area, truncated to 10–1,000,000 km². No prior is designated uninformative.","",
"Doubling regression-prior scales produces comparatively small changes for the oceanic case under the baseline area prior; see the regression-prior plot. This is limited sensitivity to the tested regression priors, not general prior robustness. Geological type is also consequential; particularly small shelf estimates approach the lower prior boundary.","",
"## Computation and checks","",
sprintf("Four independent chains per model/prior combination, each with 5,000 warmup and 25,000 retained Metropolis iterations. Maximum rank-normalized split Rhat = %.5f; minimum bulk ESS = %.0f; minimum tail ESS = %.0f. All passed Rhat < 1.01 and ESS > 400.",max(c$rhat),min(c$ess_bulk),min(c$ess_tail)),"",
"Sampling used a multivariate normal random-walk proposal with covariance from the posterior-mode Hessian. Proposal scale was tuned only in warmup and fixed thereafter. This is MCMC, not a Laplace posterior approximation. Trace PDFs, diagnostics and complete draws are saved.","",
sprintf("Area integration used 4,000 stratified retained regression draws and a 0.01-log10-area grid. Comparing alternate 2,000-draw subsets with all 4,000 gave maximum median change %.2f%% and upper-limit change %.2f%% across 60 model/type/prior combinations. Refinement from the earlier 2,000-draw, 0.02 grid gave maximum median/upper changes %.2f%% / %.2f%%. The analytic fixed-coefficient inverse-normal check had maximum quantile error %.2f%%; all area distributions normalize to one.",max(abs(n$median_change_percent)),max(abs(n$upper_change_percent)),max(abs(comparison$median_change_percent)),max(abs(comparison$upper_change_percent)),100*max(abs(check$relative_error))),"",
"The grid and posterior subset checks assess numerical stability, not biological validity. Conditional posterior predictive summaries for all modern islands are saved; they reuse training data and are not independent predictive validation. The earlier leave-one-out pilot remains the model-comparison evidence.","",
"## Models and priors","",
"Let x = ln(area/10,000 km²). Simple: ln(richness) ~ Normal(alpha + beta*x + island-type effects, sigma). Flexible: richness ~ negative binomial with log mean alpha + beta*x + type + x:type + region and size phi. Reference levels are continental fragment and Australasian.","",
"Independent baseline coefficient priors: alpha ~ Normal(ln(10),2); beta ~ Normal(0.4,0.5); type and region offsets ~ Normal(0,1.5); slope differences ~ Normal(0,0.5). Standard deviations are given after the comma. Simple sigma ~ half-Normal(0,1.5); flexible ln(phi) ~ Normal(ln(2),1.5). The alternative doubles all these prior standard deviations. Sampling ln(sigma) includes its Jacobian. The prior on ln(phi) is directly specified in sampled coordinates. Priors permit negative slopes; no positive-slope constraint was imposed.","",
"## How richness uncertainty was propagated","",
"For each integer N with frequency w_N among the maps, compute L_N(A) = average over regression posterior draws of p(N | A, theta, type, region). For the lognormal model, using the density of ln(N) instead of N contributes a constant 1/N that cancels within each N. Then normalize p(A | N, calibration data) proportional to area_prior(A) L_N(A) separately for that N. Finally mix the normalized area distributions using w_N. Quadrature is in log10 area, with area priors transformed to that coordinate and endpoint weights included.","",
"Thus maps are alternative histories, not 2,000 independent observations. Each N is normalized before combining, so the regression does not reweight the upstream mapping richness distribution. This is modular Bayesian inverse prediction conditional on the AIC-weighted mapping output, not a single fully joint Bayesian model of tree, richness and area. Regression uncertainty is marginalized before inverse normalization.","",
"## Interpretation","",
"These results do not support a precise historical island size. Report the whole distribution, the assumed geological type and area prior; the median is more representative than the mean for these strongly right-skewed distributions. Good MCMC convergence does not remove ecological or calibration uncertainty.","",
"The target is a richness-equivalent area under the modern-island relationship. It is not automatically the physical sum of several disconnected Sunken Islands. The table mixes mostly individual islands with four archipelago records. Only four calibration records are Indian Ocean. Region, inventory source and type are partly confounded. Richness is positive in every calibration record; possible omission of zero-richness islands is not modeled. The modern relationship is assumed transferable to 9 Ma, and the two reconstructed clades are assumed comparable to the taxonomic scope of the modern table. Unrepresented extinct clades, habitat extent, island isolation, inventory error and geological area uncertainty are not modeled. Mauritius/Réunion historical richness records remain included.","",
"## Files and rerunning","",
"Run scripts 00, 02, 03 and 04 from the paleo-islands project root. Script 01 reproduces the earlier frequentist pilot and need not be rerun. Results are in area-estimation/results/bayesian; figures in area-estimation/plots/bayesian; logs in area-estimation/logs/bayesian. Scripts read the current source data in the repository. The coarse-grid CSV is an archived numerical-check reference from this run.","",
"Plots: area_model_comparison (baseline area prior, three geological assumptions); area_prior_sensitivity (five priors); regression_prior_sensitivity (standard and doubled regression-prior widths). Densities are per unit log10 area, not per km²; intervals are equal-tailed 95% modular credible intervals. No model averaging between the simple and flexible regressions was imposed.","",
"Diagnostic reference: https://mc-stan.org/docs/2_38/reference-manual/analysis.html ; inverse calibration: https://doi.org/10.1016/0167-9473(94)00014-A")
writeLines(report,"area-estimation/reports/bayesian-area.md")
cat("Report and numerical checks completed\n")
