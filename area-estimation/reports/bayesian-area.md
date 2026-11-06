# Bayesian area estimation at 9 Ma

The analysis uses every map-derived richness value, without truncation to the 95% interval. Combined richness has mean 5.5025, median 4, 95% interval 2–21 and full range 2–101 across 2,000 paired maps. All 61 calibration records were retained.

## Main estimates

Baseline area prior: uniform in log10 area between 10 and 1,000,000 km². These are provisional prior bounds. The three island types are alternative assumptions, not posterior classifications. The flexible model conditions on Indian Ocean; the simple model has no region term.

| Model | Island type | Median km² | Mean km² | 95% credible interval km² |
|---|---|---:|---:|---:|
| simple_prior1 | Continental fragment | 2,497 | 28,798 | 33–281,610 |
| simple_prior1 | Continental shelf | 212 | 4,393 | 12–28,313 |
| simple_prior1 | Oceanic | 1,660 | 22,379 | 26–215,860 |
| flexible_prior1 | Continental fragment | 1,318 | 11,296 | 23–85,497 |
| flexible_prior1 | Continental shelf | 71 | 4,696 | 11–25,377 |
| flexible_prior1 | Oceanic | 815 | 17,687 | 15–173,461 |

## Prior sensitivity: oceanic assumption

| Model | Area prior | Median km² | 95% interval km² |
|---|---|---:|---:|
| simple_prior1 | log_uniform_10_1million | 1660 | 26–215860 |
| simple_prior1 | log_uniform_100_100k | 1922 | 124–65018 |
| simple_prior1 | log_uniform_1_10million | 1628 | 17–298318 |
| simple_prior1 | lognormal_10k | 3516 | 140–114646 |
| simple_prior1 | uniform_area_10_1million | 82608 | 977–874901 |
| flexible_prior1 | log_uniform_10_1million | 815 | 15–173461 |
| flexible_prior1 | log_uniform_100_100k | 1228 | 114–57805 |
| flexible_prior1 | log_uniform_1_10million | 723 | 4–257649 |
| flexible_prior1 | lognormal_10k | 2594 | 87–101301 |
| flexible_prior1 | uniform_area_10_1million | 157782 | 1191–920801 |

The area prior strongly affects the answer. Log-uniform bounds changes have a smaller effect on the oceanic median than switching to an area-uniform prior, but greatly affect the tails and arithmetic mean. An area-uniform prior concentrates much more initial mass at large areas. The lognormal alternative is Normal(4,1) in log10 area, truncated to 10–1,000,000 km². No prior is designated uninformative.

Doubling regression-prior scales produces comparatively small changes for the oceanic case under the baseline area prior; see the regression-prior plot. This is limited sensitivity to the tested regression priors, not general prior robustness. Geological type is also consequential; particularly small shelf estimates approach the lower prior boundary.

## Computation and checks

Four independent chains per model/prior combination, each with 5,000 warmup and 25,000 retained Metropolis iterations. Maximum rank-normalized split Rhat = 1.00256; minimum bulk ESS = 2735; minimum tail ESS = 4538. All passed Rhat < 1.01 and ESS > 400.

Sampling used a multivariate normal random-walk proposal with covariance from the posterior-mode Hessian. Proposal scale was tuned only in warmup and fixed thereafter. This is MCMC, not a Laplace posterior approximation. Trace PDFs, diagnostics and complete draws are saved.

Area integration used 4,000 stratified retained regression draws and a 0.01-log10-area grid. Comparing alternate 2,000-draw subsets with all 4,000 gave maximum median change 4.62% and upper-limit change 1.93% across 60 model/type/prior combinations. Refinement from the earlier 2,000-draw, 0.02 grid gave maximum median/upper changes 4.20% / 5.96%. The analytic fixed-coefficient inverse-normal check had maximum quantile error 1.15%; all area distributions normalize to one.

The grid and posterior subset checks assess numerical stability, not biological validity. Conditional posterior predictive summaries for all modern islands are saved; they reuse training data and are not independent predictive validation. The earlier leave-one-out pilot remains the model-comparison evidence.

## Models and priors

Let x = ln(area/10,000 km²). Simple: ln(richness) ~ Normal(alpha + beta*x + island-type effects, sigma). Flexible: richness ~ negative binomial with log mean alpha + beta*x + type + x:type + region and size phi. Reference levels are continental fragment and Australasian.

Independent baseline coefficient priors: alpha ~ Normal(ln(10),2); beta ~ Normal(0.4,0.5); type and region offsets ~ Normal(0,1.5); slope differences ~ Normal(0,0.5). Standard deviations are given after the comma. Simple sigma ~ half-Normal(0,1.5); flexible ln(phi) ~ Normal(ln(2),1.5). The alternative doubles all these prior standard deviations. Sampling ln(sigma) includes its Jacobian. The prior on ln(phi) is directly specified in sampled coordinates. Priors permit negative slopes; no positive-slope constraint was imposed.

## How richness uncertainty was propagated

For each integer N with frequency w_N among the maps, compute L_N(A) = average over regression posterior draws of p(N | A, theta, type, region). For the lognormal model, using the density of ln(N) instead of N contributes a constant 1/N that cancels within each N. Then normalize p(A | N, calibration data) proportional to area_prior(A) L_N(A) separately for that N. Finally mix the normalized area distributions using w_N. Quadrature is in log10 area, with area priors transformed to that coordinate and endpoint weights included.

Thus maps are alternative histories, not 2,000 independent observations. Each N is normalized before combining, so the regression does not reweight the upstream mapping richness distribution. This is modular Bayesian inverse prediction conditional on the AIC-weighted mapping output, not a single fully joint Bayesian model of tree, richness and area. Regression uncertainty is marginalized before inverse normalization.

## Interpretation

These results do not support a precise historical island size. Report the whole distribution, the assumed geological type and area prior; the median is more representative than the mean for these strongly right-skewed distributions. Good MCMC convergence does not remove ecological or calibration uncertainty.

The target is a richness-equivalent area under the modern-island relationship. It is not automatically the physical sum of several disconnected Sunken Islands. The table mixes mostly individual islands with four archipelago records. Only four calibration records are Indian Ocean. Region, inventory source and type are partly confounded. Richness is positive in every calibration record; possible omission of zero-richness islands is not modeled. The modern relationship is assumed transferable to 9 Ma, and the two reconstructed clades are assumed comparable to the taxonomic scope of the modern table. Unrepresented extinct clades, habitat extent, island isolation, inventory error and geological area uncertainty are not modeled. Mauritius/Réunion historical richness records remain included.

## Files and rerunning

Run scripts 00, 02, 03 and 04 from the paleo-islands project root. Script 01 reproduces the earlier frequentist pilot and need not be rerun. Results are in area-estimation/results/bayesian; figures in area-estimation/plots/bayesian; logs in area-estimation/logs/bayesian. Scripts read the current source data in the repository. The coarse-grid CSV is an archived numerical-check reference from this run.

Plots: area_model_comparison (baseline area prior, three geological assumptions); area_prior_sensitivity (five priors); regression_prior_sensitivity (standard and doubled regression-prior widths). Densities are per unit log10 area, not per km²; intervals are equal-tailed 95% modular credible intervals. No model averaging between the simple and flexible regressions was imposed.

Diagnostic reference: https://mc-stan.org/docs/2_38/reference-manual/analysis.html ; inverse calibration: https://doi.org/10.1016/0167-9473(94)00014-A
