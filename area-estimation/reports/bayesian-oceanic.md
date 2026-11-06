# Sunken Islands effective area: oceanic-only calibration

Calibrated on all 31 oceanic islands, without richness-based exclusions. Model: ln(N) = alpha + beta ln(A/10,000) + epsilon. This is the area-only model selected by AICc in the oceanic comparison. The sampling distribution is lognormal; inverse prediction includes residual variation among islands, not only uncertainty in the fitted mean.

Used the full saved combined 2,000-map richness distribution at 9 Ma, including the tails. For each richness value, integrated the likelihood over regression draws, multiplied by the area prior and normalized. Then mixed these conditional area distributions with the unchanged map frequencies. This modular calculation preserves upstream richness probabilities rather than updating them with the island calibration data.

Regression priors: alpha ~ Normal(ln10, SD 2), beta ~ Normal(0.4, SD 0.5), sigma ~ half-Normal(SD 1.5). The wider-prior analysis doubles all these SDs. Four chains per prior, 5,000 warmup plus 25,000 retained iterations each; fixed-covariance random-walk Metropolis with scale adaptation during warmup only. Log-sigma target includes its Jacobian. Baseline area prior: log-uniform from 10 to 1,000,000 km2. No geology was used as a calibration target.

| Regression prior | Area prior | Mean km2 | Median km2 | 95% interval km2 |
|---|---|---:|---:|---:|
| simple_prior1 | log_uniform_10_1million | 12328 | 1528 | 63–93469 |
| simple_prior1 | log_uniform_100_100k | 6661 | 1585 | 133–53256 |
| simple_prior1 | log_uniform_1_10million | 18501 | 1527 | 61–97431 |
| simple_prior1 | lognormal_10k | 9901 | 2308 | 164–64703 |
| simple_prior1 | uniform_area_10_1million | 57531 | 8213 | 306–535234 |
| simple_prior2 | log_uniform_10_1million | 12288 | 1519 | 63–93142 |
| simple_prior2 | log_uniform_100_100k | 6638 | 1577 | 133–53111 |
| simple_prior2 | log_uniform_1_10million | 18450 | 1518 | 61–97112 |
| simple_prior2 | lognormal_10k | 9868 | 2301 | 163–64472 |
| simple_prior2 | uniform_area_10_1million | 57202 | 8187 | 305–531626 |

Maximum Rhat 1.00088; minimum bulk/tail ESS 8008 / 10782. Half-draw integration changes the median/upper limit by at most 1.20% / 2.32%. All area distributions normalize to one.

Interpretation: effective richness-equivalent area conditional on the selected oceanic classification, regression, priors and upstream maps. This is not an independent measurement of emerged geological area. Model-selection uncertainty and systematic richness/inventory error are not propagated. Earlier analyses and inputs were left unchanged.
