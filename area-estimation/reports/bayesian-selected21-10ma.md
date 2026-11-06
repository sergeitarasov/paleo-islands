# Sunken Islands effective area: selected 21-island calibration

Calibrated on 16 Caribbean records plus Réunion, Mauritius, Madagascar, New Caledonia and New Zealand (18 species). NZ area is provisionally 263,310 km2; see selected-21.md for source. No mainland-distance predictor is included. Model: ln(N) = alpha + beta ln(A/10,000) + epsilon. This is the pooled area-only model preferred in the requested subset comparison. The sampling distribution is lognormal; inverse prediction includes residual variation among islands, not only uncertainty in the fitted mean.

Used the full saved combined 2,000-map richness distribution at 10 Ma, including the tails. For each richness value, integrated the likelihood over regression draws, multiplied by the area prior and normalized. Then mixed these conditional area distributions with the unchanged map frequencies. This modular calculation preserves upstream richness probabilities rather than updating them with the island calibration data.

Regression priors: alpha ~ Normal(ln10, SD 2), beta ~ Normal(0.4, SD 0.5), sigma ~ half-Normal(SD 1.5). The wider-prior analysis doubles all these SDs. Four chains per prior, 5,000 warmup plus 25,000 retained iterations each; fixed-covariance random-walk Metropolis with scale adaptation during warmup only. Log-sigma target includes its Jacobian. Baseline area prior: log-uniform from 10 to 1,000,000 km2. No geology was used as a calibration target.

| Regression prior | Area prior | Mean km2 | Median km2 | 95% interval km2 |
|---|---|---:|---:|---:|
| simple_prior1 | log_uniform_10_1million | 22076 | 1666 | 29–214257 |
| simple_prior1 | log_uniform_100_100k | 8561 | 1893 | 126–64367 |
| simple_prior1 | log_uniform_1_10million | 49277 | 1635 | 18–309641 |
| simple_prior1 | lognormal_10k | 15740 | 3363 | 151–112559 |
| simple_prior1 | uniform_area_10_1million | 214246 | 90666 | 1007–891214 |
| simple_prior2 | log_uniform_10_1million | 22386 | 1668 | 28–218283 |
| simple_prior2 | log_uniform_100_100k | 8614 | 1902 | 126–64619 |
| simple_prior2 | log_uniform_1_10million | 50327 | 1636 | 18–319103 |
| simple_prior2 | lognormal_10k | 15935 | 3398 | 151–114174 |
| simple_prior2 | uniform_area_10_1million | 217486 | 94405 | 1044–892821 |

Maximum Rhat 1.00086; minimum bulk/tail ESS 7791 / 10620. Half-draw integration changes the median/upper limit by at most 1.10% / 0.80%. All area distributions normalize to one.

Interpretation: effective richness-equivalent area conditional on the selected calibration subset, regression, priors and upstream maps. This is not an independent measurement of emerged geological area. Model-selection uncertainty and systematic richness/inventory error are not propagated. Earlier analyses and inputs were left unchanged.
