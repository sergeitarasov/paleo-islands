# Sunken Islands effective area: selected 27-island calibration

Calibrated on 16 Caribbean records plus Réunion, Mauritius, Madagascar, New Caledonia and New Zealand (18 species), Borneo, Java, Sumatra, Sulawesi, New Guinea and Sri Lanka. NZ area is provisionally 263,310 km2; see selected-27.md for source. No mainland-distance predictor is included. Model: ln(N) = alpha + beta ln(A/10,000) + epsilon. This analysis uses the pooled area-only model to retain comparability with the 21-record estimate. The sampling distribution is lognormal; inverse prediction includes residual variation among islands, not only uncertainty in the fitted mean.

Used the full saved combined 2,000-map richness distribution at 10 Ma, including the tails. For each richness value, integrated the likelihood over regression draws, multiplied by the area prior and normalized. Then mixed these conditional area distributions with the unchanged map frequencies. This modular calculation preserves upstream richness probabilities rather than updating them with the island calibration data.

Regression priors: alpha ~ Normal(ln10, SD 2), beta ~ Normal(0.4, SD 0.5), sigma ~ half-Normal(SD 1.5). The wider-prior analysis doubles all these SDs. Four chains per prior, 5,000 warmup plus 25,000 retained iterations each; fixed-covariance random-walk Metropolis with scale adaptation during warmup only. Log-sigma target includes its Jacobian. Baseline area prior: log-uniform from 10 to 1,000,000 km2. No geology was used as a calibration target.

| Regression prior | Area prior | Mean km2 | Median km2 | 95% interval km2 |
|---|---|---:|---:|---:|
| simple_prior1 | log_uniform_10_1million | 10752 | 1322 | 41–80842 |
| simple_prior1 | log_uniform_100_100k | 6103 | 1465 | 127–48120 |
| simple_prior1 | log_uniform_1_10million | 14903 | 1313 | 36–84381 |
| simple_prior1 | lognormal_10k | 9112 | 2444 | 155–58186 |
| simple_prior1 | uniform_area_10_1million | 88095 | 17576 | 466–669230 |
| simple_prior2 | log_uniform_10_1million | 10789 | 1322 | 41–81244 |
| simple_prior2 | log_uniform_100_100k | 6115 | 1468 | 127–48201 |
| simple_prior2 | log_uniform_1_10million | 14950 | 1313 | 35–84838 |
| simple_prior2 | lognormal_10k | 9146 | 2454 | 154–58425 |
| simple_prior2 | uniform_area_10_1million | 88714 | 17914 | 470–670212 |

Maximum Rhat 1.00073; minimum bulk/tail ESS 8109 / 10532. Half-draw integration changes the median/upper limit by at most 0.60% / 0.44%. All area distributions normalize to one.

Interpretation: effective richness-equivalent area conditional on the selected calibration subset, regression, priors and upstream maps. This is not an independent measurement of emerged geological area. Model-selection uncertainty and systematic richness/inventory error are not propagated. Earlier analyses and inputs were left unchanged.

Within this 27-record subset, the separate Caribbean/rest slope model is preferred to the pooled model (AICc 64.80 versus 69.59; LOO log RMSE 0.790 versus 0.825). The area estimate above deliberately retains the pooled model for direct methodological comparison with the previous 21-record analysis; it should not be described as the best-supported model for this expanded subset. The larger pooled R2 after adding large islands is not by itself evidence of improved inverse-prediction accuracy.
