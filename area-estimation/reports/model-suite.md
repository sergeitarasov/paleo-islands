# Regression model comparison on the unfiltered island data

All 61 islands retained. Sixteen candidates were fitted: seven lognormal and seven negative-binomial regressions, plus two lognormal mixed-effects models. No age predictor was used because the available ages are sparse and heterogeneous. No new exclusions, mixture components or Sunken Islands area estimates were introduced.

## What was tested

- Area only; area + mainland isolation; area + three regions + three island types; these additive effects plus isolation; the previous area-by-Caribbean/rest interaction plus isolation.
- Natural cubic spline in log area with three degrees of freedom, alone and with additive region, type and isolation. Spline knots are estimated using each training fold only.
- Random geographic-group intercept with area, type and isolation; and the same with fixed region effects. These mixed models use ML, not REML, for likelihood comparisons.

## Validation

Each candidate predicts every island twice: once with that island withheld, and once with its entire geographic block withheld. Fifteen blocks comprise nine multi-island groups and six singletons; the mapping is saved in geographic-blocks.csv. Random effects for unseen groups are integrated as zero-mean normals, adding group variance to residual variance. Singleton groups are also treated as unseen during island-level validation.

Primary score is mean negative log probability of the observed count. For lognormal models, probabilities are integrated over the integer bin [N-0.5, N+0.5]; negative-binomial models use their probability mass. This avoids comparing a continuous density directly with a discrete probability. Zero is a permitted prediction; models are not conditioned on richness being positive. Secondary log RMSE uses lognormal medians or negative-binomial count medians (zero medians floored at 0.5 for the log metric). AICc is only comparable within likelihood families; it is not used to rank count and continuous fits together.

| Model | Holdout | Negative log score | Log RMSE |
|---|---|---:|---:|
| Area x Caribbean/rest + isolation | Lognormal | Hold out geographic group | 3.267 | 0.902 |
| Area + isolation | Lognormal | Hold out geographic group | 3.298 | 0.942 |
| Curved area | Lognormal | Hold out geographic group | 3.321 | 0.944 |
| Area x Caribbean/rest + isolation | Negative binomial | Hold out geographic group | 3.327 | 0.918 |
| Area + type + isolation + group RE | Mixed lognormal | Hold out geographic group | 3.346 | 0.955 |
| Area | Lognormal | Hold out geographic group | 3.365 | 1.004 |
| Curved area + region + type + isolation | Lognormal | Hold out geographic group | 3.390 | 0.948 |
| Area + region + type + isolation | Lognormal | Hold out geographic group | 3.405 | 0.979 |
| Area + isolation | Negative binomial | Hold out geographic group | 3.415 | 0.950 |
| Area + region + type + isolation + group RE | Mixed lognormal | Hold out geographic group | 3.453 | 1.015 |
| Area + region + type | Lognormal | Hold out geographic group | 3.460 | 1.015 |
| Area + region + type + isolation | Negative binomial | Hold out geographic group | 3.554 | 1.035 |
| Area | Negative binomial | Hold out geographic group | 3.561 | 1.034 |
| Curved area | Negative binomial | Hold out geographic group | 3.577 | 1.028 |
| Area + region + type | Negative binomial | Hold out geographic group | 3.608 | 1.059 |
| Curved area + region + type + isolation | Negative binomial | Hold out geographic group | 3.643 | 1.027 |
| Area x Caribbean/rest + isolation | Lognormal | Hold out island | 3.220 | 0.865 |
| Area + isolation | Lognormal | Hold out island | 3.270 | 0.915 |
| Curved area | Lognormal | Hold out island | 3.277 | 0.906 |
| Area x Caribbean/rest + isolation | Negative binomial | Hold out island | 3.282 | 0.892 |
| Area | Lognormal | Hold out island | 3.333 | 0.969 |
| Area + type + isolation + group RE | Mixed lognormal | Hold out island | 3.350 | 0.934 |
| Area + region + type + isolation | Lognormal | Hold out island | 3.363 | 0.945 |
| Area + isolation | Negative binomial | Hold out island | 3.370 | 0.908 |
| Curved area + region + type + isolation | Lognormal | Hold out island | 3.372 | 0.935 |
| Area + region + type + isolation + group RE | Mixed lognormal | Hold out island | 3.382 | 0.943 |
| Area + region + type | Lognormal | Hold out island | 3.386 | 0.959 |
| Area + region + type + isolation | Negative binomial | Hold out island | 3.451 | 0.947 |
| Area + region + type | Negative binomial | Hold out island | 3.485 | 0.968 |
| Curved area | Negative binomial | Hold out island | 3.495 | 0.953 |
| Curved area + region + type + isolation | Negative binomial | Hold out island | 3.509 | 0.933 |
| Area | Negative binomial | Hold out island | 3.522 | 1.005 |

## Main finding

The previous lognormal area x Caribbean/rest + mainland-distance model ranks first for both held-out count probability and log RMSE under both validation schemes. Additional region/type effects, negative-binomial errors, splines and group random effects did not improve these scores. The next contenders are simpler lognormal area + isolation and the simple area spline. Ranking differences are uncertain: the previous model gains only 1.83 total blocked log-score units over area + isolation (approximate block-level SE 6.11), and 3.63 over the matching negative-binomial model (SE 3.30). This is not decisive evidence for one unique best model.

All 16 full-data fits and 1,216 validation fits completed, with 61 predictions per candidate per holdout scheme. One mixed-model blocked fit (region-adjusted, withholding Moluccas) had a singular random effect; this diagnostic is retained rather than dropping the fold. Full-data mixed-model random-intercept SDs were about 0.326 and 0.319 log-richness units.

## Inverse area validation

Six representative candidates were tested by withholding whole geographic blocks, supplying the withheld island richness and its other covariates, and estimating area. Every inverse calculation uses the same log-uniform area prior, 10–1,000,000 km2. These are approximate screening intervals: 800 parameter draws per training fit. Lognormal fixed-effect draws use the reference posterior proportional to 1/sigma2; negative-binomial draws use a joint-normal coefficient approximation and an independent log-dispersion approximation. Mixed-model fixed coefficients are sampled approximately, but variance-component uncertainty is not propagated. A new-group random intercept is integrated. These differ from the proper priors used for the final Bayesian Sunken Islands analysis.

| Model | Coverage of 95% interval | Median absolute area-error factor | Median upper/lower interval ratio |
|---|---:|---:|---:|
| lognormal:area | 93.4% | 2.54 | 1135 |
| lognormal:isolation | 93.4% | 2.79 | 770 |
| lognormal:previous_best | 95.1% | 2.83 | 976 |
| lognormal:curve | 93.4% | 3.02 | 1528 |
| NB:previous_best | 95.1% | 3.51 | 870 |
| mixed:block | 95.1% | 2.70 | 920 |
| prior_only | 98.4% | 5.18 | 56234 |

Area-only regression has the smallest median absolute inverse-area error among these six candidates (factor 2.54 versus 2.83 for the previous richness-prediction winner). Its interval coverage is 57/61 (93.4%); the previous model covers 58/61 (95.1%). Mixed effects cover 58/61 with median error factor 2.70. The ranking depends on the prediction target, so the richness-score winner should not automatically be called the best area-reconstruction model. Median interval upper/lower ratios are about 770–1,528 across candidates: uncertainty remains substantial.

Coverage is descriptive for these 61 observations, not proof of calibrated 95% intervals for Sunken Islands. Large interval widths can deliver high coverage without precise area estimation. The prior-only benchmark separates model information from coverage caused by a broad prior. Candidate selection is not nested within an outer validation loop; all rankings remain exploratory. Inventory differences, classification uncertainty and palaeogeographic uncertainty are not propagated.

## Geographic grouping provenance

The workbook had no archipelago field. The separate grouping is an analyst-defined broad geographic blocking scheme, not a newly verified formal archipelago database. It uses Philippines, Greater Sunda, Lesser Sunda, Moluccas, Greater Antilles, Lesser Antilles, Bahamas, Virgin Islands and Mascarenes, with singleton blocks for Langkawi, Tioman, New Caledonia, New Guinea, Sri Lanka and Madagascar. Wetar is assigned to Moluccas; alternative grouping definitions can change mixed-model and blocked-validation results. Geographic association is not proof of independent evolutionary histories.

Indonesian group context: [UK Permanent Committee on Geographical Names, Indonesia toponymic factfile](https://assets.publishing.service.gov.uk/media/67852932c6428e0131881700/Indonesia_Toponymic_Factfile.pdf). Existing Matthews et al. data were consulted for available archipelago names, but not blindly name-matched: its Flores entry is the Azores island, not Indonesian Flores. The original workbook was not modified.

## Reproduce

From the paleo-islands project root run scripts 23-model-comparison.R, 24-inverse-validation.R and 25-report-models.R in order. The geographic-blocks.csv input is versionable and explicit. Results and figures are in model-suite subdirectories. The source workbook is read only. Existing Bayesian area estimates remain unchanged.
