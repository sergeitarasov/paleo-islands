# Richness predicted from geological areas at 10 Ma

Scotese: approximately 32,000 km²; Straume: approximately 114,000 km², read from the supplied figure. Geological areas are treated as fixed scenarios without measurement errors. All emerged area is assumed to function as effective area for the calibrated fauna. No age, isolation or carrying-capacity mechanism is included.

Forward posterior predictive calculation using the saved 27- and 22-island pooled and separate-slope Bayesian regressions. Regional predictions use the non-Caribbean group. Each model uses all 100,000 retained posterior draws. Richness has a continuous lognormal sampling distribution: mixture CDF quantiles are obtained numerically; the predictive mean averages exp(mu+sigma²/2). Fractional summaries are expected and do not represent fractional observed species. Intervals include both regression uncertainty and residual variation among islands, not geological uncertainty.

| Dataset | Model | Geology | Mean | Median | 95% predictive interval |
|---|---|---|---:|---:|---:|
| 27 islands | simple_prior1 | Scotese | 30.7 | 21.0 | 3.8–115.7 |
| 27 islands | simple_prior1 | Straume | 63.4 | 43.0 | 7.6–241.6 |
| 27 islands | separate_prior1 | Scotese | 31.7 | 23.6 | 5.2–106.8 |
| 27 islands | separate_prior1 | Straume | 74.3 | 55.5 | 12.4–248.0 |
| 22 islands | simple_prior1 | Scotese | 24.4 | 15.8 | 2.6–97.3 |
| 22 islands | simple_prior1 | Straume | 46.7 | 29.7 | 4.6–190.2 |
| 22 islands | separate_prior1 | Scotese | 32.0 | 21.5 | 3.8–122.4 |
| 22 islands | separate_prior1 | Straume | 74.3 | 48.9 | 8.1–290.5 |

These are richness predictions for a hypothetical island with the specified effective area, not a direct reconstruction of actual Sunken Islands richness, maximum capacity, or evidence that all this land was suitable habitat. A fragmented plateau may also differ from calibration island units. The upstream map-derived richness distribution was not used in this forward calculation, allowing comparison with it. Different calibration taxonomic coverage and mapping coverage must be reconciled before interpreting any deficit as habitat fraction.

Reproduce: source("area-estimation/scripts/51-richness-geological-areas.R") from the project root. Results also contain 50% predictive intervals and intervals for the conditional mean, which are distinct from the predictive intervals reported above.
