# Age and isolation: unfiltered richness analysis

Isolation was refitted on all 61 islands, without the previous richness/geographic exclusions (script 18). Mainland distance improves the area-only model: AICc 171.81 to 164.35; leave-one-out log RMSE 0.969 to 0.915. For Caribbean/rest-specific slopes, mainland distance improves AICc 172.47 to 160.45 and log RMSE 0.955 to 0.865 (9.4% reduction). This is the best of the 12 isolation candidates by AICc and log RMSE. Distance to a landmass of at least 10,000 square km alone does not improve these baselines. Isolation does not improve held-out performance after accounting for island type. See isolation-regression.md for all candidates and uncertainty.

Age cannot be evaluated across all 61 islands: only seven have point-age proxies. The following comparisons use identical complete-case islands within each sample. The three one-sided bounds were not treated as exact ages, and missing ages were not imputed as zero. This is exploratory: volcanic construction, rock ages and reconstructed land appearance are heterogeneous covariates.

| Sample | Model | n | AICc | Adjusted R2 | LOO log RMSE |
|---|---|---:|---:|---:|---:|
| All point-age proxies | Area | 7 | 29.99 | 0.527 | 0.918 |
| All point-age proxies | Area + age | 7 | 43.24 | 0.468 | 1.811 |
| All point-age proxies | Area + mainland | 7 | 42.93 | 0.491 | 1.132 |
| All point-age proxies | Area + mainland + age | 7 | 84.34 | 0.377 | 3.031 |
| Without conflicting Sumbawa | Area | 6 | 32.63 | 0.411 | 1.049 |
| Without conflicting Sumbawa | Area + age | 6 | 61.99 | 0.294 | 2.134 |
| Without conflicting Sumbawa | Area + mainland | 6 | 61.66 | 0.331 | 1.374 |
| Without conflicting Sumbawa | Area + mainland + age | 6 | NA | 0.082 | 3.773 |

Adding age worsens AICc and leave-one-out prediction on the seven-island subset. The age coefficient in area + age is 0.496 (ordinary model-based 95% CI -1.561 to 2.554; p = 0.540). Excluding conflicting Sumbawa gives the same qualitative result. AICc is undefined for the most complex six-island model because n-k-1 is zero. Small samples and high leverage make these predictive comparisons unstable. Do not interpret failure to improve fit as proof that geological age is biologically unimportant.

For now, mainland distance is the useful new candidate predictor; age is not supported by the available data. The nearest-landmass name and database ID identify the distance target and were not fitted as high-cardinality categorical predictors. The 10,000 km2 threshold is constant and cannot explain variation. Contemporary distances do not establish palaeodistances for Sunken Islands, and no Bayesian inverse area analysis was changed. No AIC comparison is made between the seven-island age models and the 61-island isolation models.
