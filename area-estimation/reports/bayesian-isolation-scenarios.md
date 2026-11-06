# Sunken Islands area under mainland-distance scenarios

All 61 unfiltered modern islands calibrate ln(N) = alpha + beta ln(A/10,000) + Caribbean/rest intercept and slope differences + gamma [ln(1+distance_km)-ln(1001)]. Sunken Islands use the non-Caribbean group. These are illustrative coast-to-Africa distance scenarios at 9 Ma, NOT a reconstructed palaeodistance distribution. No single scenario is selected or averaged.

The full saved 2,000-map combined richness distribution at 9 Ma is retained. Regression uncertainty is marginalized first within each richness value; each conditional area distribution is normalized before mixing with the unchanged map weights. Residual variation among islands is included. The baseline area prior is log-uniform on 10–1,000,000 km2. The lognormal observation model and continuous log-richness inverse likelihood follow previous Bayesian area analyses.

| Distance to Africa (km) | Median km2 | Mean km2 | 95% credible interval km2 |
|---:|---:|---:|---:|
| 500 | 866 | 11207 | 20–90065 |
| 1000 | 1622 | 18070 | 30–160373 |
| 1500 | 2354 | 23733 | 39–219361 |
| 2000 | 3070 | 28663 | 48–269316 |
| 2500 | 3772 | 33073 | 57–312067 |

Maximum Rhat 1.00167; minimum bulk/tail ESS 4501 / 7002. Four chains per prior width; 5,000 warmup and 25,000 retained draws each. Integration uses 8,000 draws; the half-draw check changes median/upper limits by at most 1.90% / 0.83%.

Priors: reference-group intercept Normal(ln10, SD2); reference area slope Normal(0.4, SD0.5); group intercept difference Normal(0, SD1.5); slope difference Normal(0, SD0.5); distance coefficient Normal(0, SD1.5); residual SD half-Normal(SD1.5). Wider-prior fit doubles all SDs. Rate-free, fixed-covariance random-walk Metropolis uses adaptation only in warmup; the log-sigma Jacobian is included. Area-prior and wider-regression-prior sensitivity are saved in area-summary.csv.

These are effective richness-equivalent areas, conditional on geographic distance, calibration model, priors and upstream maps. They are not physical land-area measurements. Geographic distance could itself depend on island geometry/area; scenarios treat distance as externally fixed. Model-selection uncertainty and palaeodistance uncertainty are not propagated. The current model wins forward richness prediction, but the simpler area-only model had a smaller median inverse-recovery error in exploratory validation. Existing estimates were not overwritten.
