# Isolation and species-area regression

All 61 original islands are retained. The response is log richness and errors are Gaussian (a lognormal richness model). Area is log transformed; isolation uses log(1 + distance in km). Twelve exploratory candidates compare area-only, Caribbean/rest-specific slopes, and island-type intercepts, with no isolation, either distance, or both. Geographic groups/types are observed covariates, not inferred mixture components. All models use the same islands. No Bayesian area estimates were updated.

| Model | AICc | Adjusted R² | LOO log RMSE | Negative log predictive density |
|---|---:|---:|---:|---:|
| Caribbean/rest: Mainland | 160.45 | 0.703 | 0.865 | 1.309 |
| Island type: None | 162.37 | 0.686 | 0.917 | 1.413 |
| Caribbean/rest: Both | 162.40 | 0.700 | 0.881 | 1.332 |
| Island type: Mainland | 162.68 | 0.692 | 0.926 | 1.422 |
| Area only: Both | 164.27 | 0.676 | 0.909 | 1.352 |
| Area only: Mainland | 164.35 | 0.669 | 0.915 | 1.355 |
| Island type: Both | 164.55 | 0.690 | 0.939 | 1.442 |
| Island type: Large landmass | 164.63 | 0.682 | 0.930 | 1.434 |
| Area only: None | 171.81 | 0.618 | 0.969 | 1.416 |
| Caribbean/rest: None | 172.47 | 0.630 | 0.955 | 1.410 |
| Area only: Large landmass | 172.48 | 0.622 | 0.971 | 1.420 |
| Caribbean/rest: Large landmass | 174.04 | 0.629 | 0.965 | 1.425 |

Lowest LOO RMSE: Caribbean/rest: Mainland
Correlation between log mainland and log large-landmass distances: -0.050.

Coefficient intervals are ordinary model-based 95% confidence intervals, conditional on the selected model; they do not account for searching across candidates. Paired LOO log-density differences and approximate standard errors are saved separately. These are exploratory predictive comparisons, not evidence that isolation causes richness differences. Records share regions and sometimes archipelagos, so island-level cross-validation does not test transfer to a wholly new region. Contemporary isolation is not a reconstruction of isolation at 9 Ma. Applying an isolation model to Sunken Islands would require their palaeogeographic distance and its uncertainty.

## Interpretation

Adding mainland distance reduces LOO log RMSE by 5.6% for area alone and 9.4% for Caribbean/rest-specific slopes.
The Caribbean/rest plus mainland-distance model has the lowest AICc and LOO error among these candidates. Richness decreases with mainland distance after adjusting for area. Distance to any large landmass does not improve held-out predictions. Adding isolation to the island-type model does not improve its LOO performance; this is not a universal gain across model specifications.
For the Caribbean/rest model, the total held-out log-density improvement from mainland distance is about 6.2 with approximate SE 5.7. Thus predictive evidence is promising but uncertain. Ordinary coefficient confidence intervals and AICc provide complementary, conditional evidence; neither proves transferability to an extinct archipelago.
