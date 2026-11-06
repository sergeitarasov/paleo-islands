# Latent mixture versus Caribbean/rest regression

Dataset: all 61 original records (16 Caribbean, 45 other), with no exclusions. Both mixture components are power regressions with Gaussian log-richness residuals and a common variance, matching the grouped regression variance assumption. No region or island-type labels enter the mixture. Mixing weights are constant across area. A minimum component-weight setting of 0.1 and shared variance avoid tiny-variance solutions. Thus this is a restricted two-component mixture, not a search across every possible mixture.

200 independent starts were tried for the full mixture; 30 per leave-one-out fold. Known-region models use the held-out island's known region; mixture predictions use the fitted mixing weights and never use held-out richness to choose a component. Log predictive density evaluates the full mixture distribution, including its possible multimodality. All models use ML residual variance for cross-validation density and arithmetic-mean correction; count RMSE can differ slightly from previous scripts that used residual-df variance.

| Model | AICc | BIC | Log CV RMSE | Negative log predictive density |
|---|---:|---:|---:|---:|
| Pooled | 171.81 | 177.72 | 0.969 | 1.416 |
| Regions: shared slope | 171.77 | 179.50 | 0.962 | 1.411 |
| Regions: separate slopes | 172.47 | 181.93 | 0.955 | 1.410 |
| Latent mixture | 178.92 | 190.03 | 0.966 | 1.420 |

Best label-matched agreement with Caribbean/rest: 67.2%. Adjusted Rand index: -0.038 (1 means identical partitions; 0 is chance-adjusted agreement). 61 islands have membership confidence below 0.8, including 34 below 0.6. Best full-data likelihood reproduced in 1/200 starts.

Hard membership is used only for summaries and point colours. Parameter estimation uses soft probabilities. Component names refer to richness at 10,000 km²; component lines may cross elsewhere. Membership-region agreement is descriptive, evaluated on the fitting data. It does not establish independently validated regional classes. No records were excluded in this run. The fitted dataset is the original workbook. No Bayesian area estimates have been updated.

Outputs: model-comparison.csv, group/region-confusion.csv, component-parameters.csv, island-membership.csv, membership-agreement.csv, mixture-starts.csv and complete held-out predictions. Reproduce with area-estimation/scripts/15-unfiltered-mixture.R from the project root.

## Interpretation

The unfiltered mixture is not supported: AICc is 7.11 worse than pooled, while log-scale prediction error is almost identical (0.966 versus 0.969) and predictive log density is slightly worse. Log-likelihood improves by just 0.010. Every island has membership confidence below 0.8; 34/61 are below 0.6. Hard assignments (6/55) are not meaningful biological groups. Geographic adjusted Rand agreement is -0.038. The clearer filtered split is not robust to restoring excluded observations.

Ten of 200 full-data starts hit the iteration limit; the best converged fit was selected. All final likelihoods were within 0.1 log units. Slow cross-validation fits were continued for up to 20,000 iterations at unchanged 1e-8 tolerance; all final holdout fits converged. Exact component labels remain weakly identified.
