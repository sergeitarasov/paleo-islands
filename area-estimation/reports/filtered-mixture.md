# Latent mixture versus Caribbean/rest regression

Dataset: the latest 39-island filter (16 Caribbean, 23 other). Both mixture components are power regressions with Gaussian log-richness residuals and a common variance, matching the grouped regression variance assumption. No region or island-type labels enter the mixture. Mixing weights are constant across area. A minimum component-weight setting of 0.1 and shared variance avoid tiny-variance solutions. Thus this is a restricted two-component mixture, not a search across every possible mixture.

200 independent starts were tried for the full mixture; 30 per leave-one-out fold. Known-region models use the held-out island's known region; mixture predictions use the fitted mixing weights and never use held-out richness to choose a component. Log predictive density evaluates the full mixture distribution, including its possible multimodality. All models use ML residual variance for cross-validation density and arithmetic-mean correction; count RMSE can differ slightly from previous scripts that used residual-df variance.

| Model | AICc | BIC | Log CV RMSE | Negative log predictive density |
|---|---:|---:|---:|---:|
| Pooled | 79.88 | 84.19 | 0.649 | 1.016 |
| Regions: shared slope | 76.97 | 82.45 | 0.623 | 0.975 |
| Regions: separate slopes | 68.70 | 75.20 | 0.548 | 0.856 |
| Latent mixture | 83.21 | 90.56 | 0.674 | 1.167 |

Best label-matched agreement with Caribbean/rest: 64.1%. Adjusted Rand index: 0.055 (1 means identical partitions; 0 is chance-adjusted agreement). 10 islands have membership confidence below 0.8, including 0 below 0.6. Best full-data likelihood reproduced in 199/200 starts.

Hard membership is used only for summaries and point colours. Parameter estimation uses soft probabilities. Component names refer to richness at 10,000 km²; component lines may cross elsewhere. Membership-region agreement is descriptive, evaluated on the fitting data. It does not establish independently validated regional classes. The response-based prior exclusions can change apparent clusters. No Bayesian area estimates have been updated.

Outputs: model-comparison.csv, group/region-confusion.csv, component-parameters.csv, island-membership.csv, membership-agreement.csv, mixture-starts.csv and complete held-out predictions. Reproduce with area-estimation/scripts/12-filtered-mixture.R from the project root.

## Conclusion

Known regions with separate slopes outperform the latent mixture: delta AICc 14.51; log prediction error 0.548 versus 0.674 (mixture 23.1% worse). Predictive density agrees. Caribbean membership is 9 lower/7 higher; other islands 7 lower/16 higher. Thus the latent components do not closely match geography (adjusted Rand index 0.055). The result is stable across 199/200 starts.
