# Three-component mixture on all 61 islands

No islands excluded. Three latent power regressions use Gaussian log-richness errors and one common residual SD. Mixing weights are constant, with minprior=0.1, matching the previous two-component specification. Region and island type are not supplied to mixture fitting. Component names are ordered at area 10,000 km² and need not retain that ordering if curves cross.

Full fit: 200 random starts, 2,000 iterations followed by continuation up to 20,000 if necessary, at unchanged 1e-8 tolerance. Best converged three-component fit selected. Cross-validation uses 30 starts per held-out island and the same continuation rule. Any component loss is recorded; failed fits are not silently discarded. Two-component and reference scores were reused only after checking identical source island names, areas and richness. Held-out prediction integrates component weights without using held-out richness to choose a component.

| Model | AICc | BIC | Log CV RMSE | Negative log predictive density |
|---|---:|---:|---:|---:|
| Pooled | 171.81 | 177.72 | 0.969 | 1.416 |
| Regions: shared slope | 171.77 | 179.50 | 0.962 | 1.411 |
| Regions: separate slopes | 172.47 | 181.93 | 0.955 | 1.410 |
| Mixture: 2 components | 178.92 | 190.03 | 0.966 | 1.420 |
| Mixture: 3 components | 173.05 | 188.52 | 1.066 | 1.532 |

## Assignment agreement

| Classification | ARI, two components | ARI, three components |
|---|---:|---:|
| region | -0.055 | -0.010 |
| type | 0.062 | 0.079 |
| group | -0.038 | -0.015 |

28/61 islands have membership confidence below 0.8; 9 below 0.6. Best full likelihood reproduced within 1e-4 in 4/200 starts. 8 starts converged and retained three components. All 61 holdout fits converged; 37 retained three components.

Adjusted Rand index is label-invariant and corrects chance agreement (1: identical partitions, approximately 0: little agreement). Hard assignments are for visualization, not estimation; soft-count tables and individual probabilities are supplied. Associations are descriptive, in-sample, and do not establish a causal regional/geological mechanism. Shared variance/minimum weight restrict the tested mixture family. Model-selection exploration and AIC for weakly identified mixtures merit caution; predictive comparison is also reported. No Bayesian area estimates were changed.

Run area-estimation/scripts/16-three-component-mixture.R from the project root. Outputs are in results/three-component-mixture and plots/three-component-mixture.

## Interpretation

The three-component model improves AICc over two components by 5.87, but its leave-one-out log RMSE is 10.4% worse (1.066 versus 0.966) and predictive log density is worse. The pooled/known-region models remain competitive or better. This is evidence against adopting the three-component model for prediction from these data.

The fitted groups contain 5, 25 and 31 islands. The five-island component contains Aru Islands, Panay, Langkawi, Tioman and Réunion and has a negative area exponent (-0.735); the other exponents are 0.421 and 0.570. This small component accommodates contrasting observations rather than a coherent positive species-area relationship. Component labels are ordered only at 10,000 km2.

Group correspondence is weak for both three regions (ARI -0.010) and island type (ARI 0.079). Caribbean islands split 9/7 between Middle/High; Australasian islands occur in all three groups. Most full-data starts collapsed to fewer components, and only 37/61 cross-validation fits retained three components. Thus predictive scores describe the specified fitting procedure, including component loss, rather than a stable three-component partition. Further constrained-mixture variants would constitute new analyses.
