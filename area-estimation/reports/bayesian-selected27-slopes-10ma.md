# Sunken Islands effective area at 10 Ma: pooled versus separate slopes

Both analyses use the same 27-island calibration data and the full saved combined richness distribution from 2,000 stochastic maps at 10 Ma. No truncation to its 95% interval is applied. No mainland-distance predictor is included. Separate slopes also include separate intercepts for Caribbean and non-Caribbean islands; Sunken Islands predictions use the non-Caribbean group. Residual variance is shared across groups.

| Model | Mean km² | Median km² | Equal-tailed 95% credible interval km² |
|---|---:|---:|---:|
| Pooled | 10,752 | 1,322 | 41–80,842 |
| Separate slopes, non-Caribbean prediction | 7,375 | 1,702 | 56–46,819 |

The common area prior is log-uniform from 10 to 1,000,000 km². Regression coefficients and residual variation among islands are integrated into inverse prediction. For each possible richness value, the posterior-predictive likelihood is integrated over regression draws, multiplied by the area prior and normalized; those area distributions are then mixed using the unchanged map richness probabilities. Thus the regression does not reweight the upstream maps.

Regression priors: intercept Normal(ln10, SD 2), baseline slope Normal(0.4, SD 0.5), residual SD half-Normal(SD 1.5). The separate-slope model additionally uses Normal(0, SD 1.5) for the group intercept contrast and Normal(0, SD 0.5) for the slope contrast. All priors are independent on these coefficient coordinates. Predictor is ln(area/10,000); response is ln(richness). These are proper-prior Bayesian fits, unlike the reference-prior predictive calculations used in matched cross-validation.

Four chains per model ran 30,000 iterations, discarding 5,000 warmup, retaining 100,000 draws per model. Maximum Rhat was 1.0011; minimum bulk/tail ESS across models exceeded 5,000. For integration, 16,000 draws distributed across all chains were used. Halving those draws changed baseline medians and upper limits by less than 0.25%. The pooled result reproduces the previous calculation. The plotted distributions were visually checked.

Separate slopes raise the median but shorten the right tail, lowering both the mean and upper limit. The mean exceeds the median because both distributions are strongly right-skewed. A narrower interval does not establish better reliability: the matched held-out comparison favoured pooled slopes for area interval scores and found non-Caribbean misses under separate slopes. Retain pooled as the baseline and separate slopes as sensitivity.

These are richness-equivalent effective areas, not direct geological measurements. They depend on the selected islands, area and regression priors, model form, and upstream mapping assumptions. Model-selection and inventory biases are not propagated. Area-prior sensitivity calculations are saved alongside baseline results; the uniform-area prior produces considerably larger estimates.

## Reproduce

Run from the paleo-islands project root:

```r
source('area-estimation/scripts/41-fit-selected27-slopes.R')
source('area-estimation/scripts/42-estimate-selected27-slopes-area.R')
```

Results: `area-estimation/results/bayesian-selected27-slopes-10ma/`.
Plot: `area-estimation/plots/bayesian-selected27-slopes-10ma/comparison.png` (also PDF).
