# Island richness–area regression: pilot

## Recommendation

Start with a power species–area relationship that includes island type. Retain a negative-binomial model with type-specific slopes and region as a more flexible alternative. There is no single clear winner across prediction criteria, and the present data do not support adding a latent mixture simply to improve fit.

The interpretable baseline is `log(n_sp) ~ log(area_km2) + type`. Its estimated area exponent is **0.458 (95% confidence interval 0.355–0.560)**, adjusted R² = **0.686**. A tenfold area increase corresponds to approximately **2.87 times** the median richness, conditional on island type. This is an association, not a causal area effect.

## Data and procedure

The input is `data/richness/island-richness-area.xlsx`, sheet 1: 61 records, unchanged. There are 41 Australasian, 16 Caribbean and four Indian Ocean records; types comprise 21 continental fragments, nine continental-shelf islands and 31 oceanic islands. The Caribbean has no continental-shelf records. A frozen input snapshot is included as `input-data.csv`.

We fitted 19 models: six predictor structures each under lognormal, negative-binomial and Poisson responses, plus one two-component lognormal regression mixture. Structures were area alone; area + region; area + type; area + region + type; region-specific slopes + type; and type-specific slopes + region. Region and type are fixed effects, not random intercepts: each has only three levels, and the cross-classification is sparse. This is distinct from the finite mixture, whose two unobserved components have separate intercepts and slopes and a shared residual variance.

The mixture used flexmix, 50 starts for the full fit and 15 starts per holdout, a common variance and minimum component weight setting of 0.1. Those restrictions reduce unstable tiny-component solutions. This is one constrained two-component mixture test, not a rejection of every possible mixture formulation. Method: https://search.r-project.org/CRAN/refmans/flexmix/html/flexmix.html

Every model was evaluated by leaving out each island once (61 holdouts). No held-out richness was used to assign mixture membership. Lognormal predictions were corrected to arithmetic means using exp(prediction + residual variance/2). Mixture mean predictions average component means using fitted mixing weights. The principal shared score is RMSE between observed log richness and the log of predicted mean richness; count RMSE and count MAE are also reported. These criteria answer different prediction questions. The additional `log_RMSE` column uses predicted expected log richness for lognormal/mixture models but log mean for count models, so it should not be used as a common cross-family ranking.

AICc and BIC are interpreted within the log-response models or within count models, never by comparing raw AIC values across transformed continuous and discrete responses. Model selection and evaluation use this same small dataset; the winning score is exploratory, not an independent final validation.

## Main results

| Candidate | Leave-one-out log-mean RMSE | Count RMSE | AICc |
|---|---:|---:|---:|
| Negative binomial: type-specific slopes + region | 0.929 | 31.72 | 397.87 |
| Negative binomial: area + type | 0.954 | 44.35 | 402.14 |
| Lognormal: area + type | 0.981 | 39.67 | 162.37 |
| Lognormal: area + region + type | 1.022 | 39.83 | 166.11 |
| Two-component lognormal mixture | 1.064 | 28.38 | 178.93 |
| Lognormal: area only | 1.068 | 28.21 | 171.81 |

**Count-model result.** Negative binomial with type-specific slopes + region has the lowest log-mean prediction error and count-model AICc. However, the shared-slope, type-only negative binomial has lower BIC (411.60 versus 413.34). The paired squared-log-error improvement is 0.047 with descriptive SE 0.077; this is not strong evidence for a definitive prediction advantage. The SE is descriptive because leave-one-out training sets overlap. The best log score still corresponds to an error factor of roughly exp(0.929) = 2.53, not tight prediction.

**Lognormal result.** Area + type is the best AICc/BIC lognormal candidate. Adding region increases AICc by 3.74. Type-specific slopes increase AICc by 5.41. Continental-shelf richness is higher conditional on area in this model, but there are only nine shelf records and they mostly belong to one region. Do not interpret that coefficient as an isolated geological effect.

**Mixture result.** Its AICc is 16.56 higher than the type-only lognormal model and 7.12 higher than area only. Prediction is essentially unchanged from area only. Hard component assignment is 52/9, despite mixing proportions 0.548/0.452, indicating considerable overlap. There is no useful evidence here for replacing observed island groups with two latent richness regimes.

**Poisson result.** A Poisson area + region model has the lowest raw count RMSE (17.67) and MAE (9.63). Thus there is no universal predictive winner. However, even the more flexible Poisson type-slope model has Pearson dispersion 5.38, contradicting its equidispersion assumption. Its likelihood is markedly poorer than the corresponding negative binomial. Ordinary Poisson uncertainty would be too narrow; its raw-count prediction success is not sufficient reason to use its uncertainty model.

## Sensitivity and limits

- Aru Islands, Madagascar and Réunion have the largest Cook's distances in the common-slope region + type model. They remain in the main analysis. Aru has one listed species over about 8,117 km²; its richness and geographic scope merit biological review rather than automatic removal.
- Excluding four archipelago records gives a common area exponent of 0.421 (0.330–0.513), versus 0.442 (0.334–0.550) in the full region + type model. Excluding Mauritius and Réunion gives 0.426 (0.321–0.531). Positive area association is stable in these checks.
- Model preference is less stable: after excluding archipelago records, region-specific slopes and the simpler type-only lognormal model are essentially tied in AICc. See `sensitivity-model-selection.csv`. These are alternative datasets, so compare AICc only within each subset.
- Only four Indian Ocean observations are available. This does not establish a well-calibrated Indian Ocean relationship. Regional source differences and tentative geological classifications are not resolved by adding regression terms. The user-confirmed richness scope was retained.
- All records have positive richness. The pilot does not model island occupancy or selection of islands with known beetles. If islands with zero richness were systematically omitted, that sampling process must be addressed before treating this as an unconditional island model.
- Observations are treated as independent; geographic/archipelago dependence, incomplete inventories and area error are not modeled. Historical/subfossil species in Mauritius and the Réunion record were retained. Classification uncertainty was not propagated.
- Do not invert a fitted line and its confidence band to obtain a Sunken Islands area interval. That needs inverse calibration accounting for richness uncertainty, residual predictive variation, regression uncertainty and the assumed island type/region. This pilot has not estimated Sunken Islands area.

## Figures

- `negative_binomial.pdf/png`: best log-mean-score count model, by type and region. Ribbons are approximate pointwise 95% confidence intervals for expected richness, **not prediction intervals for a new island**. Type/region curves appear only inside their observed area ranges and where at least two observations exist.
- `species_area.pdf/png`: common-slope lognormal model with region and type. Curves and bands describe conditional median richness (exponentiated log-scale fitted means), not arithmetic means.
- `model_comparison.pdf/png`: common leave-one-out log-mean prediction score for all 19 models.
- `predictive_check.pdf/png`: held-out arithmetic-mean predictions for the simple type model, flexible negative binomial and mixture. Grey line denotes perfect prediction.

## Reproduce

Open the paleo-islands RStudio project and run:

```r
# Install once if needed:
# install.packages(c("readxl", "MASS", "ggplot2", "flexmix"))
source("area-estimation/scripts/01-pilot-regression.R")
```

Or, from the project root: `Rscript area-estimation/scripts/01-pilot-regression.R`.

The script reads the current workbook, writes tables/models in `area-estimation/results/pilot-regression`, figures in `area-estimation/plots/pilot-regression`, and session information in `area-estimation/logs/pilot-regression`, and records `session-info.txt`. It leaves the workbook unchanged. `run.log` records this pilot execution. Numerical tables and the source snapshot accompany the figures; models are saved in `models.rds`.


Reorganized into `area-estimation` without rerunning or changing the fitted analyses. Filenames above refer to their respective results, plots or logs directories. The historical run.log records the original execution.
