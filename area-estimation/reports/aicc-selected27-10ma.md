# AICc-averaged Sunken Islands effective area at 10 Ma

The three Gaussian log-richness versus log-area models use the same 27-island dataset. Weights are exp(-Delta AICc/2), normalized over these three candidates. For models with region effects, predictions use the non-Caribbean group.

| Model | AICc | Weight | Mean area km² | Median km² | 95% interval km² |
|---|---:|---:|---:|---:|---:|
| Pooled | 69.5916 | 0.073406 | 10,752 | 1,322 | 41–80,842 |
| Separate intercepts, common slope | 68.5965 | 0.120728 | 6,459 | 420 | 15–45,551 |
| Separate intercepts and slopes | 64.7997 | 0.805867 | 7,375 | 1,702 | 56–46,819 |
| AICc mixture | — | 1 | **7,512** | **1,484** | **35–48,771** |

The complete combined 2,000-map richness distribution at 10 Ma is propagated without truncating to its 95% interval. Each component includes regression posterior uncertainty and residual variation among islands. All components use the same log-uniform area prior of 10–1,000,000 km². Each area distribution is normalized conditional on a richness value, then mixed with unchanged map probabilities. Finally the three model-specific area distributions are mixed with the fixed AICc weights. Quantiles are calculated from that mixture, not by averaging model-specific interval endpoints.

The new separate-intercept fit uses the same priors and sampler as the previous models: baseline intercept Normal(ln10,2), slope Normal(0.4,0.5), intercept contrast Normal(0,1.5), residual SD half-Normal(1.5). The separate-slope model additionally uses slope contrast Normal(0,0.5). Four chains each retained 25,000 iterations after 5,000 warmup. The separate-intercept fit had maximum Rhat 1.000851 and minimum bulk/tail ESS 6,154; prior models passed their diagnostics. Integrations use 16,000 draws per model on a common log-area grid. All component and mixed distributions were checked to sum to one; the final plot was visually checked.

This is an AICc-weighted mixture of Bayesian conditional area distributions, not formal Bayesian model averaging: AICc weights are not posterior model probabilities. The 95% interval is the equal-tailed uncertainty interval of this hybrid mixture. Weights are treated as fixed; uncertainty in weights and candidate-set selection is not propagated. The calculation concerns effective richness-equivalent area and does not independently measure geological land area.

AICc gives separate slopes most weight, despite poorer area interval performance in matched validation. Averaging therefore does not resolve that predictive limitation. The result remains prior-sensitive: with a uniform-area rather than log-uniform prior over the same bounds, the median is 8,157 km² and the interval 429–388,090 km². All five area-prior sensitivities are saved.

Reproduce from the project root by running scripts 43-fit-selected27-intercepts.R and 44-estimate-selected27-intercepts-area.R, then 45-average-selected27-area.R in area-estimation/scripts/. Scripts 41 and 42 supply the already calculated pooled and separate-slope distributions.

Outputs: area-estimation/results/aicc-selected27-10ma/ and area-estimation/plots/aicc-selected27-10ma/.
