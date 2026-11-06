# Comparing species-area curve shapes

All four models use Gaussian residuals in natural-log richness. Thus curve shape changes while the observation model stays the same. Lines describe conditional median richness on the original scale. Arithmetic-mean predictions use the lognormal variance correction.

Models: power S=c*A^z; exponential S=c*exp(b*A); Michaelis-Menten saturation S=K*A/(H+A); and curved log-log ln(S)=a+b*ln(A)+c*ln(A)^2. Saturation has two mean parameters; the quadratic has three. Residual variance counts as an additional parameter in AICc/BIC.

| Dataset | Model | Log CV RMSE | Count CV RMSE | Delta AICc |
|---|---|---:|---:|---:|
| All | Power | 0.969 | 28.21 | 7.37 |
| All | Exponential | 1.187 | 356.22 | 31.34 |
| All | Saturating | 1.318 | 48.09 | 41.75 |
| All | Curved log-log | 0.904 | 34.68 | 0.00 |
| Filtered | Power | 0.883 | 29.38 | 1.63 |
| Filtered | Exponential | 1.180 | 310.93 | 27.21 |
| Filtered | Saturating | 1.227 | 53.98 | 28.31 |
| Filtered | Curved log-log | 0.858 | 28.86 | 0.00 |

## Interpretation

The exponential curve performs worse than the power relationship on both datasets. The curved log-log model has the best log-scale predictive score: a 6.7% RMSE reduction relative to power with all islands and 2.9% with filtered islands. AICc favours curvature by 7.37 with all islands, but only 1.63 after filtering, so the filtered power and curved models are close. The tested saturating curve is poorly supported. Count-scale prediction does not consistently favour curvature: on the full dataset its count RMSE is worse than the power model.

These are pooled relationships without region or type effects. Apparent curvature may partly reflect differing groups rather than a universal nonlinear mechanism. The previous area-plus-type lognormal model on the full dataset had AICc 162.37, slightly better than this pooled curved model (164.44). No prior Bayesian area estimates have been replaced by these descriptive fits.

Compare models within each dataset: holdouts are identical within a dataset but the full and filtered datasets have different test populations. AICc is comparable across these four shapes within a dataset only. Log-scale RMSE gives small and large islands more balanced influence; count RMSE is sensitive to species-rich islands. Different scores can prefer different curves.

Saturation uses profile least squares over H=1e-6 to 1e12 km², profiling log K analytically. Full and leave-one-out fits record whether the numerical boundary was approached. This tests the Michaelis-Menten curve, not every possible saturating relationship. Quadratic log-log curves can turn down or up outside the data range; they are descriptive and should not be extrapolated blindly for inverse area estimation.

The filter removes 16 Australasian records with richness <6 plus Reunion. Richness-based selection can change the fitted relationship and is not evidence of incomplete inventories. No source records or earlier Bayesian area results were changed.

Run area-estimation/scripts/07-compare-relationships.R from the repository root. Models, full predictions and scores are saved in results/curve-comparison; plots are in plots/curve-comparison.
