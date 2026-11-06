# Caribbean versus other islands

Three lognormal power models were compared: one pooled relationship; two group intercepts with a common slope; separate group intercepts and slopes. Grouping is Caribbean versus Australasian + Indian Ocean. Leave-one-island-out tests use the same observations for all models within each dataset.

| Data | Model | Adjusted R² | Log CV RMSE | Count CV RMSE | Delta AICc |
|---|---|---:|---:|---:|---:|
| All | Pooled | 0.618 | 0.969 | 28.21 | 0.04 |
| All | Two intercepts | 0.626 | 0.962 | 28.35 | 0.00 |
| All | Two slopes | 0.630 | 0.955 | 24.33 | 0.70 |
| Filtered | Pooled | 0.691 | 0.883 | 29.38 | 15.11 |
| Filtered | Two intercepts | 0.787 | 0.742 | 33.95 | 0.00 |
| Filtered | Two slopes | 0.783 | 0.761 | 33.41 | 2.42 |

## Interpretation

For the filtered data, a common slope with two intercepts is preferred: log-scale prediction error improves by 16.0%, adjusted R² increases from 0.691 to 0.787, and AICc improves by 15.11. Separate slopes are not needed (AICc 2.42 worse than shared slope). However, count-scale RMSE increases from 29.38 to 33.95, so improvement is specific to proportional/log-scale prediction rather than every criterion.

For all islands, group models improve log-scale RMSE by only 0.7–1.4%, and AICc differences are below one. Thus the split is useful for the filtered sample but is not strongly supported as a universal difference in the complete data.

Full data: 16 Caribbean and 45 other islands. Filtered data: the same 16 Caribbean and 28 other islands. The filter removes 16 Australasian records below six species and Reunion, leaving 44 records. Compare AICc only within a dataset. All models use Gaussian errors on log richness. Plot bands describe uncertainty in fitted median richness, not new-island prediction intervals.

The filtered comparison deliberately retains low-richness Caribbean records while excluding low-richness Australasian records. This response-based selection can create or exaggerate a group contrast. Differences may also reflect geology, source/inventory differences and isolation. The split is not proof of a causal regional effect; assess the all-data result alongside the filtered one. These exploratory comparisons reuse the dataset to select hypotheses and are not independent validation.

Previous Bayesian area estimates are unchanged. Scripts use the authoritative workbook without editing it. Run area-estimation/scripts/08-caribbean-versus-rest.R from the project root.
