# Additional island exclusions

Starting with the previous exclusions (Australasian richness below six species and Reunion), additionally excluded Tioman, Langkawi, Sibuyan, Lombok and Bali. This leaves 39 records: 16 Caribbean and 23 other islands. The source workbook and earlier results are unchanged.

Compared pooled power regression, common slope with group-specific intercepts, and separate slopes. All use Gaussian residuals on natural-log richness. Leave-one-out predictions use identical 39-island test sets for the three candidates.

| Model | Adjusted R² | Log CV RMSE | Count CV RMSE | Delta AICc |
|---|---:|---:|---:|---:|
| Pooled | 0.847 | 0.649 | 30.66 | 11.18 |
| Two intercepts | 0.863 | 0.623 | 32.98 | 8.28 |
| Two slopes | 0.893 | 0.548 | 23.22 | 0.00 |

Compare models on these same 39 records; do not compare AICc with fits on the earlier 44 or 61 records. Successive exclusions based on inspecting richness and the fitted plot can make apparent fit and predictive scores optimistic for the original population. An independent inventory-quality justification is needed before adopting this filtered calibration as primary. Low-richness Caribbean records remain included. This analysis does not revise Bayesian area estimates.

Run area-estimation/scripts/09-additional-exclusions.R from the project root. Outputs: results/additional-exclusions and plots/additional-exclusions.

Separate slopes are preferred on these 39 selected islands: log-scale prediction error is 15.7% lower than pooled and AICc improves by 11.18. Estimated area exponents are 0.358 for Caribbean and 0.646 for other islands.
