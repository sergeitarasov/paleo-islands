# Richness-based exclusion sensitivity

Excluded Australasian records with fewer than six species and Reunion: 16 + 1 records; 44 of 61 remain. The source workbook and previous analyses are unchanged.

Adjusted R²: 0.618 to 0.691. Log residual SD: 0.954 to 0.862.
Area exponent: 0.498 (95% CI 0.397–0.598) before; 0.492 (0.391–0.593) after filtering.

Predictive comparison uses the SAME 44 held-out retained islands, with each target removed from both training sets. This avoids comparing errors on different test islands. Selection is still based on observed richness, so the comparison is conditional on passing the filter.

Log-scale leave-one-out RMSE: 0.913 (full training) versus 0.883 (filtered training), a 3.3% reduction.
Count-scale RMSE: 32.87 versus 29.38; MAE: 16.97 versus 16.29. Held-out 95% prediction-interval coverage: 95.5% versus 93.2%.

A better fit after this filter does not demonstrate that the excluded records are erroneous or understudied. Removing low response values selectively from one region changes the sample and can bias its richness–area relationship. Confirm inventory completeness independently before using the filtered fit as the primary calibration. AIC values are not compared across different datasets. These confidence intervals do not account for data-driven selection of the exclusion rule.

This run refits the pooled power regression requested for the single-panel plot. It does not refit the grouped Bayesian area-estimation models; their previous results remain unchanged. The reconstructed Sunken Islands richness is concentrated near the exclusion threshold, making extrapolation of this filtered calibration particularly sensitive to the inventory assumption.

Run: source("area-estimation/scripts/06-filtered-regression.R") from the repository root. Results and figures are in results/filtered-regression and plots/filtered-regression.
