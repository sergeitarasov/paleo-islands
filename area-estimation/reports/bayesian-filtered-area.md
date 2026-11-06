# Filtered non-Caribbean area reconstruction at 9 Ma

Used 39 islands (16 Caribbean, 23 other), retaining the agreed exclusions. The Bayesian lognormal regression has separate intercepts and slopes for Caribbean and non-Caribbean groups and a shared residual standard deviation. Inverse prediction uses only the non-Caribbean relationship. It does not include island-type or separate Indian Ocean effects.

The full 2,000-map combined richness distribution at exactly 9 Ma was retained: median 4, mean 5.5025, 95% interval 2–21, range 2–101. No map was discarded or reweighted. Regression coefficients and residual variation are integrated out before normalizing each richness-specific area distribution and combining by its map frequency. This is the same modular inverse-prediction procedure as before.

## Baseline comparison

Area prior: uniform on log10 area, 10–1,000,000 km². Previous estimates below use their oceanic-island setting; the new group model has no geological-type setting. Both data selection and model structure change, so their effects cannot be attributed separately.

| Analysis | Mean km² | Median km² | 95% interval km² |
|---|---:|---:|---:|
| Previous simple: 61 islands | 22379 | 1660 | 26–215860 |
| Previous flexible: 61 islands | 17687 | 815 | 15–173461 |
| Filtered grouped: 39 islands | 6094 | 1766 | 135–35555 |

## Area-prior sensitivity: new model

| Prior | Median km² | 95% interval km² |
|---|---:|---:|
| log_uniform_10_1million | 1766 | 135–35555 |
| log_uniform_100_100k | 1782 | 176–32951 |
| log_uniform_1_10million | 1766 | 134–35569 |
| lognormal_10k | 2262 | 255–32023 |
| uniform_area_10_1million | 4048 | 391–88791 |

## Checks and interpretation

Four independent chains for each of two regression-prior widths; 5,000 warmup and 25,000 retained iterations per chain. Max Rhat 1.00092; minimum bulk/tail ESS 5277 / 7267. All passed Rhat <1.01 and ESS >400. Area integration uses 4,000 posterior draws and a 0.01 log10-area grid; half-draw checks have maximum median/upper-limit change 1.01% / 2.53%. All distributions normalize to one.

Centered predictor x=ln(area/10,000). Priors: Caribbean intercept Normal(ln10,2), Caribbean slope Normal(0.4,0.5), non-Caribbean intercept difference Normal(0,1.5), slope difference Normal(0,0.5), residual sigma half-Normal(0,1.5); values after commas are SDs. The wider alternative doubles every prior SD. Sampling log sigma includes the Jacobian. Proposal covariance is fixed from the posterior-mode Hessian and proposal scale is tuned only during warmup. Both chains and the wider-prior area estimates are saved.

Important extrapolation: retained non-Caribbean islands have at least 6 species and area at least 1811 km². 69.9% of mapped richness values are below that minimum richness; 50.2% of the baseline inferred area distribution lies below the minimum non-Caribbean calibration area. The Caribbean data help estimate shared residual scatter but do not replace missing low-richness non-Caribbean calibration islands.

Successive data-dependent exclusions reduce observed scatter and may narrow intervals without establishing better accuracy. Interpretation remains conditional on the excluded inventories being unsuitable calibration records. Coefficient intervals and area intervals do not include uncertainty about choosing exclusions or grouping. This is effective richness-equivalent area, not necessarily total exposed land area. Previous Bayesian results and original data remain unchanged.

Reproduce from the project root: source('area-estimation/scripts/10-fit-filtered-bayesian.R'), then source('area-estimation/scripts/11-estimate-filtered-area.R'). Results, plots and logs are in the corresponding bayesian-filtered subdirectories.
