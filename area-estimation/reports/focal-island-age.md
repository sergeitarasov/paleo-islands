# Does assigned island age improve richness–area regression?

Updated workbook analysis, 2026-10-28. Source: data/richness/island-richness-area.xlsx, sheet Focal islands ages. All 21 focal records have positive age, area and richness; no imputation was performed. Sri Lanka and the Malay Archipelago/New Guinea additions are not in this tab. The workbook was not modified; a source snapshot and checksum are saved.

## Main result

Adding age produces a modest increase in in-sample R² but no clear AICc improvement. Raw age is effectively tied with area alone (Delta AICc 0.20). Log age slightly improves leave-one-out prediction (RMSE 0.857 versus 0.878; about 2.3%), but this advantage is not robust to removing the five ages explicitly marked Low-confidence. No model has a statistically clear age coefficient under the stated assumptions. This does not demonstrate that biological time has no effect.

All models use Gaussian errors for ln(richness), with ln(area) as the area predictor. Age is tested either in Ma or as ln(age). Regional models allow separate Caribbean/non-Caribbean intercepts and slopes, with a common additive age effect and shared residual variance. AICc counts regression coefficients plus residual variance. Each LOO fit excludes the test island and predicts its log-richness.

| Model | R² | Adjusted R² | AICc | LOO log-RMSE |
|---|---:|---:|---:|---:|
| Area only | 0.728 | 0.714 | 55.88 | 0.878 |
| Area + log age | 0.753 | 0.725 | 56.99 | 0.857 |
| Area + age | 0.763 | 0.737 | 56.08 | 0.869 |
| Regional slopes | 0.782 | 0.743 | 57.88 | 1.018 |
| Regional slopes + log age | 0.796 | 0.744 | 60.50 | 1.002 |
| Regional slopes + age | 0.784 | 0.730 | 61.66 | 1.025 |

The pooled log-age coefficient is 0.154 (95% coefficient CI -0.089 to 0.397; p=0.199). The raw-age coefficient is 0.0153 per Ma (-0.0044 to 0.0349; p=0.120). With regional slopes, the age-test p-values are 0.311 for log age and 0.684 for raw age. These are exploratory, unadjusted tests conditional on assigned ages.

## Sensitivity

| Subset | n | Area-only RMSE | Area + log-age RMSE | Area + raw-age RMSE |
|---|---:|---:|---:|---:|
| Without five low-confidence ages | 16 | 1.012 | 1.087 | 1.053 |
| Without separation-age proxies | 19 | 0.751 | 0.735 | 0.775 |
| Without low-confidence and regional 35 Ma proxies | 12 | 1.113 | 1.087 | 1.002 |

The five explicitly low-confidence ages are Eleuthera, New Providence, Saint John, Saint Thomas and Tortola. The 19-island subset removes Madagascar and New Zealand, whose ages represent separation rather than emergence. The 12-island subset excludes the five low-confidence entries and every remaining assigned 35 Ma regional proxy. Area-only has the lowest AICc in each sensitivity set among the three pooled alternatives. Compare models within each subset, not AICc between different datasets. Removing Madagascar and New Zealand leaves only three non-Caribbean records, making regional slope LOO prediction particularly unstable.

## Interpretation and limits

Age is not uniformly defined here: construction/rock ages, emergence, landscape-deposit ages, regional land persistence, and continental separation are mixed. Several 35 Ma assignments are shared regional proxies rather than independent island dates. Treating them as fixed, equally reliable measurements ignores age error; the coefficient CIs do not propagate geological uncertainty. Archipelagos also receive single age proxies. These issues can obscure or create associations. No new source validation was performed; the supplied source notes were preserved and used to define sensitivity exclusions.

The results justify keeping age as an exploratory sensitivity covariate rather than replacing the area-only calibration on present evidence. This test assesses contemporary richness prediction; it does not validate an age-dependent inverse estimate for Sunken Islands. A comparable definition of habitat persistence time would be needed for that application.

## Outputs

- scripts/53-test-focal-island-age.R: readable analysis, runnable from the project root.
- results/focal-island-age/: input snapshot, coefficients, AICc/LOO results, nested tests, leave-one-out coefficient diagnostics and fitted models.
- plots/focal-island-age/age-effect.png and .pdf: log-age partial residual plot adjusted for area. Shading is the coefficient interval for the log-age contribution, not a predictive interval.
- plots/focal-island-age/prediction-comparison.png and .pdf: matched LOO errors.

Both PNG figures were visually checked. Previous area estimates were not changed.
