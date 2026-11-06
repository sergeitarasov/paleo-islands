# Exploratory area estimate from higher-richness mixture islands

Selected all 23 pink/red points assigned to the higher-richness component in the previous mixture fit: 7 Caribbean and 16 other islands. Refit a single Bayesian power regression to these observations, ignoring region because the requested subset defines this exploratory calibration. This is not a Bayesian mixture fit: hard assignments are treated as fixed, and membership/selection uncertainty is not propagated.

Used the unchanged complete 2,000-map combined richness distribution at 9 Ma (mean 5.5025, median 4, 95% interval 2–21, full range 2–101). All richness values and their frequencies were retained. The same modular inverse calibration marginalizes regression uncertainty, normalizes the area distribution conditional on each richness value, then combines by map frequencies.

## Baseline comparison

Area prior: uniform on log10 area between 10 and 1,000,000 km².

| Analysis | Mean km² | Median km² | 95% interval km² |
|---|---:|---:|---:|
| Filtered non-Caribbean relationship | 6094 | 1766 | 135–35555 |
| Selected higher-richness islands | 2666 | 697 | 93–15055 |

## Area-prior sensitivity

| Prior | Median km² | 95% interval km² |
|---|---:|---:|
| log_uniform_10_1million | 697 | 93–15055 |
| log_uniform_100_100k | 700 | 124–15025 |
| log_uniform_1_10million | 697 | 93–15055 |
| lognormal_10k | 835 | 127–14855 |
| uniform_area_10_1million | 1028 | 141–22908 |

Four chains for each of two regression-prior widths, 5,000 warmup + 25,000 retained iterations each. Maximum Rhat 1.00063; minimum bulk/tail ESS 8021 / 10261. All passed Rhat <1.01 and ESS >400. Integration uses 4,000 draws and a 0.01 log10-area grid; half-draw changes in median/upper limit are at most 0.53% / 0.32%. All area distributions normalize to one.

Model: ln(N)=alpha+beta*ln(A/10,000)+epsilon. Priors alpha ~ Normal(ln10,2), beta ~ Normal(0.4,0.5), sigma ~ half-Normal(0,1.5), with SDs after commas. The wider-prior fit doubles these SDs. Log-sigma sampling includes its Jacobian. MCMC proposal covariance comes from posterior-mode curvature; proposal scale is tuned only during warmup.

Any increased precision is conditional on selecting the higher-richness pattern using the same response data. Hard selection removes low-richness residuals and ignores uncertain memberships, so it can understate population scatter and uncertainty. This experiment does not validate the selected subset as the correct analogue for Sunken Islands. The mixture itself predicted less well than the known-region model. Results are exploratory effective richness-equivalent area, not a primary geological estimate. No original data, prior analyses or stochastic maps were changed.

Scripts: 13-fit-mixture-high-bayesian.R and 14-estimate-mixture-high-area.R. Results/plots/logs are in bayesian-mixture-high subdirectories. calibration-data.csv lists all included islands. Baseline and wider-regression-prior results remain separate.
