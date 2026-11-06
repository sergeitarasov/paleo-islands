# Excluding Malay Archipelago islands: effective area at 10 Ma

The 23-island analysis excludes Borneo, Java, Sumatra and Sulawesi from the 27-island calibration. The 22-island sensitivity also excludes New Guinea. Sri Lanka remains in both; these are not the original 21-island subset. The non-Caribbean sample sizes are 11, 7 and 6, respectively.

Both pooled and separate-intercept/separate-slope models retain the same regression priors and full 2,000-map richness distribution at 10 Ma. Regional predictions use the non-Caribbean group. The baseline area prior remains log-uniform over 10–1,000,000 km². Full posterior-predictive likelihoods are integrated over regression draws, normalized conditional on each richness value and then mixed with unchanged map frequencies. Intervals are equal-tailed Bayesian credible intervals.

| Islands | Model | Mean km² | Median km² | 95% interval km² |
|---|---|---:|---:|---:|
| 27 | simple_prior1 | 10,752 | 1,322 | 41–80,842 |
| 27 | separate_prior1 | 7,375 | 1,702 | 56–46,819 |
| 23 | simple_prior1 | 15,936 | 1,461 | 33–137,545 |
| 23 | separate_prior1 | 9,832 | 1,774 | 39–67,104 |
| 22 | simple_prior1 | 18,968 | 1,524 | 29–175,314 |
| 22 | separate_prior1 | 10,838 | 1,738 | 33–76,139 |

Exclusion changes medians modestly but increases the mean and upper tail. This is not evidence that the reduced datasets are better: they contain fewer large-area calibration observations. Separate slopes have only seven or six non-Caribbean records after exclusion. All four new fits passed Rhat <1.01 and bulk/tail ESS >400; maximum Rhat was 1.00113 and minimum ESS was 4,912.

Reproduce with area-estimation/scripts/47-fit-no-malay-23.R and 48-estimate-no-malay-23.R, or the corresponding -22.R scripts. Results and plots are in bayesian-no-malay-23-10ma and bayesian-no-malay-22-10ma under area-estimation/results and area-estimation/plots. Prior sensitivity results are also retained.
