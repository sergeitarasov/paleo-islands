# Matched validation of the 21- and 27-island regressions

Date: 2026-10-27

## Main finding

Adding six islands does not materially improve overall richness prediction on the same held-out islands. The 27-island pooled model provides the best area interval score of the four candidates, but all models give broad area distributions. Separate regional slopes improve some point estimates without improving reliable area uncertainty: they miss Réunion and New Zealand, two of the five non-Caribbean test islands.

## Design

Each of the original 21 islands (16 Caribbean and five non-Caribbean) was held out in turn. The two training datasets then contained 20 or 26 islands. Borneo, Java, Sumatra, Sulawesi, New Guinea and Sri Lanka were always additional training observations for the expanded dataset. No held-out richness or area was used to fit its regression.

Four candidates used natural-log richness and natural-log area: pooled intercept/slope, or separate Caribbean/non-Caribbean intercepts and slopes, for each dataset. There was no mainland-distance covariate. Forward validation predicts richness from known area. Inverse validation hides area and infers it from richness and region.

Predictive distributions are exact Student-t distributions for Gaussian linear regression under reference prior p(beta,sigma²) proportional to 1/sigma². They include regression-coefficient and residual-variance uncertainty. Inverse area distributions normalize that predictive likelihood against a common log-uniform area prior of 10–1,000,000 km². Numerical integration uses log10-area spacing 0.0025. These are validation calculations, distinct from the proper-prior MCMC used in the earlier Sunken Islands estimates; those estimates were not changed.

## Same 21 test islands

| Training dataset | Slopes | Richness log-RMSE | Median area error factor | Area 95% coverage | Mean area interval score |
|---|---|---:|---:|---:|---:|
| 21 | Pooled | 0.878 | 2.74 | 21/21 | 2.815 |
| 21 | Separate | 1.018 | 2.50 | 19/21 | 4.085 |
| 27 | Pooled | 0.877 | 2.59 | 20/21 | 2.529 |
| 27 | Separate | 0.874 | 2.26 | 19/21 | 4.425 |

Lower errors and interval scores are better. Area error factor is 10 raised to the median absolute log10 error of inferred median area; 2.74 means a typical multiplicative error of about 2.74-fold. The 95% interval score is calculated on log10 area, penalizing both width and misses: (U-L) + 40(L-truth) if truth<L + 40(truth-U) if truth>U. Coverage alone rewards very broad intervals and should not determine the model choice. Median upper/lower interval ratios were 553, 1618, 297 and 664, respectively. All models therefore retain substantial inverse uncertainty.

Mean absolute log10-area errors were 0.495, 0.523, 0.515 and 0.508, respectively. Thus the smaller median error for separate slopes does not imply uniformly better point recovery. Mean negative log density of the true log10 area was 0.973, 1.224, 0.962 and 1.226; this also favours pooled slopes.

## Non-Caribbean held-out islands

These are Réunion, Mauritius, Madagascar, New Caledonia and New Zealand (18 species).

| Training dataset | Slopes | Richness log-RMSE | Median area error factor | Area 95% coverage |
|---|---|---:|---:|---:|
| 21 | Pooled | 1.466 | 6.95 | 5/5 |
| 21 | Separate | 1.870 | 9.33 | 3/5 |
| 27 | Pooled | 1.241 | 6.98 | 4/5 |
| 27 | Separate | 1.533 | 5.74 | 3/5 |

Adding the six islands improves richness prediction for this group, but not enough to establish dependable inverse area recovery. In the 27-island separate-slope model, Réunion is inferred at median 59 km² (95% interval 11–899), versus known area 2,508 km²; New Zealand at 17,993 km² (2,399–107,470), versus 263,310 km². These misses illustrate the limits of transferring richness to area, rather than a reason to remove these islands after inspecting the results. There are only five non-Caribbean test cases.

## Separate slopes

On all 27 observations, the Caribbean slope is 0.358 (95% coefficient CI 0.200–0.515), and the non-Caribbean slope is 0.682 (0.474–0.889). The interaction test has p=0.017. With 21 observations the slopes are 0.358 and 0.648, interaction p=0.103. These intervals assume a shared residual variance across groups.

Evidence of differing slopes in the fitted data does not establish improved prediction. On the matched 21 test cases, the 27-island separate-slope model barely changes richness RMSE, and worsens area interval and density scores relative to pooled slopes. It improves Caribbean richness predictions but performs worse on the five other islands.

## Paired comparisons and limitations

Five thousand paired bootstrap resamples retained the 16/5 group composition. For 27 pooled versus 21 pooled, the mean squared-log-richness error change was -0.0015 (descriptive 95% bootstrap interval -0.3745 to 0.3274), while the area interval-score change was -0.286 (-0.474 to -0.032). For 27 separate versus 27 pooled, richness squared-error change was -0.0047 (-0.2944 to 0.2684), while area interval-score change was +1.896 (+0.422 to +3.410). Adding six islands did help the separate-slope model compared with its 21-island version: squared-richness-error change -0.273 (-0.497 to -0.072).

These bootstrap intervals describe variability across these test islands, not a new independent experiment. Cross-validation fits share training observations; islands can share geography, history and inventory biases. Dataset selection was exploratory. The common test set does not validate prediction on the six added islands or on all islands worldwide. The finite area prior influences inverse estimates and constrains large-island upper intervals, especially Madagascar. Coverage fractions based on 21, and particularly five, cases are imprecise.

For now, the 27-island pooled model is a reasonable simple baseline for area estimation, with separate slopes retained as a sensitivity analysis. These results do not support presenting a precise Sunken Islands area estimate, or choosing the interaction model solely because it improves in-sample fit.

## Reproduce and outputs

From the paleo-islands project root:

```r
source('area-estimation/scripts/40-matched-21-27-validation.R')
```

- Results: `area-estimation/results/matched-21-27/`
- Plots: `area-estimation/plots/matched-21-27/`
- `held-out-predictions.csv`: all 84 individual predictions and intervals.
- `summary.csv`: overall and regional metrics.
- `direct-paired-contrasts.csv`: direct dataset and slope comparisons.
- `full-data-slopes.csv`: descriptive slopes and interaction tests.
- `session-info.txt`: R environment.

Both PNG figures were visually checked. PDF export only reported a harmless dash-substitution warning in the area-prior caption.
