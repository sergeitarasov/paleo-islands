# Summary of the 9 Ma biogeographic analyses

Source: `results.csv`  
Compiled: 2026-09-22

## Analysis design and completeness

The compiled results contain 140 maximum-likelihood fits:

- four biogeographic scenarios;
- seven parameter-grouping models per scenario; and
- five optimization starts per parameter grouping.

This gives 4 scenarios x 7 parameter groups x 5 starts = 140 fits. All 140 fits have convergence code 0, finite AIC and log-likelihood values, and a successful parameter-consistency check. The analysis grid is complete.

The five runs within each scenario and parameter group are repeated optimization starts, not independent biological models. Model comparisons below therefore retain only the converged run with the lowest AIC for each parameterization.

## Best fit for each biogeographic scenario

| Scenario | Interpretation | Best group | Best run | Parameters | log likelihood | AIC | Delta AIC | Akaike weight |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| AfrSunken16 | Africa to Sunken Islands to Mauritius | 7 | 1 | 10 | -1035.007901 | **2090.015802** | 0.000000 | approximately 1.000000 |
| OrientalSunken16 | Oriental region to Sunken Islands to Mauritius | 6 | 4 | 9 | -1074.100669 | 2166.201337 | 76.185535 | 2.86 x 10^-17 |
| Afr16 | Direct Africa to Mauritius | 5 | 3 | 9 | -1076.636366 | 2171.272731 | 81.256929 | 2.27 x 10^-18 |
| Oriental16 | Direct Oriental region to Mauritius | 5 | 2 | 9 | -1076.757613 | 2171.515227 | 81.499425 | 2.01 x 10^-18 |

The Africa-Sunken-Islands scenario is overwhelmingly favored in this four-scenario comparison. Its AIC is 76.19 units lower than that of the second-best scenario. The corresponding evidence ratio is approximately 3.5 x 10^16 to 1 relative to the Oriental-Sunken-Islands scenario.

These weights quantify relative support within this candidate set. They do not establish the historical route independently of the model structure, geographic coding, root assumption, sampling fractions, phylogeny, or the 9 Ma epoch boundary.

## Comparison under matched parameter groupings

Africa-Sunken Islands has the lowest AIC in every matched parameter group. This shows that its overall advantage is not produced solely by selecting a more flexible parameterization.

| Parameter group | Winning scenario | Winning AIC | Runner-up | Delta AIC to runner-up | Winner weight within group |
|---:|---|---:|---|---:|---:|
| 1 | AfrSunken16 | 2208.048371 | Oriental16 | 3.005308 | 0.805570 |
| 2 | AfrSunken16 | 2165.467382 | Afr16 | 7.615392 | 0.952831 |
| 3 | AfrSunken16 | 2163.859328 | Afr16 | 7.847194 | 0.969134 |
| 4 | AfrSunken16 | 2154.420822 | Afr16 | 19.296691 | 0.999895 |
| 5 | AfrSunken16 | 2153.706172 | OrientalSunken16 | 16.208728 | 0.999409 |
| 6 | AfrSunken16 | 2090.461493 | OrientalSunken16 | 75.739844 | approximately 1.000000 |
| 7 | AfrSunken16 | 2090.015802 | OrientalSunken16 | 76.996860 | approximately 1.000000 |

Support is moderate in group 1, strong in groups 2 and 3, and overwhelming in groups 4-7.

## Parameterization within the Africa-Sunken-Islands scenario

| Group | Best run | Parameters | log likelihood | AIC | Delta AIC | Conditional weight |
|---:|---:|---:|---:|---:|---:|---:|
| 7 | 1 | 10 | -1035.007901 | **2090.015802** | 0.000000 | 0.555482 |
| 6 | 5 | 9 | -1036.230747 | **2090.461493** | 0.445691 | 0.444518 |
| 5 | 3 | 9 | -1067.853086 | 2153.706172 | 63.690370 | 8.21 x 10^-15 |
| 4 | 4 | 8 | -1069.210411 | 2154.420822 | 64.405020 | 5.75 x 10^-15 |
| 3 | 3 | 6 | -1075.929664 | 2163.859328 | 73.843526 | 5.13 x 10^-17 |
| 2 | 5 | 4 | -1078.733691 | 2165.467382 | 75.451580 | 2.29 x 10^-17 |
| 1 | 4 | 3 | -1101.024186 | 2208.048371 | 118.032569 | 1.30 x 10^-26 |

Groups 7 and 6 are statistically indistinguishable by the usual Delta AIC <= 2 criterion. Group 7 has a slightly lower AIC, but its conditional weight is only 0.555 versus 0.445 for group 6. Both should be retained when interpreting individual rate parameters. Groups 1-5 have essentially no support relative to groups 6 and 7 within this scenario.

## Optimization stability

Repeated starting values produced highly consistent optima. The largest AIC range among the five starts of any parameterization was 0.08323, for Oriental16 group 7. All other within-parameterization ranges were smaller. This indicates that the completed searches converged to effectively the same likelihood within each parameterization.

## Main conclusion

Under the current 16-state, two-epoch analysis with Mauritius emerging at 9 Ma, the results strongly favor colonization from Africa through the Sunken Islands over direct colonization from Africa, direct colonization from the Oriental region, or colonization from the Oriental region through the Sunken Islands.

