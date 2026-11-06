# Oceanic-only richness regression

All 31 oceanic islands retained; no richness-based exclusions. Region counts:

Australasian    Caribbean Indian Ocean 
          18           11            2 

Response: log richness; Gaussian errors. Area uses natural log; mainland distance uses log(1 + km). Group is Caribbean versus other islands. Age is omitted because of sparse, heterogeneous proxies. Geographic classifications were not changed.

| Model | AICc | Adjusted R2 | LOO log RMSE |
|---|---:|---:|---:|
| Area | 54.76 | 0.763 | 0.551 |
| Area + isolation | 56.29 | 0.763 | 0.563 |
| Area + group + isolation | 57.31 | 0.769 | 0.557 |
| Area + group | 57.37 | 0.755 | 0.565 |
| Area x group + isolation | 58.45 | 0.774 | 0.545 |
| Area x group | 59.66 | 0.750 | 0.563 |

## Matched predictive comparison

Both approaches predict the same 31 held-out oceanic islands. Oceanic-only training contains the other 30 oceanic islands; full-data training contains the other 60 islands. No held-out observation is used in fitting. AICc is compared only within the oceanic subset, never across datasets.

| Model | Training data | Oceanic LOO log RMSE | Negative log predictive density |
|---|---|---:|---:|
| Area | Oceanic only | 0.551 | 0.909 |
| Area | All island types | 0.539 | 1.027 |
| Area + isolation | Oceanic only | 0.563 | 0.934 |
| Area + isolation | All island types | 0.557 | 0.988 |
| Area + group | Oceanic only | 0.565 | 0.946 |
| Area + group | All island types | 0.578 | 1.041 |
| Area + group + isolation | Oceanic only | 0.557 | 0.932 |
| Area + group + isolation | All island types | 0.543 | 0.946 |
| Area x group | Oceanic only | 0.563 | 0.947 |
| Area x group | All island types | 0.555 | 1.017 |
| Area x group + isolation | Oceanic only | 0.545 | 0.916 |
| Area x group + isolation | All island types | 0.518 | 0.917 |
| Area + island type | All island types | 0.547 | 0.962 |

These are exploratory candidate comparisons, not nested validation of the model-selection procedure. Island-level validation may benefit from shared regional/archipelago history and does not establish transfer to an extinct archipelago. The prediction score uses plug-in Gaussian residual variance, matching preceding analyses. Curves show conditional median richness, not inverse area estimates; no Sunken Islands area estimates were updated. Age and the large-landmass alternative were not searched in this focused comparison.

## Interpretation

Area alone has the lowest AICc (54.76), adjusted R2 = 0.763 and LOO log RMSE = 0.551. Adding isolation is competitive (delta AICc = 1.53), but gives worse LOO RMSE (0.563). The most complex model reduces RMSE to 0.545, only about 1.1%, with delta AICc = 3.69. Prefer area alone as the parsimonious oceanic-only candidate.
Restricting training to oceanic islands does not improve point prediction on the same held-out oceanic records: area-only RMSE is 0.551 versus 0.539 using all island types, and the previously selected full-data area x group + isolation model yields 0.518. However oceanic area-only predicts uncertainty better than full-data area-only by the plug-in log-density score (0.909 versus 1.027; lower is better). Thus different metrics need not agree; there is no clear general predictive benefit from discarding the other islands. Oceanic-only remains a useful biological sensitivity analysis.
