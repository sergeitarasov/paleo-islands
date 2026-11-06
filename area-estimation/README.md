# Sunken Islands area estimation

Open `paleo-islands.Rproj` and run scripts from the repository root. Paths are relative to the project, not this folder.

## Layout

- `scripts/00-prepare-inputs.R`: reads the current island table and saved stochastic-map richness counts; extracts all paired 9 Ma draws.
- `scripts/01-pilot-regression.R`: compares the 19 pilot regression candidates using leave-one-island-out prediction.
- `results/inputs/`: derived calibration table, complete 9 Ma richness draws, frequencies, summaries and source checksums.
- `results/pilot-regression/`: fitted models, coefficients, cross-validation predictions, sensitivity checks and the original pilot data snapshot.
- `plots/pilot-regression/`: PDF and PNG figures.
- `reports/pilot-regression.md`: pilot methods, findings and limitations.
- `logs/pilot-regression/`: original execution log and R session information.

## Source data

The authoritative modern-island data remain in `data/richness/island-richness-area.xlsx`.

The authoritative map-derived richness data remain in `smaps/summaries/posterior-aic-saved-asr/richness-through-time.rds`. Counts pool states A.S, S and S.R, including observed backbone and extinct side lineages. `smaps/lineages-through-time-simple.R` documents calculation and pairing by model and posterior parameter row. We reuse that combined distribution, not a sum of separate interval endpoints. The preparation script checks that the saved combined counts equal the paired clade sums.

No large map files, MCMC chains or ASR files are duplicated here. Their uncertainty has already been propagated into the saved mapping distribution. The preparation script records checksums of its direct source files. The old regression folder has been moved here; existing fits were preserved, not recalculated during the move.

## Run in RStudio

```r
source("area-estimation/scripts/00-prepare-inputs.R")
# To refit the regression pilot against the current workbook:
source("area-estimation/scripts/01-pilot-regression.R")
```

Dependencies: readxl, MASS, ggplot2, flexmix. The first script needs only readxl. Existing outputs are replaced when scripts are rerun. The pilot's `input-data.csv` is its historical input snapshot; the preparation script's `island-calibration.csv` is extracted from the current workbook.

## Next analysis: Bayesian inverse prediction

Bayesian area estimation is complete. See `reports/bayesian-area.md` and `plots/bayesian/`. Both models, three geological assumptions, five area priors and two regression-prior widths are included. All four-chain fits passed the convergence checks; inferred area remains broad and strongly sensitive to the area prior.

Planned primary regression: Bayesian lognormal power relationship with island type. Sensitivity model: negative-binomial regression with type-specific slopes and region. Use the entire 9 Ma richness distribution, preserving its AIC-weighted map frequencies. Treat maps as alternative histories, not independent richness observations. Propagate regression uncertainty and residual predictive variation when inferring area, then average conditional area distributions over the richness frequencies.

A proposed pilot prior is uniform on log10 area from 1 to 6 (10–1,000,000 square kilometres). These bounds remain provisional and need prior sensitivity checks. Geological type and whether the target is a single island or the aggregate exposed Sunken Islands area must be specified before interpreting the result. The comparison also assumes the reconstructed two-clade richness is comparable to the calibration table; unrepresented extinct clades are not accounted for.

## Run Bayesian analysis

```r
source("area-estimation/scripts/02-fit-bayesian.R")
source("area-estimation/scripts/03-estimate-area.R")
source("area-estimation/scripts/04-report.R")
```

Additional dependency: `posterior`. Script 02 refreshes inputs, fits four chains for each model/prior setting and checks convergence. Script 03 integrates area distributions and plots them. Script 04 performs numerical checks and writes the report. `results/bayesian/area-summary-coarse-grid.csv` preserves the original coarse integration for comparison; it is not a new input dataset.

## Filtered Caribbean-versus-rest reconstruction

The additional-exclusion analysis retains 39 islands and fits separate Caribbean/non-Caribbean slopes. The new area calculation uses the non-Caribbean relationship with the complete 9 Ma richness distribution. It preserves previous estimates for comparison.

```r
source("area-estimation/scripts/10-fit-filtered-bayesian.R")
source("area-estimation/scripts/11-estimate-filtered-area.R")
```

See `reports/bayesian-filtered-area.md`, `results/bayesian-filtered/` and `plots/bayesian-filtered/`. The narrower interval is conditional on the exclusions; much of the reconstructed richness lies below the retained non-Caribbean calibration range.

## Exploratory higher-richness mixture subset

Scripts `13-fit-mixture-high-bayesian.R` and `14-estimate-mixture-high-area.R` refit a pooled power model to the 23 pink higher-richness mixture assignments and propagate the full 9 Ma richness distribution. Outputs use `bayesian-mixture-high` subdirectories; see `reports/bayesian-mixture-high.md`. This response-selected subset analysis ignores assignment uncertainty and is exploratory, not a replacement for the main calibration.

## Three-component mixture on the unfiltered data

Run `Rscript area-estimation/scripts/16-three-component-mixture.R` from the project root after script 15. This compares three versus two latent components on all 61 islands, including leave-one-out prediction and agreement with regions/island types. Results and PNG/PDF plots are in `results/three-component-mixture/` and `plots/three-component-mixture/`; see `reports/three-component-mixture.md` for interpretation and fitting instability.

## Isolation predictors

The `Island isolation` workbook tab adds coastline distance to the nearest continental mainland and nearest other landmass >=10,000 km², with names. All 61 records are retained. `scripts/17-island-distances.R` reproduces distances from archived geometry; `scripts/18-isolation-regression.R` runs the 12-model predictive comparison and creates plots. See `reports/isolation-methods.md` and `reports/isolation-regression.md`. Geometry is archived in `results/isolation/geography/` and is deliberately not downloaded again by the analysis scripts.
