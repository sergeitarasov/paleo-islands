# Posterior stochastic maps

Run from the paleo-islands project root. Use a fresh terminal R process after updating DivExtra.

## Terminal commands

First run the small pilot:

```sh
bash smaps/run-smaps-local.sh pilot
```

Then, after the pilot passes, start or resume the full analysis:

```sh
bash smaps/run-smaps-local.sh full
```

To run both sequentially, with the full analysis starting only if the pilot succeeds:

```sh
bash smaps/run-smaps-local.sh pilot && bash smaps/run-smaps-local.sh full
```

Progress is written to `runs/posterior-aic-saved-asr-pilot/run.log` or `runs/posterior-aic-saved-asr/run.log`. Keep the Mac plugged in with its lid open. The launcher prevents idle sleep, limits numeric libraries to one thread, and runs six mapping workers. It refuses a simultaneous second launcher for the same mode. After an abrupt interruption, check that all processes have stopped before removing a stale `launcher.lock` directory.

## Pilot and full design

The pilot was run successfully after commit 24fabe0: all 12 maps passed, including all three previously failing draws. Observed tip states and stored ASR priors were independently verified. Measured times were 57–72 seconds for Nesosisyphus and 17–24 seconds for Nesovinsonia. The full analysis remains stopped for terminal launch.

The pilot freshly generates **12 maps**: six parameter vectors, each mapped on both clades. They are gr6 maps 1, 2 and 20, and gr7 maps 1, 18 and 46. This includes all three previously failing draws, with their original posterior rows and seeds. Pilot outputs are separate and are not added to the full sample.

A successful pilot writes `runs/posterior-aic-saved-asr-pilot/pilot-success.rds` and `summaries/posterior-aic-saved-asr-pilot/pilot-results.csv`. The full runner checks that the installed package and scientific inputs still match this pilot.

The full analysis produces **2,000 maps per clade**: 889 gr6 and 1,111 gr7, using the saved two-model AIC weights. Complete posterior rate vectors come from four chains per model, removing the first 1,000 of 5,000 sweeps. Each clade pair shares its model and posterior vector, but uses separate mapping seeds. The full run resumes existing validated map files; failed draws are never replaced by different draws.

The original 9 Ma epoch boundary, extinction-only direct sampler and complete focal sampling are retained. Hidden side lineages are expanded explicitly. The observed backbone is conditioned on all supplied clade tips; focus selection then retains the Nesovinsonia path or four Nesosisyphus paths and side clades born on them.

## Root priors and numerical checks

Fixed ML-ASR root priors are read directly from `asr/results/sisyphini.gr6.rds`, `sisyphini.gr7.rds`, `nesovin.gr6.rds` and `nesovin.gr7.rds`, identifying roots by topology and aligning state names. They are not recalculated per posterior draw. This is the requested two-stage plug-in analysis: the subclade likelihood updates the supplied prior. It is not a joint full-tree posterior and does not propagate uncertainty in the root-prior vectors.

DivExtra commit **24fabe0** fixes rounded interpolation probes at epoch boundaries without weakening accuracy thresholds. The installed cache function has been compared with this source, and the analytic regression test passes. The numerical investigation is preserved under `runs/posterior-aic-saved-asr/interpolation-fix/`. Previously completed maps remain usable because their caches passed the unchanged thresholds; the runner accepts only the explicitly validated prior package fingerprint or the current installation.

Every new map is checked for extinct generated tips, correct retained observed tips and mapped branch durations. Files are written atomically. Errors are logged per worker, and final collection requires every worker to succeed. The package and input hashes are saved with each map.

## Files

- `setup-smaps.R`: data, posterior draws, model weights and saved ASR root priors.
- `run-smaps-posterior-parallel.R`: pilot/full modes, six workers and collection.
- `plot-sunken-richness-posterior.R`: richness summaries and figures after full collection.
- `maps/<run-name>/`: stochastic maps only.
- `runs/<run-name>/`: logs, priors, draw selections and provenance.
- `plots/<run-name>/`: PDF and PNG figures.
- `summaries/<run-name>/`: tables, count matrices and reports.
- `examples/make-smaps.R`: your original ML example.

The full run saves separate model pools, AIC-weighted combined pools and observed-only copies for arrival-time analyses. Richness plots count all living lineages in A.S, S or S.R, including observed backbone and extinct side lineages, on a 0.1 Ma grid. Matched clade counts are added before calculating means, medians and pointwise 95% distribution intervals (quantile type 7). Combined curves cover only ages represented by both clades, 0–22.73 Ma. These intervals describe variation across histories, not Monte Carlo uncertainty of the mean.
