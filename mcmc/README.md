# Empirical MCMC: Afr–Sunken–Mauritius models

Open the **paleo-islands project root** as your working directory.
This folder contains the completed group-6 and group-7 MCMC run and the inputs
needed to reproduce it, without relying on the old OneDrive project.

## Layout

- `mcmc-simple-gr6.R`, `mcmc-simple-gr7.R`: readable, standalone run scripts.
- `inputs/`: exact tree, geographic observations, model YAML files, saved ML
  fits and fixed pilot transformations used for this analysis.
- `results/roihu-5000-20261023/`: eight completed chain files from Roihu.
- `provenance/original-scripts/`: unmodified scripts used for the submitted run.
- `provenance/roihu/`: LOCAL ONLY, excluded from Git. Original submission/check scripts, logs, package-session
  information and local-versus-remote likelihood validation. These archived
  scripts retain the old project paths; they are records, not new entry points.
- `download-roihu.sh`: optional re-download into this new project's layout.

## Completed analysis

Four independent chains per model, 5,000 slice sweeps per chain, first 1,000
discarded for summaries: 16,000 retained joint draws per model. Every sweep is
saved. Group 6 has nine free rates and group 7 ten. Priors are independent
exponentials with mean 1 on free rates; the log-rate Jacobian is included.
The pilot-based fixed transformation changes sampling coordinates, not the
likelihood or prior. Group 6's transformation used matching parameter columns
from the group-7 pilot and a group-6 ML-based center.

The original full tree, sampling fractions, African root, two-epoch model
and 9 Ma boundary are preserved. Validation job 2499294 and all eight tasks of
2499295 completed successfully. Maximum local-versus-Roihu validation difference
was approximately 9.1e-13. Runtime was 20.7–22.4 hours per group-6 chain and
23.6–27.4 hours per group-7 chain. All eight result matrices have 5,000 rows
and finite positive rates.

Convergence was assessed on 25 October 2026: all 9 group-6 and all 10 group-7
parameters passed the diagnostic screen. Maximum R-hat was 1.00143 and 1.00104;
minimum bulk/tail ESS was 2684/1271 and 3980/3501, respectively.
Numerical diagnostics support adequate mixing but do not prove convergence or
validate the model assumptions. See `summaries/convergence-report.txt`.

## Analyze and plot

Install `posterior` once if needed: `install.packages("posterior")`.
Then run from the project root:

```r
source("mcmc/analyse-mcmc.R")
```

This readable script makes separate posterior density and post-warmup trace
figures for each model, in vector PDF and 300 dpi PNG, under `mcmc/plots/`.
No new MCMC is run. Summary tables in `mcmc/summaries/` include means, medians,
equal-tailed 95% credible intervals, ML, R-hat, ESS and Monte Carlo errors.
The shaded density interval is based on raw posterior quantiles; the KDE is
only a visual aid. These are posterior credible intervals, not ML confidence
intervals. The figures retain exact parameter identifiers and omit an overall
title so the interpretation can go in the publication caption.

## Read a chain

```r
x <- readRDS("mcmc/results/roihu-5000-20261023/gr6-chain-1.rds")
rates <- x$rates[1001:5000, ]
```

Each row is a complete joint parameter vector. The RDS also contains warmup,
sampler coordinates, runtime, model/chain IDs, seed, prior rate, transformation,
  ML estimates and sanitized session information. Keep model-specific summaries; do not
treat group-6 and group-7 parameter matrices as interchangeable.

## Rerun, if wanted

Requires `divextra`, `diversitree` and their dependencies. See the recorded
session information for versions. No analysis is launched simply by opening
this folder.

```sh
Rscript --vanilla mcmc/mcmc-simple-gr6.R
Rscript --vanilla mcmc/mcmc-simple-gr7.R
```

Each script runs four chains sequentially. To run a single chain:

```sh
MCMC_CHAIN=1 Rscript --vanilla mcmc/mcmc-simple-gr6.R
```

Reruns save to `mcmc/results/local-rerun/`, protecting the completed Roihu
results. Repeating a rerun overwrites its own matching files. The simple
scripts save after each chain, not midway through it, and do not auto-resume.

## Retrieve from Roihu again

```sh
bash mcmc/download-roihu.sh
```

Requires your CSC access and valid SSH certificate. Set `ROIHU_HOST` and
`ROIHU_REMOTE` in your terminal; there are no personal defaults in the script.
Keep authentication in your private SSH configuration or `RSYNC_RSH` environment
variable. No private keys or certificates are included in this folder.
Re-downloads go into the Git-ignored `provenance/roihu/raw-results/` directory
because original RDS session metadata includes machine paths. They do not replace
the sanitized publication copies in `results/`.

## Provenance of this copy

Copied from the `divextra-sunken` project's `asr` outputs on 25 October 2026.
Paths, the rerun output directory and saved session metadata were adjusted in
the two top-level R scripts; scientific settings and sampling code are unchanged.
Input snapshots intentionally duplicate some project data to preserve the
exact completed analysis. Future data changes should start a separate run.
Original files were not moved or deleted. No Git commit or push was made.

## Privacy

Shared chain objects retain R/package versions but omit full session metadata
containing personal or system installation paths. All other chain fields,
including posterior draws, are unchanged. Raw cluster provenance stays local
and is excluded by the repository's `.gitignore`; do not force-add it. These
rules do not protect files manually attached elsewhere. No passwords or private
keys were found in this MCMC collection. This is not a guarantee about unrelated
datasets or future files added to the repository.
