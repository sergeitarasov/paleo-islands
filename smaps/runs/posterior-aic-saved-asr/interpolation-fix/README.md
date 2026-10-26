# Extinction interpolation-check correction

The original failures were numerical false positives near the 9 Ma epoch boundary. The ODE evaluated rounded, representable probe ages, but the check computed interpolation at nominal fractions 0.25, 0.5 and 0.75. At very short intervals, these no longer describe the same time. Refining the grid reached machine precision without eliminating the apparent error.

The DivExtra correction calculates each interpolation weight as `(probe - left) / width`. It retains the original ODE tolerances, the 1e-6 interpolation threshold and the 30-refinement limit. It does not floor extinction probabilities or substitute posterior draws.

## Validation

- Both original Nesosisyphus failures (gr6 map 20, posterior row 5645; gr7 map 18, row 1775) passed with the corrected check using default, tighter and still-tighter ODE tolerances.
- The maximum absolute difference between default and tightest E(t) calculations over 1,100 sampled ages was 3.50e-14 for gr6 and 2.15e-14 for gr7.
- A single-state pure-death model with zero recent extinction and old extinction 0.0155, separated at 9 Ma, reproduces the failure before the fix. The new regression test compares the corrected cache with its analytic solution.
- Installed-code regression tests passed: rounded-probe regression, direct extinction, tree-conditioned mapping and no-side-expansion invariance. See tests.log.
- An additional draw (gr7 map 46, posterior row 1237) had failed with the same message before the restart. It is retried by its original worker.

Completed maps from the old installation retain their original package hashes: their caches already passed the unchanged numerical checks. The runner explicitly accepts that one validated previous installation and the new installation. Scientific input hashes must still match. Prior logs and error files are preserved; no failed posterior draws were replaced.

All six workers were resumed. The original failed maps must finish and pass structural validation before their recovery is considered complete.

## Confirmed recovery

All three previously failed Nesosisyphus maps completed and passed map validation after the fix: gr6 map 20 (68.062 seconds), gr7 map 18 (65.337 seconds), and gr7 map 46 (66.996 seconds). Their original posterior rows and seeds are preserved. The full six-worker run continues.
