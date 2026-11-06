# Rectangular biogeographic ingroup tree

The vector PDF shows the complete 271-tip ingroup from `asr/results/full.avg.rds` as a conventional rectangular timetree. All original specimen identifiers are displayed, including STL, ST, COL and UCE identifiers where supplied. The page is 17 x 40 inches to keep all names readable; it is intended for a zoomable supplementary figure, not reduction onto one A4 page.

Run from the GitHub repository root:

```sh
Rscript asr/plot-asr-rectangular.R
```

Required packages: `ape` and `readxl`. The script also works when invoked by its absolute path. Optional command-line arguments specify the project directory and a separate output directory.

## Approved names

The author approved these changes on 2026-11-03:

| Specimen | Original name | Displayed name |
|---|---|---|
| STL18 | Sisyphini sp morpho1 | Sisyphus sp1 |
| STL23 | Sisyphini sp morpho6 | Sisyphus sp6 |
| STL40 | Sisyphini sp morpho8 | Sisyphus sp8 |
| STL42 | Sisyphini sp morpho 10 | Sisyphus sp10 |
| STL9 | Sisyphus gazanu s | Sisyphus gazanus |
| STL10223 | Tanzanolus sp1 | Tanzanolus leleupi |

The separately approved specimen STL10224 (Tanzanolus sp2 -> Tanzanolus leleupi) is absent from this 271-tip dataset; it is included in the updated 316-tip age-bar figure. Changes in this figure affect labels only. Other morphospecies names are preserved. The original workbook and saved reconstruction are unchanged.

## Reading the figure

- Internal pies show the saved AIC-weighted marginal ancestral-range probabilities (gr6 weight 0.444515; gr7 weight 0.555485).
- Tip circles show observed geographic states.
- All 16 model states are retained. Joint states have distinct colours and are not pooled.
- The palette matches the current `plot-asr-fan.R` and `plot-asr-clades-combined.R` scripts. S denotes paleoislands and R denotes Mauritius.
- The horizontal axis is age before present in Ma. The saved ingroup root age is approximately 96.98 Ma.
- Pie radii are slightly reduced where needed to prevent overlapping symbols. Probabilities are not thresholded or recalculated.

## Outputs

- `full-asr-rectangular.pdf`: complete single-page vector figure.
- `tip-labels.csv`: original names, approved display labels, specimen IDs, observed states and plotted rows.
- `node-probabilities.csv`: unchanged saved probabilities, node ages, and stable node tags.
- `full-asr-rectangular-Ma.nex` and `.tre`: renamed tree, branch lengths in Ma. Internal `ASR_node_*` labels connect nodes to the probability table; these are identifiers, not support values. Plain Newick labels are quoted to preserve parentheses; `ape::read.tree()` retains these outer quotes when reading, so remove them before matching labels if necessary.
- `figure-legend.tex`: legend text ready to paste into the supplement.
- `validation.txt`, `input-checksums.csv`, `session-info.txt`, `node-pie-sizes.csv`: checks and reproducibility information.

For LaTeX, put the legend on a normal page and insert the large figure with `pdfpages`:

```latex
% Preamble: \usepackage{pdfpages}
\input{figs-suppl/figure-legend.tex}
\includepdf[pages=1,fitpaper=true]{figs-suppl/full-asr-rectangular.pdf}
```

Copy the PDF and legend into your LaTeX project's `figs-suppl/` directory, or adjust these paths.
