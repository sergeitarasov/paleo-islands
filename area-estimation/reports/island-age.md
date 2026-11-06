# Island-age compilation: initial source audit

The `Island isolation` worksheet contains age information for records currently classified as Oceanic. The original richness worksheet, island classifications, richness, areas and isolation values are unchanged. This is a partial literature compilation, not a complete or harmonized emergence-age database. No new age regression has been fitted.

## Coverage and definitions

Of 31 oceanic records, seven have published numerical point values and three have one-sided bounds. Twenty-one have no numerical value verified here. Thirty non-oceanic records are marked `not_applicable` for this compilation. Missing values are not zeros.

| Island | Reported age (Ma) | Interpretation |
|---|---:|---|
| Mauritius | 8.9 | Geological age adopted by Valente et al.; exposed-rock-age proxy |
| Réunion | 5 | Geological construction age adopted by Valente et al.; not necessarily emergence |
| Lombok | 11 | Appearance in a geological reconstruction |
| Sumbawa | 11 | Same reconstruction; conflicting/younger volcanic history explicitly flagged |
| Flores | 9 | Appearance in a geological reconstruction |
| Saint Lucia | 18.3 | Older-arc volcanic age; emergence not verified |
| Saint Vincent | 3.5 | Younger-arc volcanic age; emergence not verified |
| Dominica | <5 | Published upper bound |
| Grenada | >5 | Published lower bound; volcanic-age caveat |
| Martinique | <34 | Bound for composite island with differently aged precursors |

Bounds are **not confidence intervals**. Blank lower or upper bounds mean no bound was supplied. Inequality direction and strictness are recorded in notes. Point values have no inferred uncertainty interval. The `island_age_ma` column records reported estimates/proxies; its name does not make them interchangeable emergence dates.

For Guadeloupe, the published 4–6 Ma estimate concerns Basse Terre. It was not assigned to the entire archipelago represented in the richness table. Similarly, no Philippine archipelago-wide geological age was assigned to individual islands. Ages inferred from beetle or bird divergences were not used as geological covariates.

Unresolved numerical ages: Adonara, Basilan, Bohol, Cebu, Halmahera, Komodo, Leyte, Luzon, Marinduque, Mindanao, Negros, Rinca, Samar, Ternate, Wetar, Eleuthera, Guadeloupe, New Providence, Saint John, Saint Thomas and Tortola. This means no suitable value was verified in this source audit, not that no relevant literature exists.

## Verified sources

1. [Valente et al. (2020), Nature](https://doi.org/10.1038/s41586-020-2022-5), Extended Data Table 1 and Supplementary Data 2 (`41586_2020_2022_MOESM4_ESM.xlsx`). Mauritius: Archipelagos E28 = 8.9; Réunion: E36 = 5. Extended table references 81 (Montaggioni & Nativel 1988) and 89 (Gillot et al. 1994), respectively. Values checked against the published table image and spreadsheet. This is an archipelago-level dataset; only geographically matching records were used.
2. [Tänzler et al. (2016), Scientific Reports](https://doi.org/10.1038/srep18793), Discussion. Reports the Spakman & Hall (2010) reconstruction: Lombok and Sumbawa appear at 11 Ma, Flores at 9 Ma. The same discussion describes much younger Sumbawa volcanism (Barberi et al. 1987). These values are model-derived geological proxies, not precise dated emergence events.
3. [Thorpe et al. (2018), Molecular Phylogenetics and Evolution](https://doi.org/10.1016/j.ympev.2018.06.014), Table 1, printed pp. 689–690. Island-age column only; molecular ages in adjacent columns were not extracted. Table and superscripts visually checked: Saint Lucia is 18.3 with footnote 6, not 18.36. Footnotes warn about precursor islands, arc superposition and uncertain emergence. Original geological references are retained in workbook notes.

## Other suggested sources screened

- Matthews et al. (2026), [10.1073/pnas.2518902123](https://doi.org/10.1073/pnas.2518902123): downloaded island table has area, geography, type and isolation but no island-age column.
- Weigelt et al. (2016), [10.1038/nature17443](https://doi.org/10.1038/nature17443): downloaded supplementary Excel has historical area/isolation and environmental covariates but no general island-age column.
- Weigelt et al. (2013), [10.1073/pnas.1306309110](https://doi.org/10.1073/pnas.1306309110): useful physical/environmental database; no ages were imported from it. Dryad file download was inaccessible during this audit.
- Tonelli et al. (2023), [10.1111/jbi.14715](https://doi.org/10.1111/jbi.14715): Mediterranean/Macaronesian scope does not supply matching records for our oceanic-island subset. Dryad metadata were inspected; the spreadsheet download was inaccessible.
- Whittaker et al. (2017), [10.1126/science.aam8326](https://doi.org/10.1126/science.aam8326): conceptual geological-history reference; no numerical ages imported.

## Before regression

The desired biological covariate is preferably time since establishment of persistent subaerial habitat, not the age of underlying crust or the youngest eruption. Current coverage and mixed definitions are insufficient for a defensible general age-adjusted regression. Prioritize harmonizing emergence definitions and resolving missing records before fitting. Do not impute non-oceanic islands as age zero or substitute bounds as exact values. If enough compatible ages become available, compare area-only and area-plus-isolation/age models on the **same complete-case islands**, and examine age uncertainty separately. This avoids confusing a change in sample composition with improved model fit.

The workbook values are also saved in `area-estimation/results/island-age/island-ages.csv`; a JSON copy is provided for exact reapplication. Source-table snapshots are archived alongside them. Isolation-only regression results predate the request to add age; they remain exploratory and have not been expanded with these heterogeneous age proxies.
