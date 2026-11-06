# Requested 21-island subset

Included all 16 Caribbean records plus Réunion, Mauritius, Madagascar, New Caledonia and New Zealand. Mauritius retains the workbook value of 6 species. NZ uses the supplied 18 species and 263,310 km2 land area from https://www.dfat.gov.au/geo/new-zealand/new-zealand-country-brief. NZ is a country/island-group observation; area source is not yet harmonized with Global Islands v3. The workbook is unchanged.

Exploratory lognormal regression. AICc comparisons only use models fitted to the same 21 observations. LOO validates islands, not entire geographic regions. The five non-Caribbean observations cover different climates/geological histories and do not define a single homogeneous reference group.

             model  n        R2 adjusted_R2     AICc LOO_log_RMSE
1           pooled 21 0.7282620   0.7139600 55.88154    0.8778554
2 group_intercepts 21 0.7434560   0.7149511 57.76147    0.9435128
3     group_slopes 21 0.7815643   0.7430168 57.88449    1.0184845

               subset  n     slope      lower     upper        R2            p
x                 all 21 0.4724615  0.3338832 0.6110398 0.7282620 8.764095e-07
x1 without_Madagascar 20 0.3917391  0.2547258 0.5287524 0.6671709 1.111384e-05
x2          Caribbean 16 0.3578301  0.2398087 0.4758515 0.7512719 1.395374e-05
x3         other_five  5 0.6478092 -0.2428723 1.5384907 0.6410460 1.035892e-01


	Spearman's rank correlation rho

data:  d$area_km2 and d$n_sp
S = 315.92, p-value = 1.671e-05
alternative hypothesis: true rho is not equal to 0
sample estimates:
      rho 
0.7948589 


A stronger association within this selected subset does not establish that it is the correct calibration for Sunken Islands. Selection and inventory uncertainty are not propagated, and forward association alone does not establish accurate inverse area recovery. No Sunken Islands area estimate was recalculated.
