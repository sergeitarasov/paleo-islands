# Requested 27-island subset

Included all 16 Caribbean records plus Réunion, Mauritius, Madagascar, New Caledonia, New Zealand, Borneo, Java, Sumatra, Sulawesi, New Guinea and Sri Lanka. Mauritius retains the workbook value of 6 species. NZ uses the supplied 18 species and 263,310 km2 land area from https://www.dfat.gov.au/geo/new-zealand/new-zealand-country-brief. NZ is a country/island-group observation; area source is not yet harmonized with Global Islands v3. The workbook is unchanged.

Exploratory lognormal regression. AICc comparisons only use models fitted to the same 27 observations. LOO validates islands, not entire geographic regions. The eleven non-Caribbean observations cover different climates/geological histories and do not define a single homogeneous reference group.

             model  n        R2 adjusted_R2     AICc LOO_log_RMSE
1           pooled 27 0.8381566   0.8316829 69.59156    0.8251855
2 group_intercepts 27 0.8592469   0.8475174 68.59649    0.8315027
3     group_slopes 27 0.8907288   0.8764761 64.79973    0.7904459

               subset  n     slope     lower     upper        R2            p
x                 all 27 0.5668834 0.4642761 0.6694908 0.8381566 2.226358e-11
x1 without_Madagascar 26 0.5484591 0.4424161 0.6545021 0.8260196 1.353209e-10
x2          Caribbean 16 0.3578301 0.2398087 0.4758515 0.7512719 1.395374e-05
x3       other_eleven 11 0.6815364 0.3826103 0.9804625 0.7471972 5.970087e-04


	Spearman's rank correlation rho

data:  d$area_km2 and d$n_sp
S = 359.17, p-value = 5.025e-10
alternative hypothesis: true rho is not equal to 0
sample estimates:
      rho 
0.8903646 


A stronger association within this selected subset does not establish that it is the correct calibration for Sunken Islands. Selection and inventory uncertainty are not propagated, and forward association alone does not establish accurate inverse area recovery. No Sunken Islands area estimate was recalculated.
