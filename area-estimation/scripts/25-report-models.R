# Summarize the model screen after scripts 23 and 24.
library(ggplot2)
out <- 'area-estimation/results/model-suite'
s <- read.csv(file.path(out,'model-comparison.csv'))
i <- read.csv(file.path(out,'inverse-summary.csv'))
paired <- read.csv(file.path(out,'paired-block-comparison.csv'))
issues <- read.csv(file.path(out,'fit-issues.csv'))
labels <- c(area='Area',isolation='Area + isolation',region_type='Area + region + type',
 region_type_isolation='Area + region + type + isolation',previous_best='Area x Caribbean/rest + isolation',
 curve='Curved area',curve_adjusted='Curved area + region + type + isolation',
 block='Area + type + isolation + group RE',block_region='Area + region + type + isolation + group RE')
s$label <- unname(labels[sub('.*:','',s$id)])
s$family <- factor(s$family,levels=c('lognormal','NB','mixed'),labels=c('Lognormal','Negative binomial','Mixed lognormal'))
s$label <- paste(s$label,s$family,sep=' | ')
s$validation <- factor(s$validation,levels=c('block','island'),labels=c('Hold out geographic group','Hold out island'))
p <- ggplot(s,aes(negative_log_score,reorder(label,-negative_log_score),colour=family))+
 geom_point(size=2.5)+facet_wrap(~validation)+
 scale_colour_manual(values=c('#0072B2','#D55E00','#009E73'))+
 labs(x='Held-out negative log probability (lower is better)',y=NULL,colour=NULL)+
 theme_classic(base_size=10)+theme(legend.position='bottom')
ggsave('area-estimation/plots/model-suite/predictive-comparison.png',p,width=13,height=7,dpi=180)
ggsave('area-estimation/plots/model-suite/predictive-comparison.pdf',p,width=13,height=7)
report <- c('# Regression model comparison on the unfiltered island data','',
'All 61 islands retained. Sixteen candidates were fitted: seven lognormal and seven negative-binomial regressions, plus two lognormal mixed-effects models. No age predictor was used because the available ages are sparse and heterogeneous. No new exclusions, mixture components or Sunken Islands area estimates were introduced.','',
'## What was tested','',
'- Area only; area + mainland isolation; area + three regions + three island types; these additive effects plus isolation; the previous area-by-Caribbean/rest interaction plus isolation.','- Natural cubic spline in log area with three degrees of freedom, alone and with additive region, type and isolation. Spline knots are estimated using each training fold only.','- Random geographic-group intercept with area, type and isolation; and the same with fixed region effects. These mixed models use ML, not REML, for likelihood comparisons.','',
'## Validation','',
'Each candidate predicts every island twice: once with that island withheld, and once with its entire geographic block withheld. Fifteen blocks comprise nine multi-island groups and six singletons; the mapping is saved in geographic-blocks.csv. Random effects for unseen groups are integrated as zero-mean normals, adding group variance to residual variance. Singleton groups are also treated as unseen during island-level validation.','',
'Primary score is mean negative log probability of the observed count. For lognormal models, probabilities are integrated over the integer bin [N-0.5, N+0.5]; negative-binomial models use their probability mass. This avoids comparing a continuous density directly with a discrete probability. Zero is a permitted prediction; models are not conditioned on richness being positive. Secondary log RMSE uses lognormal medians or negative-binomial count medians (zero medians floored at 0.5 for the log metric). AICc is only comparable within likelihood families; it is not used to rank count and continuous fits together.','',
'| Model | Holdout | Negative log score | Log RMSE |', '|---|---|---:|---:|')
for(j in seq_len(nrow(s))) report <- c(report,sprintf('| %s | %s | %.3f | %.3f |',s$label[j],s$validation[j],s$negative_log_score[j],s$log_RMSE[j]))
report <- c(report,'','## Main finding','',
'The previous lognormal area x Caribbean/rest + mainland-distance model ranks first for both held-out count probability and log RMSE under both validation schemes. Additional region/type effects, negative-binomial errors, splines and group random effects did not improve these scores. The next contenders are simpler lognormal area + isolation and the simple area spline. Ranking differences are uncertain: the previous model gains only 1.83 total blocked log-score units over area + isolation (approximate block-level SE 6.11), and 3.63 over the matching negative-binomial model (SE 3.30). This is not decisive evidence for one unique best model.','',
'All 16 full-data fits and 1,216 validation fits completed, with 61 predictions per candidate per holdout scheme. One mixed-model blocked fit (region-adjusted, withholding Moluccas) had a singular random effect; this diagnostic is retained rather than dropping the fold. Full-data mixed-model random-intercept SDs were about 0.326 and 0.319 log-richness units.','',
'## Inverse area validation','',
'Six representative candidates were tested by withholding whole geographic blocks, supplying the withheld island richness and its other covariates, and estimating area. Every inverse calculation uses the same log-uniform area prior, 10–1,000,000 km2. These are approximate screening intervals: 800 parameter draws per training fit. Lognormal fixed-effect draws use the reference posterior proportional to 1/sigma2; negative-binomial draws use a joint-normal coefficient approximation and an independent log-dispersion approximation. Mixed-model fixed coefficients are sampled approximately, but variance-component uncertainty is not propagated. A new-group random intercept is integrated. These differ from the proper priors used for the final Bayesian Sunken Islands analysis.','',
'| Model | Coverage of 95% interval | Median absolute area-error factor | Median upper/lower interval ratio |','|---|---:|---:|---:|')
for(j in seq_len(nrow(i))) report <- c(report,sprintf('| %s | %.1f%% | %.2f | %.0f |',i$id[j],100*i$coverage[j],i$median_error_factor[j],i$median_interval_factor[j]))
report <- c(report,'',
'Area-only regression has the smallest median absolute inverse-area error among these six candidates (factor 2.54 versus 2.83 for the previous richness-prediction winner). Its interval coverage is 57/61 (93.4%); the previous model covers 58/61 (95.1%). Mixed effects cover 58/61 with median error factor 2.70. The ranking depends on the prediction target, so the richness-score winner should not automatically be called the best area-reconstruction model. Median interval upper/lower ratios are about 770–1,528 across candidates: uncertainty remains substantial.', '',
'Coverage is descriptive for these 61 observations, not proof of calibrated 95% intervals for Sunken Islands. Large interval widths can deliver high coverage without precise area estimation. The prior-only benchmark separates model information from coverage caused by a broad prior. Candidate selection is not nested within an outer validation loop; all rankings remain exploratory. Inventory differences, classification uncertainty and palaeogeographic uncertainty are not propagated.','',
'## Geographic grouping provenance','',
'The workbook had no archipelago field. The separate grouping is an analyst-defined broad geographic blocking scheme, not a newly verified formal archipelago database. It uses Philippines, Greater Sunda, Lesser Sunda, Moluccas, Greater Antilles, Lesser Antilles, Bahamas, Virgin Islands and Mascarenes, with singleton blocks for Langkawi, Tioman, New Caledonia, New Guinea, Sri Lanka and Madagascar. Wetar is assigned to Moluccas; alternative grouping definitions can change mixed-model and blocked-validation results. Geographic association is not proof of independent evolutionary histories.','',
'Indonesian group context: [UK Permanent Committee on Geographical Names, Indonesia toponymic factfile](https://assets.publishing.service.gov.uk/media/67852932c6428e0131881700/Indonesia_Toponymic_Factfile.pdf). Existing Matthews et al. data were consulted for available archipelago names, but not blindly name-matched: its Flores entry is the Azores island, not Indonesian Flores. The original workbook was not modified.','',
'## Reproduce','',
'From the paleo-islands project root run scripts 23-model-comparison.R, 24-inverse-validation.R and 25-report-models.R in order. The geographic-blocks.csv input is versionable and explicit. Results and figures are in model-suite subdirectories. The source workbook is read only. Existing Bayesian area estimates remain unchanged.')
writeLines(report,'area-estimation/reports/model-suite.md')
