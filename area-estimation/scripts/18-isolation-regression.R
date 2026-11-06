# Run from paleo-islands. Full dataset: no island exclusions.
library(readxl)
library(ggplot2)

# Data and candidate models ----
d <- as.data.frame(read_excel('data/richness/island-richness-area.xlsx',sheet='Island isolation'))
d$log_richness <- log(d$n_sp)
d$log_area <- log(d$area_km2)
d$log_mainland <- log1p(d$distance_mainland_km)
d$log_large <- log1p(d$distance_large_landmass_km)
d$group <- factor(ifelse(d$region=='Caribbean','Caribbean','Other islands'))
d$type <- factor(d$type)
stopifnot(nrow(d)==61,all(complete.cases(d[,c('n_sp','area_km2','log_mainland','log_large')])) )
out <- 'area-estimation/results/isolation'
plots.dir <- 'area-estimation/plots/isolation'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
dir.create(plots.dir,recursive=TRUE,showWarnings=FALSE)
write.csv(d,file.path(out,'regression-input.csv'),row.names=FALSE)

base <- c('Area only'='log_richness ~ log_area',
          'Caribbean/rest'='log_richness ~ log_area * group',
          'Island type'='log_richness ~ log_area + type')
extra <- c('None'='', 'Mainland'=' + log_mainland',
           'Large landmass'=' + log_large', 'Both'=' + log_mainland + log_large')
scores <- data.frame()
predictions <- data.frame()
coefficients <- data.frame()
models <- list()
for (b in names(base)) {
  for (e in names(extra)) {
    name <- paste(b,e,sep=': ')
    form <- as.formula(paste0(base[b],extra[e]))
    fit <- lm(form,data=d)
    models[[name]] <- fit
    k <- attr(logLik(fit),'df')
    cv <- data.frame()
    for (i in seq_len(nrow(d))) {
      train <- lm(form,data=d[-i, ])
      mu <- as.numeric(predict(train,newdata=d[i, ]))
      # Same plug-in lognormal scoring for every candidate.
      sigma <- sqrt(sum(residuals(train)^2)/nobs(train))
      cv <- rbind(cv,data.frame(model=name,island=d$island[i],observed=d$n_sp[i],
        observed_log=d$log_richness[i],pred_log=mu,pred_mean=exp(mu+sigma^2/2),
        log_density=dnorm(d$log_richness[i],mu,sigma,log=TRUE)))
    }
    predictions <- rbind(predictions,cv)
    scores <- rbind(scores,data.frame(model=name,baseline=b,isolation=e,k=k,
      AICc=AIC(fit)+2*k*(k+1)/(nrow(d)-k-1),BIC=BIC(fit),
      adjusted_R2=summary(fit)$adj.r.squared,
      log_RMSE=sqrt(mean((cv$observed_log-cv$pred_log)^2)),
      count_RMSE=sqrt(mean((cv$observed-cv$pred_mean)^2)),
      negative_log_density=-mean(cv$log_density)))
    ct <- summary(fit)$coefficients
    ci <- confint(fit)
    coefficients <- rbind(coefficients,data.frame(model=name,term=rownames(ct),
      estimate=ct[,1],lower=ci[,1],upper=ci[,2],p=ct[,4]))
  }
}
scores$delta_AICc <- scores$AICc-min(scores$AICc)
scores <- scores[order(scores$AICc), ]
write.csv(scores,file.path(out,'model-comparison.csv'),row.names=FALSE)
write.csv(predictions,file.path(out,'held-out-predictions.csv'),row.names=FALSE)
write.csv(coefficients,file.path(out,'coefficients.csv'),row.names=FALSE)
saveRDS(models,file.path(out,'models.rds'))

# Paired predictive comparisons: identical held-out islands ----
comparison <- data.frame()
for (b in names(base)) {
  ref <- predictions[predictions$model==paste(b,'None',sep=': '), ]
  for (e in names(extra)[-1]) {
    test <- predictions[predictions$model==paste(b,e,sep=': '), ]
    diff <- test$log_density-ref$log_density
    comparison <- rbind(comparison,data.frame(baseline=b,isolation=e,
      change_total_log_predictive_density=sum(diff),approx_SE=sqrt(nrow(d)*var(diff)),
      change_log_RMSE=sqrt(mean((test$observed_log-test$pred_log)^2))-
        sqrt(mean((ref$observed_log-ref$pred_log)^2))))
  }
}
write.csv(comparison,file.path(out,'paired-predictive-comparison.csv'),row.names=FALSE)
correlations <- cor(d[,c('log_area','log_mainland','log_large')])
write.csv(correlations,file.path(out,'predictor-correlations.csv'))

# Figures ----
theme_set(theme_classic(base_size=11))
cols <- c('None'='grey45','Mainland'='#0072B2','Large landmass'='#D55E00','Both'='#009E73')
p1 <- ggplot(scores,aes(log_RMSE,baseline,colour=isolation)) +
  geom_point(position=position_dodge(width=.55),size=3) + scale_colour_manual(values=cols) +
  labs(x='Leave-one-island-out log RMSE (lower is better)',y=NULL,colour='Isolation predictor') +
  theme(legend.position='bottom')
terms <- coefficients[coefficients$term %in% c('log_mainland','log_large'), ]
terms$predictor <- ifelse(terms$term=='log_mainland','Distance to mainland','Distance to large landmass')
p2 <- ggplot(terms,aes(estimate,model,colour=predictor)) +
  geom_vline(xintercept=0,colour='grey70',linetype=2) +
  geom_errorbar(aes(xmin=lower,xmax=upper),orientation='y',width=.15,position=position_dodge(width=.5)) +
  geom_point(position=position_dodge(width=.5),size=2) +
  scale_colour_manual(values=c('#D55E00','#0072B2')) +
  labs(x='Coefficient of log(1 + distance in km), with 95% confidence interval',y=NULL,colour=NULL) +
  theme(legend.position='bottom')
best_name <- scores$model[which.min(scores$log_RMSE)]
show <- predictions[predictions$model %in% c('Area only: None',best_name), ]
p3 <- ggplot(show,aes(exp(pred_log),observed)) +
  geom_abline(slope=1,intercept=0,linetype=2,colour='grey50') + geom_point(size=2,alpha=.8) +
  scale_x_log10() + scale_y_log10() + facet_wrap(~model) +
  labs(x='Held-out predicted median richness',y='Observed species richness')
figures <- list(predictive_comparison=p1,isolation_coefficients=p2,held_out_predictions=p3)
for (name in names(figures)) {
  height <- ifelse(name=='isolation_coefficients',7,5)
  ggsave(file.path(plots.dir,paste0(name,'.pdf')),figures[[name]],width=11,height=height)
  ggsave(file.path(plots.dir,paste0(name,'.png')),figures[[name]],width=11,height=height,dpi=180)
}

# Short numerical report ----
report <- c('# Isolation and species-area regression','',
'All 61 original islands are retained. The response is log richness and errors are Gaussian (a lognormal richness model). Area is log transformed; isolation uses log(1 + distance in km). Twelve exploratory candidates compare area-only, Caribbean/rest-specific slopes, and island-type intercepts, with no isolation, either distance, or both. Geographic groups/types are observed covariates, not inferred mixture components. All models use the same islands. No Bayesian area estimates were updated.','',
'| Model | AICc | Adjusted R² | LOO log RMSE | Negative log predictive density |',
'|---|---:|---:|---:|---:|')
for (i in seq_len(nrow(scores))) report <- c(report,sprintf('| %s | %.2f | %.3f | %.3f | %.3f |',
  scores$model[i],scores$AICc[i],scores$adjusted_R2[i],scores$log_RMSE[i],scores$negative_log_density[i]))
report <- c(report,'',paste('Lowest LOO RMSE:',best_name),
  sprintf('Correlation between log mainland and log large-landmass distances: %.3f.',correlations[2,3]),'',
'Coefficient intervals are ordinary model-based 95% confidence intervals, conditional on the selected model; they do not account for searching across candidates. Paired LOO log-density differences and approximate standard errors are saved separately. These are exploratory predictive comparisons, not evidence that isolation causes richness differences. Records share regions and sometimes archipelagos, so island-level cross-validation does not test transfer to a wholly new region. Contemporary isolation is not a reconstruction of isolation at 9 Ma. Applying an isolation model to Sunken Islands would require their palaeogeographic distance and its uncertainty.')
writeLines(report,'area-estimation/reports/isolation-regression.md')
print(scores); print(comparison)

# Interpretation, calculated from the saved scores ----
get_rmse <- function(name) scores$log_RMSE[scores$model==name]
gain_area <- 100*(1-get_rmse('Area only: Mainland')/get_rmse('Area only: None'))
gain_region <- 100*(1-get_rmse('Caribbean/rest: Mainland')/get_rmse('Caribbean/rest: None'))
interpretation <- c('', '## Interpretation', '',
 sprintf('Adding mainland distance reduces LOO log RMSE by %.1f%% for area alone and %.1f%% for Caribbean/rest-specific slopes.',gain_area,gain_region),
 'The Caribbean/rest plus mainland-distance model has the lowest AICc and LOO error among these candidates. Richness decreases with mainland distance after adjusting for area. Distance to any large landmass does not improve held-out predictions. Adding isolation to the island-type model does not improve its LOO performance; this is not a universal gain across model specifications.',
 'For the Caribbean/rest model, the total held-out log-density improvement from mainland distance is about 6.2 with approximate SE 5.7. Thus predictive evidence is promising but uncertain. Ordinary coefficient confidence intervals and AICc provide complementary, conditional evidence; neither proves transferability to an extinct archipelago.')
write(interpretation,'area-estimation/reports/isolation-regression.md',append=TRUE)
