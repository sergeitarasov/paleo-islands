# Run from the paleo-islands project root.
library(readxl)
library(flexmix)
library(ggplot2)
library(ggrepel)
d <- as.data.frame(read_excel("data/richness/island-richness-area.xlsx"))
d$x <- log(d$area_km2/10000)
d$y <- log(d$n_sp)
d$group <- factor(ifelse(d$region == "Caribbean", "Caribbean", "Other islands"))
stopifnot(nrow(d)==61)
out <- "area-estimation/results/unfiltered-mixture"
plots.dir <- "area-estimation/plots/unfiltered-mixture"
dir.create(out,recursive=TRUE,showWarnings=FALSE)
dir.create(plots.dir,recursive=TRUE,showWarnings=FALSE)
write.csv(d,file.path(out,"input-data.csv"),row.names=FALSE)

# Fit known groups and pooled model ----
fits <- list(Pooled=lm(y~x,data=d),
  "Regions: shared slope"=lm(y~x+group,data=d),
  "Regions: separate slopes"=lm(y~x*group,data=d))

# Two latent power regressions; common variance prevents variance collapse.
# Region is not supplied to the mixture. Minimum mixing proportion = 0.1.
driver <- FLXMRglmfix(family="gaussian",varFix=TRUE)
set.seed(261028)
starts <- data.frame()
best <- NULL
for (i in 1:200) {
  m <- flexmix(y~x,data=d,k=2,model=driver,
    control=list(iter.max=2000,minprior=.1,tolerance=1e-8))
  starts <- rbind(starts,data.frame(start=i,logLik=as.numeric(logLik(m)),
    components=m@k,converged=m@converged))
  if (m@converged && (is.null(best) || logLik(m)>logLik(best))) best <- m
}
stopifnot(!is.null(best),best@k==2)
fits[["Latent mixture"]] <- best
write.csv(starts,file.path(out,"mixture-starts.csv"),row.names=FALSE)
saveRDS(fits,file.path(out,"models.rds"))

# Leave one island out. Latent membership must NOT use held-out richness.
scores <- data.frame(); predictions <- data.frame()
for (name in names(fits)) {
  f <- fits[[name]]
  k <- attr(logLik(f),"df"); n <- nrow(d)
  row <- data.frame(model=name,k=k,AICc=AIC(f)+2*k*(k+1)/(n-k-1),BIC=BIC(f))
  for (i in seq_len(n)) {
    if (name != "Latent mixture") {
      cv <- lm(formula(f),data=d[-i, ])
      mu <- as.numeric(predict(cv,d[i, ,drop=FALSE]))
      sd_ml <- sqrt(mean(residuals(cv)^2))
      pred_log <- mu
      pred_mean <- exp(mu+sd_ml^2/2)
      log_density <- dnorm(d$y[i],mu,sd_ml,log=TRUE)
      converged <- TRUE; components <- 1
    } else {
      set.seed(261028+i)
      cv <- stepFlexmix(y~x,data=d[-i, ],k=2,nrep=30,model=driver,
        control=list(iter.max=2000,minprior=.1,tolerance=1e-8),verbose=FALSE)
      if (inherits(cv,"stepFlexmix")) cv <- getModel(cv,"BIC")
      # Nearly coincident components can make EM converge slowly.
      # Continue the selected fit at the SAME tolerance, without new random starts.
      if (!cv@converged) cv <- flexmix(y~x,data=d[-i, ],cluster=posterior(cv),model=driver,
        control=list(iter.max=20000,minprior=.1,tolerance=1e-8))
      stopifnot(cv@converged)
      pars <- parameters(cv)
      mu <- pars[1, ]+pars[2, ]*d$x[i]
      sd_ml <- pars[3, ]; weights <- prior(cv)
      pred_log <- sum(weights*mu)
      pred_mean <- sum(weights*exp(mu+sd_ml^2/2))
      terms <- log(weights)+dnorm(d$y[i],mu,sd_ml,log=TRUE)
      log_density <- max(terms)+log(sum(exp(terms-max(terms))))
      converged <- cv@converged; components <- cv@k
    }
    predictions <- rbind(predictions,data.frame(model=name,island=d$island[i],observed=d$n_sp[i],
      pred_log=pred_log,pred_mean=pred_mean,log_density=log_density,converged=converged,components=components))
  }
  p <- predictions[predictions$model==name, ]
  row$log_RMSE <- sqrt(mean((log(p$observed)-p$pred_log)^2))
  row$count_RMSE <- sqrt(mean((p$observed-p$pred_mean)^2))
  row$mean_negative_log_density <- -mean(p$log_density)
  scores <- rbind(scores,row)
}
scores$delta_AICc <- scores$AICc-min(scores$AICc)
write.csv(scores,file.path(out,"model-comparison.csv"),row.names=FALSE)
write.csv(predictions,file.path(out,"held-out-predictions.csv"),row.names=FALSE)

# Order component names by fitted richness at 10,000 km2.
pars <- parameters(best)
order_components <- order(pars[1, ])
prob <- posterior(best)[,order_components]
pars <- pars[,order_components]
components <- data.frame(component=c("Lower-richness component","Higher-richness component"),
  intercept=pars[1, ],area_exponent=pars[2, ],sigma=pars[3, ],weight=prior(best)[order_components])
write.csv(components,file.path(out,"component-parameters.csv"),row.names=FALSE)
d$p_high <- prob[,2]
d$component <- ifelse(d$p_high>.5,"Higher-richness component","Lower-richness component")
d$confidence <- pmax(d$p_high,1-d$p_high)
write.csv(d,file.path(out,"island-membership.csv"),row.names=FALSE)
counts <- table(d$group,d$component)
write.csv(as.data.frame(counts),file.path(out,"group-confusion.csv"),row.names=FALSE)
write.csv(as.data.frame(table(d$region,d$component)),file.path(out,"region-confusion.csv"),row.names=FALSE)
# Label-invariant assignment accuracy and adjusted Rand index.
accuracy <- max(sum(diag(counts)),sum(diag(counts[,2:1])))/nrow(d)
choose2 <- function(x) x*(x-1)/2
index <- sum(choose2(counts))
expected <- sum(choose2(rowSums(counts)))*sum(choose2(colSums(counts)))/choose2(sum(counts))
maximum <- (sum(choose2(rowSums(counts)))+sum(choose2(colSums(counts))))/2
ari <- (index-expected)/(maximum-expected)
agreement <- data.frame(best_matching_accuracy=accuracy,adjusted_Rand_index=ari,
  uncertain_below_0.8=sum(d$confidence<.8),uncertain_below_0.6=sum(d$confidence<.6),
  best_start_repeats=sum(abs(starts$logLik-as.numeric(logLik(best)))<1e-4 & starts$converged))
write.csv(agreement,file.path(out,"membership-agreement.csv"),row.names=FALSE)

# Plots ----
theme_set(theme_classic(base_size=12))
colors <- c("Lower-richness component"="#0072B2","Higher-richness component"="#CC79A7")
grid <- data.frame()
for (j in 1:2) {
  x <- seq(min(d$x),max(d$x),length.out=300)
  grid <- rbind(grid,data.frame(area_km2=10000*exp(x),richness=exp(pars[1,j]+pars[2,j]*x),component=components$component[j]))
}
p1 <- ggplot(d,aes(area_km2,n_sp)) +
  geom_line(data=grid,aes(area_km2,richness,colour=component),linewidth=.8) +
  geom_point(aes(colour=component,shape=group),size=2.5) +
  geom_text_repel(aes(label=island,colour=component),seed=26,max.overlaps=Inf,max.time=10,size=3,show.legend=FALSE) +
  scale_x_log10() + scale_y_log10() + scale_colour_manual(values=colors) +
  labs(x=expression("Island area (km"^2*")"),y="Species richness",colour=NULL,shape=NULL,
    caption="Component curves nearly overlap; every island has membership confidence below 80%.") +
  theme(legend.position="bottom",legend.box="vertical",plot.caption=element_text(size=9))
d$island_order <- factor(d$island,levels=d$island[order(d$group,d$p_high)])
p2 <- ggplot(d,aes(p_high,island_order,colour=group)) + geom_vline(xintercept=.5,linetype=2,colour="grey65") +
  geom_segment(aes(x=0,xend=p_high,yend=island_order),colour="grey85") + geom_point(size=2.4) +
  scale_x_continuous(limits=c(0,1)) + scale_colour_manual(values=c(Caribbean="#D55E00","Other islands"="#0072B2")) +
  labs(x="Probability of higher-richness component",y=NULL,colour=NULL) + theme(legend.position="bottom")
metrics <- rbind(data.frame(model=scores$model,metric="Log-scale RMSE",value=scores$log_RMSE),
  data.frame(model=scores$model,metric="Negative log predictive density",value=scores$mean_negative_log_density))
p3 <- ggplot(metrics,aes(value,model)) + geom_point(size=3,colour="#0072B2") + facet_wrap(~metric,scales="free_x") +
  labs(x="Leave-one-out error (lower is better)",y=NULL)
ggsave(file.path(plots.dir,"mixture-labelled.pdf"),p1,width=12,height=9)
ggsave(file.path(plots.dir,"mixture-labelled.png"),p1,width=12,height=9,dpi=180)
ggsave(file.path(plots.dir,"membership-by-island.pdf"),p2,width=8,height=16)
ggsave(file.path(plots.dir,"membership-by-island.png"),p2,width=8,height=16,dpi=180)
ggsave(file.path(plots.dir,"prediction-comparison.pdf"),p3,width=11,height=5)
ggsave(file.path(plots.dir,"prediction-comparison.png"),p3,width=11,height=5,dpi=180)

report <- c("# Latent mixture versus Caribbean/rest regression", "",
"Dataset: all 61 original records (16 Caribbean, 45 other), with no exclusions. Both mixture components are power regressions with Gaussian log-richness residuals and a common variance, matching the grouped regression variance assumption. No region or island-type labels enter the mixture. Mixing weights are constant across area. A minimum component-weight setting of 0.1 and shared variance avoid tiny-variance solutions. Thus this is a restricted two-component mixture, not a search across every possible mixture.", "",
"200 independent starts were tried for the full mixture; 30 per leave-one-out fold. Known-region models use the held-out island's known region; mixture predictions use the fitted mixing weights and never use held-out richness to choose a component. Log predictive density evaluates the full mixture distribution, including its possible multimodality. All models use ML residual variance for cross-validation density and arithmetic-mean correction; count RMSE can differ slightly from previous scripts that used residual-df variance.","",
"| Model | AICc | BIC | Log CV RMSE | Negative log predictive density |", "|---|---:|---:|---:|---:|")
for(i in seq_len(nrow(scores))) report <- c(report,sprintf("| %s | %.2f | %.2f | %.3f | %.3f |",scores$model[i],scores$AICc[i],scores$BIC[i],scores$log_RMSE[i],scores$mean_negative_log_density[i]))
report <- c(report,"",
"## Interpretation", "",
"The two-component mixture is not supported on the unfiltered data. AICc is 7.11 worse than pooled and 7.16 worse than the shared-slope regional model. Its log-scale prediction error (0.966) is almost identical to pooled (0.969), with slightly worse predictive log density. The fitted log-likelihood gain over pooled is only 0.010. Component curves largely coincide.", "",
"Every island has assignment confidence below 0.8; 34/61 are below 0.6. Hard counts of 6 versus 55 therefore should not be interpreted as biological classes. Adjusted Rand agreement with geography is -0.038. The apparently clearer mixture split in the 39-island subset is not robust to restoring excluded data.", "",
"Ten of 200 full-data starts hit the 2,000-iteration limit; the best converged fit was selected. All 200 final likelihoods were within 0.1 log units, consistent with a very flat, poorly separated mixture. Where a cross-validation fit hit its iteration limit, the selected fit was continued for up to 20,000 iterations at the same 1e-8 tolerance. All final holdout fits converged. Exact component memberships remain weakly identified despite these small likelihood differences.", "",
sprintf("Best label-matched agreement with Caribbean/rest: %.1f%%. Adjusted Rand index: %.3f (1 means identical partitions; 0 is chance-adjusted agreement). %d islands have membership confidence below 0.8, including %d below 0.6. Best full-data likelihood reproduced in %d/200 starts.",100*accuracy,ari,agreement$uncertain_below_0.8,agreement$uncertain_below_0.6,agreement$best_start_repeats),"",
"Hard membership is used only for summaries and point colours. Parameter estimation uses soft probabilities. Component names refer to richness at 10,000 km²; component lines may cross elsewhere. Membership-region agreement is descriptive, evaluated on the fitting data. It does not establish independently validated regional classes. No records were excluded in this run. The fitted dataset is the original workbook. No Bayesian area estimates have been updated.","",
"Outputs: model-comparison.csv, group/region-confusion.csv, component-parameters.csv, island-membership.csv, membership-agreement.csv, mixture-starts.csv and complete held-out predictions. Reproduce with area-estimation/scripts/15-unfiltered-mixture.R from the project root.")
writeLines(report,"area-estimation/reports/unfiltered-mixture.md")
print(scores); print(counts); print(agreement); print(components)
