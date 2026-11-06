# Original 61-island dataset; no exclusions. Run from the project root.
library(readxl)
library(flexmix)
library(ggplot2)
library(ggrepel)
d <- as.data.frame(read_excel("data/richness/island-richness-area.xlsx"))
d$x <- log(d$area_km2/10000); d$y <- log(d$n_sp)
d$group <- ifelse(d$region=="Caribbean","Caribbean","Other islands")
old_data <- read.csv("area-estimation/results/unfiltered-mixture/input-data.csv")
stopifnot(identical(d$island,old_data$island), isTRUE(all.equal(d$n_sp,old_data$n_sp)),
          isTRUE(all.equal(d$area_km2,old_data$area_km2)))
out <- "area-estimation/results/three-component-mixture"
plots.dir <- "area-estimation/plots/three-component-mixture"
dir.create(out,recursive=TRUE,showWarnings=FALSE)
dir.create(plots.dir,recursive=TRUE,showWarnings=FALSE)
write.csv(d,file.path(out,"input-data.csv"),row.names=FALSE)
driver <- FLXMRglmfix(family="gaussian",varFix=TRUE)

# Match the 2-component settings: shared variance, minimum mixing weight 0.1.
# Save all starts, including nonconvergence/component loss.
set.seed(261030)
starts <- data.frame(); best <- NULL
for(i in 1:200) {
  fit <- flexmix(y~x,data=d,k=3,model=driver,
    control=list(iter.max=2000,minprior=.1,tolerance=1e-8))
  initial_convergence <- fit@converged
  if (!fit@converged) fit <- flexmix(y~x,data=d,cluster=posterior(fit),model=driver,
    control=list(iter.max=20000,minprior=.1,tolerance=1e-8))
  starts <- rbind(starts,data.frame(start=i,logLik=as.numeric(logLik(fit)),
    components=fit@k,initial_convergence=initial_convergence,converged=fit@converged))
  if (fit@converged && fit@k==3 && (is.null(best) || logLik(fit)>logLik(best))) best <- fit
  if (i %% 20==0) cat("Full-data starts:",i,"/200\n")
}
write.csv(starts,file.path(out,"starts.csv"),row.names=FALSE)
stopifnot(!is.null(best),best@k==3)
saveRDS(best,file.path(out,"model.rds"))

# Leave one island out; constant mixing weights predict membership.
# Observed held-out richness is never used to select a component.
predictions <- data.frame()
for(i in seq_len(nrow(d))) {
  set.seed(261030+i)
  fit <- stepFlexmix(y~x,data=d[-i, ],k=3,nrep=30,model=driver,
    control=list(iter.max=2000,minprior=.1,tolerance=1e-8),verbose=FALSE)
  if (inherits(fit,"stepFlexmix")) fit <- getModel(fit,"BIC")
  initial_convergence <- fit@converged
  if (!fit@converged) fit <- flexmix(y~x,data=d[-i, ],cluster=posterior(fit),model=driver,
    control=list(iter.max=20000,minprior=.1,tolerance=1e-8))
  stopifnot(fit@converged)
  pars <- parameters(fit)
  mu <- pars[1, ]+pars[2, ]*d$x[i]
  weights <- prior(fit); sigma <- pars[3, ]
  terms <- log(weights)+dnorm(d$y[i],mu,sigma,log=TRUE)
  predictions <- rbind(predictions,data.frame(island=d$island[i],observed=d$n_sp[i],
    pred_log=sum(weights*mu),pred_mean=sum(weights*exp(mu+sigma^2/2)),
    log_density=max(terms)+log(sum(exp(terms-max(terms)))),components=fit@k,
    initial_convergence=initial_convergence,converged=fit@converged))
  write.csv(predictions,file.path(out,"held-out-predictions.csv"),row.names=FALSE)
  cat("Holdout:",i,"/61\n")
}
previous <- read.csv("area-estimation/results/unfiltered-mixture/model-comparison.csv")
previous$model[previous$model=="Latent mixture"] <- "Mixture: 2 components"
k <- attr(logLik(best),"df"); n <- nrow(d)
score <- data.frame(model="Mixture: 3 components",k=k,
  AICc=AIC(best)+2*k*(k+1)/(n-k-1),BIC=BIC(best),
  log_RMSE=sqrt(mean((log(predictions$observed)-predictions$pred_log)^2)),
  count_RMSE=sqrt(mean((predictions$observed-predictions$pred_mean)^2)),
  mean_negative_log_density=-mean(predictions$log_density),delta_AICc=NA_real_)
scores <- rbind(previous,score)
scores$delta_AICc <- scores$AICc-min(scores$AICc)
write.csv(scores,file.path(out,"model-comparison.csv"),row.names=FALSE)

# Component labels ordered by richness at 10,000 km2, not by region/type.
pars <- parameters(best); ordering <- order(pars[1, ])
pars <- pars[,ordering]; prob <- posterior(best)[,ordering]
labels <- c("Low","Middle","High")
d$component <- factor(labels[max.col(prob,ties.method="first")],levels=labels)
d$confidence <- pmax(prob[,1],prob[,2],prob[,3])
d$p_low <- prob[,1]; d$p_middle <- prob[,2]; d$p_high <- prob[,3]
write.csv(d,file.path(out,"island-membership.csv"),row.names=FALSE)
component_table <- data.frame(component=labels,intercept=pars[1, ],area_exponent=pars[2, ],
  sigma=pars[3, ],weight=prior(best)[ordering],hard_count=as.integer(table(d$component)))
write.csv(component_table,file.path(out,"component-parameters.csv"),row.names=FALSE)

# Label-invariant region/type agreement, for both mixture sizes.
adjusted_rand <- function(a,b) {
  t <- table(a,b)
  choose2 <- function(x) x*(x-1)/2
  observed <- sum(choose2(t))
  expected <- sum(choose2(rowSums(t)))*sum(choose2(colSums(t)))/choose2(sum(t))
  maximum <- (sum(choose2(rowSums(t)))+sum(choose2(colSums(t))))/2
  (observed-expected)/(maximum-expected)
}
old <- read.csv("area-estimation/results/unfiltered-mixture/island-membership.csv")
stopifnot(identical(old$island,d$island))
agreement <- data.frame()
for(variable in c("region","type","group")) {
  agreement <- rbind(agreement,data.frame(classification=variable,
    ARI_two=adjusted_rand(old$component,d[[variable]]),ARI_three=adjusted_rand(d$component,d[[variable]])))
}
write.csv(agreement,file.path(out,"group-agreement.csv"),row.names=FALSE)
heatmap <- data.frame(); soft_counts <- data.frame()
for(variable in c("region","type")) {
  tab <- as.data.frame(table(class=d[[variable]],component=d$component))
  tab$classification <- variable
  totals <- table(d[[variable]])
  tab$fraction <- tab$Freq/as.numeric(totals[as.character(tab$class)])
  heatmap <- rbind(heatmap,tab)
  write.csv(tab,file.path(out,paste0(variable,"-confusion.csv")),row.names=FALSE)
  for(group in unique(d[[variable]])) {
    for(j in 1:3) soft_counts <- rbind(soft_counts,data.frame(classification=variable,class=group,
      component=labels[j],expected_count=sum(prob[d[[variable]]==group,j])))
  }
}
write.csv(soft_counts,file.path(out,"soft-membership-counts.csv"),row.names=FALSE)

# Plots ----
theme_set(theme_classic(base_size=11))
cols <- c(Low="#0072B2",Middle="#E69F00",High="#CC79A7")
grid <- data.frame()
for(j in 1:3) {
  x <- seq(min(d$x),max(d$x),length.out=300)
  grid <- rbind(grid,data.frame(area_km2=10000*exp(x),richness=exp(pars[1,j]+pars[2,j]*x),component=labels[j]))
}
p1 <- ggplot(d,aes(area_km2,n_sp)) +
  geom_line(data=grid,aes(area_km2,richness,colour=component),linewidth=.8) +
  geom_point(aes(colour=component,shape=region),size=2.4) +
  geom_text_repel(aes(label=island,colour=component),seed=26,max.overlaps=Inf,max.time=10,size=2.9,show.legend=FALSE) +
  scale_x_log10() + scale_y_log10() + coord_cartesian(ylim=c(0.7,max(d$n_sp)*1.3)) + scale_colour_manual(values=cols) +
  labs(x=expression("Island area (km"^2*")"),y="Species richness",colour="Inferred component",shape="Region",
    caption="Component curves are fitted medians; names refer to richness at 10,000 km². Regions/types were not used to fit components.") +
  theme(legend.position="bottom",plot.caption=element_text(size=9))
p2 <- ggplot(heatmap,aes(component,class,fill=fraction)) + geom_tile(colour="white") +
  geom_text(aes(label=Freq),size=4) + facet_wrap(~classification,scales="free_y",ncol=1) +
  scale_fill_gradient(low="white",high="#0072B2",limits=c(0,1)) +
  labs(x="Inferred component",y=NULL,fill="Fraction within\nregion/type") + theme_minimal(base_size=12)
long <- data.frame()
for(j in 1:3) long <- rbind(long,data.frame(island=d$island,region=d$region,component=labels[j],probability=prob[,j]))
long$island <- factor(long$island,levels=d$island[order(d$region,d$component,d$confidence)])
long$component <- factor(long$component,levels=labels)
p3 <- ggplot(long,aes(probability,island,fill=component)) + geom_col(width=.8) +
  facet_grid(region~.,scales="free_y",space="free_y") + scale_fill_manual(values=cols) +
  labs(x="Component-membership probability",y=NULL,fill=NULL) + theme(legend.position="bottom")
metrics <- rbind(data.frame(model=scores$model,metric="Log-scale RMSE",value=scores$log_RMSE),
  data.frame(model=scores$model,metric="Negative log predictive density",value=scores$mean_negative_log_density))
p4 <- ggplot(metrics,aes(value,model)) + geom_point(size=3,colour="#0072B2") + facet_wrap(~metric,scales="free_x") +
  labs(x="Leave-one-out error (lower is better)",y=NULL)
figures <- list(mixture_labelled=p1,group_agreement=p2,membership_probabilities=p3,prediction_comparison=p4)
for(name in names(figures)) {
  width <- ifelse(name=="membership_probabilities",9,12)
  height <- ifelse(name=="membership_probabilities",16,ifelse(name=="prediction_comparison",5,8))
  ggsave(file.path(plots.dir,paste0(name,".pdf")),figures[[name]],width=width,height=height)
  ggsave(file.path(plots.dir,paste0(name,".png")),figures[[name]],width=width,height=height,dpi=180)
}

report <- c("# Three-component mixture on all 61 islands","",
"No islands excluded. Three latent power regressions use Gaussian log-richness errors and one common residual SD. Mixing weights are constant, with minprior=0.1, matching the previous two-component specification. Region and island type are not supplied to mixture fitting. Component names are ordered at area 10,000 km² and need not retain that ordering if curves cross.","",
"Full fit: 200 random starts, 2,000 iterations followed by continuation up to 20,000 if necessary, at unchanged 1e-8 tolerance. Best converged three-component fit selected. Cross-validation uses 30 starts per held-out island and the same continuation rule. Any component loss is recorded; failed fits are not silently discarded. Two-component and reference scores were reused only after checking identical source island names, areas and richness. Held-out prediction integrates component weights without using held-out richness to choose a component.","",
"| Model | AICc | BIC | Log CV RMSE | Negative log predictive density |","|---|---:|---:|---:|---:|")
for(i in seq_len(nrow(scores))) report <- c(report,sprintf("| %s | %.2f | %.2f | %.3f | %.3f |",scores$model[i],scores$AICc[i],scores$BIC[i],scores$log_RMSE[i],scores$mean_negative_log_density[i]))
report <- c(report,"","## Assignment agreement","","| Classification | ARI, two components | ARI, three components |","|---|---:|---:|")
for(i in seq_len(nrow(agreement))) report <- c(report,sprintf("| %s | %.3f | %.3f |",agreement$classification[i],agreement$ARI_two[i],agreement$ARI_three[i]))
report <- c(report,"",
sprintf("%d/61 islands have membership confidence below 0.8; %d below 0.6. Best full likelihood reproduced within 1e-4 in %d/200 starts. %d starts converged and retained three components. All %d holdout fits converged; %d retained three components.",sum(d$confidence<.8),sum(d$confidence<.6),sum(starts$converged & abs(starts$logLik-as.numeric(logLik(best)))<1e-4),sum(starts$converged & starts$components==3),nrow(predictions),sum(predictions$components==3)),"",
"Adjusted Rand index is label-invariant and corrects chance agreement (1: identical partitions, approximately 0: little agreement). Hard assignments are for visualization, not estimation; soft-count tables and individual probabilities are supplied. Associations are descriptive, in-sample, and do not establish a causal regional/geological mechanism. Shared variance/minimum weight restrict the tested mixture family. Model-selection exploration and AIC for weakly identified mixtures merit caution; predictive comparison is also reported. No Bayesian area estimates were changed.","",
"Run area-estimation/scripts/16-three-component-mixture.R from the project root. Outputs are in results/three-component-mixture and plots/three-component-mixture.")
writeLines(report,"area-estimation/reports/three-component-mixture.md")
print(scores); print(component_table); print(agreement)

# Interpretation of this saved run ----
interpretation <- c("", "## Interpretation", "",
"The three-component model improves AICc over two components by 5.87, but its leave-one-out log RMSE is 10.4% worse (1.066 versus 0.966) and predictive log density is worse. The pooled/known-region models remain competitive or better. This is evidence against adopting the three-component model for prediction from these data.", "",
"The fitted groups contain 5, 25 and 31 islands. The five-island component contains Aru Islands, Panay, Langkawi, Tioman and Réunion and has a negative area exponent (-0.735); the other exponents are 0.421 and 0.570. This small component accommodates contrasting observations rather than a coherent positive species-area relationship. Component labels are ordered only at 10,000 km2.", "",
"Group correspondence is weak for both three regions (ARI -0.010) and island type (ARI 0.079). Caribbean islands split 9/7 between Middle/High; Australasian islands occur in all three groups. Most full-data starts collapsed to fewer components, and only 37/61 cross-validation fits retained three components. Thus predictive scores describe the specified fitting procedure, including component loss, rather than a stable three-component partition. Further constrained-mixture variants would constitute new analyses.")
write(interpretation,"area-estimation/reports/three-component-mixture.md",append=TRUE)
