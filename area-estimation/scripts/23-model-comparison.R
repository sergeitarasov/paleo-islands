# Run from paleo-islands. All 61 islands; no richness-based filtering.
library(readxl)
library(MASS)
library(lme4)
library(splines)
library(ggplot2)

# Inputs and candidates ----
d <- as.data.frame(read_excel('data/richness/island-richness-area.xlsx',sheet='Island isolation'))
blocks <- read.csv('area-estimation/results/model-suite/geographic-blocks.csv')
d$block <- factor(blocks$block[match(d$island, blocks$island)])
d$region <- factor(d$region)
d$type <- factor(d$type)
d$group <- factor(ifelse(d$region=='Caribbean','Caribbean','Other'))
d$log_area <- log(d$area_km2/10000)
d$log_mainland <- log1p(d$distance_mainland_km)
d$log_richness <- log(d$n_sp)
stopifnot(nrow(d)==61, !anyNA(d$block), all(d$n_sp>0))
out <- 'area-estimation/results/model-suite'
plots.dir <- 'area-estimation/plots/model-suite'
dir.create(plots.dir,recursive=TRUE,showWarnings=FALSE)
write.csv(d,file.path(out,'input-data.csv'),row.names=FALSE)
write.csv(as.data.frame(table(d$region,d$type)),file.path(out,'region-type-counts.csv'),row.names=FALSE)
forms <- list(area=~log_area, isolation=~log_area+log_mainland,
 region_type=~log_area+region+type,
 region_type_isolation=~log_area+region+type+log_mainland,
 previous_best=~log_area*group+log_mainland,
 curve=~ns(log_area,df=3),
 curve_adjusted=~ns(log_area,df=3)+region+type+log_mainland)
candidates <- data.frame()
for(family in c('lognormal','NB')) {
 for(model in names(forms)) candidates <- rbind(candidates,data.frame(id=paste(family,model,sep=':'),family,model))
}
forms$block <- ~log_area+type+log_mainland+(1|block)
forms$block_region <- ~log_area+region+type+log_mainland+(1|block)
candidates <- rbind(candidates,data.frame(id=c('mixed:block','mixed:block_region'),family='mixed',model=c('block','block_region')))

# Repeated fitting/prediction helpers keep validation identical across models.
fit_model <- function(family, form, train) {
 if(family=='NB') return(glm.nb(update(form,n_sp~.),data=train,control=glm.control(maxit=200)))
 if(family=='mixed') return(lmer(update(form,log_richness~.),data=train,REML=FALSE))
 lm(update(form,log_richness~.),data=train)
}
forecast <- function(fit,family,test,new_block=FALSE) {
 if(family=='NB') {
  mu <- as.numeric(predict(fit,test,type='response'))
  med <- qnbinom(.5,mu=mu,size=fit$theta)
  return(data.frame(pred_mean=mu,pred_log=log(pmax(med,.5)),log_probability=dnbinom(test$n_sp,mu=mu,size=fit$theta,log=TRUE)))
 }
 if(family=='mixed') {
  eta <- as.numeric(predict(fit,test,allow.new.levels=TRUE,re.form=if(new_block) NA else NULL))
  variance <- sigma(fit)^2
  if(new_block) variance <- variance+as.numeric(VarCorr(fit)$block[1])
 } else {
  eta <- as.numeric(predict(fit,test))
  variance <- sum(residuals(fit)^2)/nobs(fit)
 }
 prob <- plnorm(test$n_sp+.5,eta,sqrt(variance))-plnorm(pmax(0,test$n_sp-.5),eta,sqrt(variance))
 data.frame(pred_mean=exp(eta+variance/2),pred_log=eta,log_probability=log(pmax(prob,1e-300)))
}
# Lognormal likelihoods are evaluated on integer bins for predictive scoring;
# AICc remains comparable only within a response family.
models <- list(); fit_stats <- data.frame(); predictions <- data.frame(); issues <- data.frame()
for(j in seq_len(nrow(candidates))) {
 id <- candidates$id[j]; family <- candidates$family[j]; form <- forms[[candidates$model[j]]]
 messages <- character()
 fit <- withCallingHandlers(fit_model(family,form,d),warning=function(w){messages <<- c(messages,conditionMessage(w));invokeRestart('muffleWarning')})
 models[[id]] <- fit
 k <- attr(logLik(fit),'df')
 singular <- family=='mixed' && isSingular(fit)
 fit_stats <- rbind(fit_stats,data.frame(id,family,k,AICc=AIC(fit)+2*k*(k+1)/(nrow(d)-k-1),singular))
 if(length(messages)) issues <- rbind(issues,data.frame(id,validation='full',fold='all',message=paste(unique(messages),collapse=' | ')))
 for(validation in c('island','block')) {
  folds <- if(validation=='island') d$island else levels(d$block)
  for(fold in folds) {
   hold <- if(validation=='island') d$island==fold else d$block==fold
   train <- d[!hold, ]; test <- d[hold, ]
   messages <- character()
   fitted <- tryCatch(withCallingHandlers(fit_model(family,form,train),warning=function(w){messages <<- c(messages,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e)e)
   if(inherits(fitted,'error')) {
    issues <- rbind(issues,data.frame(id,validation,fold,message=conditionMessage(fitted))); next
   }
   # Singleton blocks are unseen even in leave-one-island-out validation.
   unseen <- validation=='block' || !as.character(test$block[1]) %in% as.character(train$block)
   p <- forecast(fitted,family,test,new_block=unseen)
   predictions <- rbind(predictions,data.frame(id,validation,fold,island=test$island,
      observed=test$n_sp,observed_log=test$log_richness,p))
   if(family=='mixed' && isSingular(fitted)) messages <- c(messages,'singular random effect')
   if(family=='NB' && !fitted$converged) messages <- c(messages,'GLM not converged')
   if(length(messages)) issues <- rbind(issues,data.frame(id,validation,fold,message=paste(unique(messages),collapse=' | ')))
  }
 }
 cat(id,'complete\n')
}
comparison <- data.frame()
for(id in candidates$id) {
 for(validation in c('island','block')) {
  z <- predictions[predictions$id==id & predictions$validation==validation, ]
  comparison <- rbind(comparison,data.frame(id,validation,n=nrow(z),
    negative_log_score=-mean(z$log_probability),log_RMSE=sqrt(mean((z$observed_log-z$pred_log)^2)),
    count_RMSE=sqrt(mean((z$observed-z$pred_mean)^2))))
 }
}
comparison <- merge(comparison,fit_stats,by='id')
comparison <- comparison[order(comparison$validation,comparison$negative_log_score), ]
write.csv(comparison,file.path(out,'model-comparison.csv'),row.names=FALSE)
write.csv(predictions,file.path(out,'held-out-predictions.csv'),row.names=FALSE)
write.csv(issues,file.path(out,'fit-issues.csv'),row.names=FALSE)
saveRDS(models,file.path(out,'models.rds'))
saveRDS(list(data=d,candidates=candidates,forms=forms),file.path(out,'analysis-input.rds'))

# Paired block-score comparisons with previous model ----
paired <- data.frame()
ref <- predictions[predictions$id=='lognormal:previous_best' & predictions$validation=='block', ]
for(id in candidates$id) {
 z <- predictions[predictions$id==id & predictions$validation=='block', ]
 difference <- z$log_probability-ref$log_probability[match(z$island,ref$island)]
 totals <- aggregate(difference,list(block=z$fold),sum)$x
 paired <- rbind(paired,data.frame(id,log_score_gain=sum(totals),SE=sqrt(length(totals)*var(totals))))
}
write.csv(paired,file.path(out,'paired-block-comparison.csv'),row.names=FALSE)

# Figures ----
theme_set(theme_classic(base_size=11))
p <- ggplot(comparison,aes(negative_log_score,reorder(id,negative_log_score),colour=family)) +
 geom_point(size=2.5)+facet_wrap(~validation)+labs(x='Held-out negative log probability (lower is better)',y=NULL,colour='Model family')+
 theme(legend.position='bottom')
ggsave(file.path(plots.dir,'predictive-comparison.png'),p,width=11,height=7,dpi=180)
ggsave(file.path(plots.dir,'predictive-comparison.pdf'),p,width=11,height=7)
capture.output(sessionInfo(),file=file.path(out,'session-info.txt'))
print(comparison)
