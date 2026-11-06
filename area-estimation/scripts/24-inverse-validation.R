# Pilot inverse validation: hold out entire geographic blocks.
# Approximate parameter uncertainty; not a new final Bayesian area analysis.
library(MASS)
library(lme4)
library(splines)
library(ggplot2)
a <- readRDS('area-estimation/results/model-suite/analysis-input.rds')
d <- a$data
ids <- c('lognormal:area','lognormal:isolation','lognormal:previous_best',
         'lognormal:curve','NB:previous_best','mixed:block')
out <- 'area-estimation/results/model-suite'
results <- data.frame()
z <- seq(1,6,by=.01)
area <- 10^z
prior <- rep(1,length(z)); prior[c(1,length(z))] <- .5
prior <- prior/sum(prior)
set.seed(271026)
ndraw <- 800
for(id in ids) {
 candidate <- a$candidates[a$candidates$id==id, ]
 family <- candidate$family
 form <- a$forms[[candidate$model]]
 for(block in levels(d$block)) {
  train <- d[d$block!=block, ]; test <- d[d$block==block, ]
  if(family=='NB') fit <- glm.nb(update(form,n_sp~.),data=train,control=glm.control(maxit=200))
  if(family=='lognormal') fit <- lm(update(form,log_richness~.),data=train)
  if(family=='mixed') fit <- lmer(update(form,log_richness~.),data=train,REML=FALSE)
  if(family=='lognormal') {
   # Exact normal-regression reference-posterior draws under p(beta,sigma2) proportional to 1/sigma2.
   # These differ from the proper priors used for the final Sunken Islands estimate.
   sig <- sqrt(sum(residuals(fit)^2)/rchisq(ndraw,df.residual(fit)))
   C <- vcov(fit)/sigma(fit)^2
   beta <- mvrnorm(ndraw,rep(0,length(coef(fit))),C)*sig
   beta <- sweep(beta,2,coef(fit),'+')
  }
  if(family=='NB') {
   # Wald approximation; beta and log-dispersion are approximated independently.
   beta <- mvrnorm(ndraw,coef(fit),vcov(fit))
   dispersion <- exp(rnorm(ndraw,log(fit$theta),fit$SE.theta/fit$theta))
  }
  if(family=='mixed') {
   beta <- mvrnorm(ndraw,fixef(fit),as.matrix(vcov(fit)))
   # Predict a new geographic group: integrate its normal random intercept.
   # Variance components held at ML estimates in this screening approximation.
   sig <- rep(sqrt(sigma(fit)^2+as.numeric(VarCorr(fit)$block[1])),ndraw)
  }
  for(i in seq_len(nrow(test))) {
   grid <- test[rep(i,length(z)), ]
   grid$log_area <- log(area/10000)
   fixed_terms <- delete.response(terms(fit))
   X <- model.matrix(fixed_terms,grid,contrasts.arg=if(family=='mixed') attr(getME(fit,'X'),'contrasts') else fit$contrasts)
   if(family=='mixed') X <- X[,names(fixef(fit)),drop=FALSE]
   eta <- X %*% t(beta)
   if(family=='NB') {
    L <- dnbinom(test$n_sp[i],mu=exp(eta),size=rep(dispersion,each=length(z)))
   } else {
    # Rounding continuous lognormal predictions to the observed integer count.
    L <- plnorm(test$n_sp[i]+.5,eta,rep(sig,each=length(z)))-
         plnorm(max(0,test$n_sp[i]-.5),eta,rep(sig,each=length(z)))
   }
   L <- matrix(L,nrow=length(z),ncol=ndraw)
   mass <- rowMeans(L)*prior
   mass <- mass/sum(mass)
   qs <- 10^approx(cumsum(mass),z,xout=c(.025,.5,.975),ties='ordered')$y
   stopifnot(all(is.finite(qs)),abs(sum(mass)-1)<1e-8)
   results <- rbind(results,data.frame(id,block,island=test$island[i],truth=test$area_km2[i],
     lower=qs[1],median=qs[2],upper=qs[3],
     covered=test$area_km2[i]>=qs[1] & test$area_km2[i]<=qs[3],
     absolute_log10_error=abs(log10(qs[2]/test$area_km2[i])),
     interval_log10_width=log10(qs[3]/qs[1]),
     mass_near_lower=sum(mass[z<1.1]),mass_near_upper=sum(mass[z>5.9])))
  }
 }
 cat(id,'inverse complete\n')
}
# Prior-only benchmark shows how much recovery exceeds an uninformative guess.
qs <- 10^approx(cumsum(prior),z,xout=c(.025,.5,.975),ties='ordered')$y
results <- rbind(results,data.frame(id='prior_only',block=d$block,island=d$island,truth=d$area_km2,
 lower=qs[1],median=qs[2],upper=qs[3],covered=d$area_km2>=qs[1]&d$area_km2<=qs[3],
 absolute_log10_error=abs(log10(qs[2]/d$area_km2)),interval_log10_width=log10(qs[3]/qs[1]),
 mass_near_lower=sum(prior[z<1.1]),mass_near_upper=sum(prior[z>5.9])))
summary <- data.frame()
for(id in unique(results$id)) {
 r <- results[results$id==id, ]
 summary <- rbind(summary,data.frame(id,n=nrow(r),coverage=mean(r$covered),
   median_absolute_log10_error=median(r$absolute_log10_error),
   median_error_factor=10^median(r$absolute_log10_error),
   median_interval_factor=10^median(r$interval_log10_width)))
}
write.csv(results,file.path(out,'inverse-held-out-islands.csv'),row.names=FALSE)
write.csv(summary,file.path(out,'inverse-summary.csv'),row.names=FALSE)
p <- ggplot(results[results$id!='prior_only', ],aes(truth,median))+
 geom_abline(slope=1,intercept=0,colour='grey60',linetype=2)+
 geom_linerange(aes(ymin=lower,ymax=upper),alpha=.2)+geom_point(size=1.4,colour='#0072B2')+
 scale_x_log10()+scale_y_log10()+facet_wrap(~id,ncol=3)+theme_classic(base_size=10)+
 labs(x=expression('Known area (km'^2*')'),y=expression('Inferred median and 95% interval (km'^2*')'))
ggsave('area-estimation/plots/model-suite/inverse-validation.png',p,width=11,height=7,dpi=180)
ggsave('area-estimation/plots/model-suite/inverse-validation.pdf',p,width=11,height=7)
print(summary)
