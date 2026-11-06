# Run from paleo-islands. Identical 21 test islands for all four candidates.
library(ggplot2)
d21 <- read.csv('area-estimation/results/selected-21/input.csv')
d27 <- read.csv('area-estimation/results/selected-27/input.csv')
stopifnot(nrow(d21)==21,nrow(d27)==27,all(d21$island %in% d27$island))
j <- match(d21$island,d27$island)
stopifnot(all(d21$n_sp==d27$n_sp[j]),all(d21$area_km2==d27$area_km2[j]))
samples <- list('21 islands'=d21,'27 islands'=d27)
forms <- list('Pooled slope'=y~x,'Separate slopes'=y~x*group)
out <- 'area-estimation/results/matched-21-27'
plots.dir <- 'area-estimation/plots/matched-21-27'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
dir.create(plots.dir,recursive=TRUE,showWarnings=FALSE)
write.csv(d21,file.path(out,'held-out-islands.csv'),row.names=FALSE)
write.csv(d27[!d27$island %in% d21$island, ],file.path(out,'six-additional-islands.csv'),row.names=FALSE)

# Exact Student-t predictive distribution for Gaussian linear regression under
# reference prior p(beta,sigma2) proportional to 1/sigma2. Includes coefficient
# and residual-variance uncertainty; no MCMC or asymptotic parameter draws.
# This is a validation calculation, not a replacement of the proper-prior
# Bayesian Sunken Islands estimates.
z <- seq(1,6,by=.0025)
area <- 10^z
results <- data.frame()
curves <- data.frame()
for(training in names(samples)) {
 data <- samples[[training]]
 data$x <- log(data$area_km2/10000)
 data$y <- log(data$n_sp)
 data$group <- factor(ifelse(data$region=='Caribbean','Caribbean','Other'),levels=c('Caribbean','Other'))
 for(model in names(forms)) {
  for(i in seq_len(nrow(d21))) {
   test <- data[data$island==d21$island[i], ]
   train <- data[data$island!=d21$island[i], ]
   fit <- lm(forms[[model]],data=train)
   # Richness prediction with known area, on the common held-out island.
   pred <- predict(fit,test,se.fit=TRUE)
   mu <- as.numeric(pred$fit)
   scale <- sqrt(sigma(fit)^2+as.numeric(pred$se.fit)^2)
   df <- df.residual(fit)
   lower <- mu+qt(.025,df)*scale
   upper <- mu+qt(.975,df)*scale
   probability <- pt((log(test$n_sp+.5)-mu)/scale,df)-pt((log(test$n_sp-.5)-mu)/scale,df)
   # Hide area: evaluate the same posterior predictive at all possible areas.
   grid <- test[rep(1,length(area)), ]
   grid$x <- log(area/10000)
   gp <- predict(fit,grid,se.fit=TRUE)
   gs <- sqrt(sigma(fit)^2+as.numeric(gp$se.fit)^2)
   loglik <- dt((test$y-as.numeric(gp$fit))/gs,df,log=TRUE)-log(gs)
   # Uniform prior on log10(area) in [1,6]. Normalize by trapezoidal integration.
   density <- exp(loglik-max(loglik))
   increments <- diff(z)*(head(density,-1)+tail(density,-1))/2
   density <- density/sum(increments)
   cdf <- c(0,cumsum(increments)/sum(increments))
   q <- approx(cdf,z,xout=c(.025,.5,.975),ties='ordered')$y
   truth <- log10(test$area_km2)
   area_log_score <- log(approx(z,density,xout=truth)$y)
   interval_score <- q[3]-q[1]+40*(q[1]-truth)*(truth<q[1])+40*(truth-q[3])*(truth>q[3])
   results <- rbind(results,data.frame(training,model,island=test$island,group=as.character(test$group),
    n_training=nrow(train),n_sp=test$n_sp,area_km2=test$area_km2,
    log_observed=test$y,log_prediction=mu,richness_lower=exp(lower),richness_upper=exp(upper),
    log_count_probability=log(probability),richness_covered=test$y>=lower & test$y<=upper,
    area_median=10^q[2],area_lower=10^q[1],area_upper=10^q[3],
    area_covered=truth>=q[1] & truth<=q[3],area_absolute_log10_error=abs(q[2]-truth),
    area_interval_log10_width=q[3]-q[1],area_interval_score=interval_score,
    area_log_density=area_log_score))
   if(test$group=='Other') curves <- rbind(curves,data.frame(training,model,island=test$island,z,density))
  }
 }
}
stopifnot(nrow(results)==84,all(is.finite(results$area_log_density)),all(is.finite(results$log_count_probability)))
summary <- data.frame()
for(training in names(samples)) for(model in names(forms)) for(group in c('All','Caribbean','Other')) {
 r <- results[results$training==training & results$model==model, ]
 if(group!='All') r <- r[r$group==group, ]
 summary <- rbind(summary,data.frame(training,model,group,n=nrow(r),
  richness_log_RMSE=sqrt(mean((r$log_observed-r$log_prediction)^2)),
  richness_negative_log_score=-mean(r$log_count_probability),richness_coverage=mean(r$richness_covered),
  area_median_error_factor=10^median(r$area_absolute_log10_error),
  area_mean_absolute_log10_error=mean(r$area_absolute_log10_error),
  area_coverage=mean(r$area_covered),area_median_interval_factor=10^median(r$area_interval_log10_width),
  area_mean_interval_score=mean(r$area_interval_score),area_negative_log_score=-mean(r$area_log_density)))
}
write.csv(results,file.path(out,'held-out-predictions.csv'),row.names=FALSE)
write.csv(summary,file.path(out,'summary.csv'),row.names=FALSE)
write.csv(curves,file.path(out,'non-caribbean-area-densities.csv'),row.names=FALSE)

# Paired descriptive bootstrap: common islands, stratified to retain 16/5 mix.
# Shared training data and archipelago dependence limit interpretation of CIs.
reference <- results[results$training=='21 islands' & results$model=='Pooled slope', ]
set.seed(271027)
boot <- matrix(0,5000,21)
for(b in 1:5000) boot[b, ] <- c(sample(which(reference$group=='Caribbean'),16,replace=TRUE),sample(which(reference$group=='Other'),5,replace=TRUE))
paired <- data.frame()
for(training in names(samples)) for(model in names(forms)) {
 r <- results[results$training==training & results$model==model, ]
 r <- r[match(reference$island,r$island), ]
 differences <- data.frame(richness_squared_log_error=(r$log_observed-r$log_prediction)^2-(reference$log_observed-reference$log_prediction)^2,
  richness_negative_log_score=-r$log_count_probability+reference$log_count_probability,
  area_absolute_log10_error=r$area_absolute_log10_error-reference$area_absolute_log10_error,
  area_interval_score=r$area_interval_score-reference$area_interval_score,
  area_negative_log_score=-r$area_log_density+reference$area_log_density)
 for(metric in names(differences)) {
  delta <- differences[[metric]]
  means <- rowMeans(matrix(delta[boot],5000,21))
  paired <- rbind(paired,data.frame(training,model,metric,mean_change=mean(delta),
   lower=quantile(means,.025),upper=quantile(means,.975)))
 }
}
write.csv(paired,file.path(out,'paired-bootstrap.csv'),row.names=FALSE)

# Plots ----
results$panel <- paste(results$training,results$model,sep=': ')
colours <- c(Caribbean='#0072B2',Other='#D55E00')
p <- ggplot(results,aes(area_km2,area_median,colour=group))+
 geom_abline(slope=1,intercept=0,linetype=2,colour='grey60')+
 geom_linerange(aes(ymin=area_lower,ymax=area_upper),alpha=.3)+geom_point(size=2)+
 scale_x_log10()+scale_y_log10()+facet_wrap(~panel,ncol=2)+scale_colour_manual(values=colours)+
 theme_classic(base_size=11)+theme(legend.position='bottom')+
 labs(x=expression('Known area (km'^2*')'),y=expression('Inferred median and 95% interval (km'^2*')'),colour=NULL,
 caption='Same 21 held-out islands in every panel; area prior: log-uniform, 10–1,000,000 km².')
ggsave(file.path(plots.dir,'matched-area-recovery.png'),p,width=11,height=9,dpi=180)
ggsave(file.path(plots.dir,'matched-area-recovery.pdf'),p,width=11,height=9)
p <- ggplot(results,aes(n_sp,exp(log_prediction),colour=group))+
 geom_abline(slope=1,intercept=0,linetype=2,colour='grey60')+
 geom_linerange(aes(ymin=richness_lower,ymax=richness_upper),alpha=.3)+geom_point(size=2)+
 scale_x_log10()+scale_y_log10()+facet_wrap(~panel,ncol=2)+scale_colour_manual(values=colours)+
 theme_classic(base_size=11)+theme(legend.position='bottom')+
 labs(x='Observed richness',y='Predicted median and 95% interval',colour=NULL)
ggsave(file.path(plots.dir,'matched-richness-prediction.png'),p,width=11,height=9,dpi=180)
ggsave(file.path(plots.dir,'matched-richness-prediction.pdf'),p,width=11,height=9)
print(summary[summary$group=='All', ])
print(summary[summary$group=='Other', ])
print(paired)

# Direct paired contrasts (negative change means improved prediction) ----
contrasts <- data.frame(a_training=c('27 islands','27 islands','21 islands','27 islands'),
 a_model=c('Pooled slope','Separate slopes','Separate slopes','Separate slopes'),
 b_training=c('21 islands','21 islands','21 islands','27 islands'),
 b_model=c('Pooled slope','Separate slopes','Pooled slope','Pooled slope'))
direct <- data.frame()
for(i in 1:nrow(contrasts)) {
 a <- results[results$training==contrasts$a_training[i] & results$model==contrasts$a_model[i], ]
 b <- results[results$training==contrasts$b_training[i] & results$model==contrasts$b_model[i], ]
 b <- b[match(a$island,b$island), ]
 differences <- data.frame(richness_squared_log_error=(a$log_observed-a$log_prediction)^2-(b$log_observed-b$log_prediction)^2,
 area_absolute_log10_error=a$area_absolute_log10_error-b$area_absolute_log10_error,
 area_interval_score=a$area_interval_score-b$area_interval_score,
 area_negative_log_score=-a$area_log_density+b$area_log_density)
 for(metric in names(differences)) {
  delta <- differences[[metric]]
  means <- rowMeans(matrix(delta[boot],5000,21))
  direct <- rbind(direct,data.frame(contrasts[i, ],metric,change=mean(delta),lower=unname(quantile(means,.025)),upper=unname(quantile(means,.975))))
 }
}
write.csv(direct,file.path(out,'direct-paired-contrasts.csv'),row.names=FALSE)

# Full-data slope estimates: descriptive, distinct from held-out validation.
slopes <- data.frame()
for(training in names(samples)) {
 data <- samples[[training]]
 data$x <- log(data$area_km2/10000)
 data$y <- log(data$n_sp)
 data$group <- factor(ifelse(data$region=='Caribbean','Caribbean','Other'),levels=c('Caribbean','Other'))
 fit <- lm(y~x*group,data=data)
 for(group in c('Caribbean','Other')) {
  contrast <- c(0,1,0,as.numeric(group=='Other'))
  estimate <- sum(contrast*coef(fit))
  se <- sqrt(as.numeric(t(contrast)%*%vcov(fit)%*%contrast))
  slopes <- rbind(slopes,data.frame(training,group,slope=estimate,
   lower=estimate-qt(.975,df.residual(fit))*se,upper=estimate+qt(.975,df.residual(fit))*se,
   interaction_p=coef(summary(fit))['x:groupOther','Pr(>|t|)']))
 }
}
write.csv(slopes,file.path(out,'full-data-slopes.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,'session-info.txt'))
