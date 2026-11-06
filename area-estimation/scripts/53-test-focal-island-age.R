# Updated age tab: compare models on identical islands within each analysis.
library(readxl)
library(ggplot2)
library(ggrepel)
d <- as.data.frame(read_excel('data/richness/island-richness-area.xlsx',sheet='Focal islands ages'))
stopifnot(nrow(d)==21,all(d$area_km2>0),all(d$age_Ma>0),all(d$n_sp>0))
d$x <- log(d$area_km2)
d$y <- log(d$n_sp)
d$log_age <- log(d$age_Ma)
d$group <- factor(ifelse(d$region=='Caribbean','Caribbean','Other'),levels=c('Caribbean','Other'))
d$low_confidence <- grepl('Low-confidence',d$age_source,fixed=TRUE)
out <- 'area-estimation/results/focal-island-age'
plots <- 'area-estimation/plots/focal-island-age'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
dir.create(plots,recursive=TRUE,showWarnings=FALSE)
write.csv(d,file.path(out,'input.csv'),row.names=FALSE)
writeLines(unname(tools::md5sum('data/richness/island-richness-area.xlsx')),file.path(out,'workbook-md5.txt'))
forms <- list('Area only'=y~x,'Area + log age'=y~x+log_age,'Area + age'=y~x+age_Ma,
 'Regional slopes'=y~x*group,'Regional slopes + log age'=y~x*group+log_age,'Regional slopes + age'=y~x*group+age_Ma)
samples <- list('All 21'=d,'Without five low-confidence ages'=d[!d$low_confidence, ],
 'Without separation-age proxies'=d[!d$island %in% c('Madagascar','New Zealand'), ],
 'Without low-confidence and regional 35 Ma proxies'=d[!d$low_confidence & d$age_Ma!=35, ])
summary <- data.frame()
coefficients <- data.frame()
predictions <- data.frame()
fits <- list()
for(subset in names(samples)) {
 data <- samples[[subset]]
 for(model in names(forms)) {
  fit <- lm(forms[[model]],data=data)
  if(subset=='All 21') fits[[model]] <- fit
  k <- attr(logLik(fit),'df')
  n <- nrow(data)
  loo <- numeric(n)
  lp <- lower <- upper <- numeric(n)
  for(i in 1:n) {
   train <- lm(forms[[model]],data=data[-i, ])
   pred <- predict(train,newdata=data[i, ],se.fit=TRUE)
   loo[i] <- pred$fit
   scale <- sqrt(sigma(train)^2+pred$se.fit^2)
   lp[i] <- dt((data$y[i]-loo[i])/scale,df.residual(train),log=TRUE)-log(scale)
   lower[i] <- loo[i]+qt(.025,df.residual(train))*scale
   upper[i] <- loo[i]+qt(.975,df.residual(train))*scale
  }
  summary <- rbind(summary,data.frame(subset,model,n,k,R2=summary(fit)$r.squared,adjusted_R2=summary(fit)$adj.r.squared,
   AICc=AIC(fit)+2*k*(k+1)/(n-k-1),LOO_log_RMSE=sqrt(mean((data$y-loo)^2)),
   LOO_negative_log_density=-mean(lp),LOO_coverage=mean(data$y>=lower & data$y<=upper)))
  cf <- coef(summary(fit)); ci <- confint(fit)
  coefficients <- rbind(coefficients,data.frame(subset,model,term=rownames(cf),estimate=cf[,1],se=cf[,2],p=cf[,4],lower=ci[,1],upper=ci[,2]))
  predictions <- rbind(predictions,data.frame(subset,model,island=data$island,group=data$group,observed=data$n_sp,log_observed=data$y,log_prediction=loo,lower=exp(lower),upper=exp(upper)))
 }
}
write.csv(summary,file.path(out,'model-comparison.csv'),row.names=FALSE)
write.csv(coefficients,file.path(out,'coefficients.csv'),row.names=FALSE)
write.csv(predictions,file.path(out,'held-out-predictions.csv'),row.names=FALSE)
saveRDS(fits,file.path(out,'full-data-fits.rds'))
# Nested tests: age conditional on area, and on area plus regional slopes.
tests <- data.frame()
for(pair in list(c('Area only','Area + log age'),c('Area only','Area + age'),c('Regional slopes','Regional slopes + log age'),c('Regional slopes','Regional slopes + age'))) {
 a <- anova(fits[[pair[1]]],fits[[pair[2]]])
 tests <- rbind(tests,data.frame(baseline=pair[1],expanded=pair[2],F=a$F[2],p=a$`Pr(>F)`[2]))
}
write.csv(tests,file.path(out,'nested-tests.csv'),row.names=FALSE)
# Leave-one-out age coefficient shows whether one island determines the effect.
influence <- data.frame()
for(i in 1:nrow(d)) {
 fit <- lm(y~x+log_age,data=d[-i, ])
 ci <- confint(fit)['log_age', ]
 influence <- rbind(influence,data.frame(omitted=d$island[i],age_effect=coef(fit)['log_age'],lower=ci[1],upper=ci[2]))
}
write.csv(influence,file.path(out,'leave-one-out-age-effect.csv'),row.names=FALSE)

# Plot held-out performance and age effect ----
p <- ggplot(summary[summary$subset=='All 21', ],aes(LOO_log_RMSE,reorder(model,LOO_log_RMSE)))+
 geom_point(size=3,colour='#0072B2')+theme_classic(base_size=12)+labs(x='Leave-one-out log-richness RMSE (lower is better)',y=NULL)
ggsave(file.path(plots,'prediction-comparison.png'),p,width=9,height=4,dpi=180)
ggsave(file.path(plots,'prediction-comparison.pdf'),p,width=9,height=4)
fit <- fits[['Area + log age']]
# Partial residuals retain the log-age contribution after adjusting for area.
d$partial <- residuals(fit)+coef(fit)['log_age']*d$log_age
age_grid <- seq(min(d$log_age),max(d$log_age),length.out=200)
ci <- confint(fit)['log_age', ]
line <- data.frame(log_age=age_grid,estimate=coef(fit)['log_age']*age_grid,
 lower=pmin(ci[1]*age_grid,ci[2]*age_grid),upper=pmax(ci[1]*age_grid,ci[2]*age_grid))
p <- ggplot(d,aes(log_age,partial))+
 geom_ribbon(data=line,inherit.aes=FALSE,mapping=aes(x=log_age,ymin=lower,ymax=upper),fill='#0072B2',alpha=.12)+
 geom_line(data=line,aes(log_age,estimate),inherit.aes=FALSE,colour='#0072B2')+
 geom_point(aes(colour=group,shape=low_confidence),size=2.5)+geom_text_repel(aes(label=island),size=3,max.overlaps=30,seed=2)+
 scale_x_continuous(breaks=log(c(.1,1,10,100)),labels=c('0.1','1','10','100'))+
 scale_colour_manual(values=c(Caribbean='#0072B2',Other='#D55E00'))+
 theme_classic(base_size=11)+theme(legend.position='bottom')+
 labs(x='Assigned age (Ma; logarithmic axis)',y='Log-richness partial residual, adjusted for area',colour=NULL,shape='Low-confidence age',
 caption='Shading: 95% coefficient interval for log-age contribution. Ages are heterogeneous geological proxies.')
ggsave(file.path(plots,'age-effect.png'),p,width=10,height=6,dpi=180)
ggsave(file.path(plots,'age-effect.pdf'),p,width=10,height=6)
capture.output(sessionInfo(),file=file.path(out,'session-info.txt'))
print(summary)
print(coefficients[coefficients$term %in% c('log_age','age_Ma'), ])
print(tests)
