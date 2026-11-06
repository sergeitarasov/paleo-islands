# Posterior predictive richness at the two fixed geological areas.
library(ggplot2)
fits <- list('27 islands'=readRDS('area-estimation/results/bayesian-selected27-slopes-10ma/regression-posteriors.rds'),
 '22 islands'=readRDS('area-estimation/results/bayesian-no-malay-22-10ma/regression-posteriors.rds'))
s <- read.csv('area-estimation/results/richness-geological-areas/predicted-richness.csv')
plots <- 'area-estimation/plots/richness-geological-areas'
dir.create(plots,recursive=TRUE,showWarnings=FALSE)
curves <- data.frame()
x <- exp(seq(log(.05),log(3000),length.out=1200))
for(dataset in names(fits)) for(model in c('simple_prior1','separate_prior1')) {
 fit <- fits[[dataset]][[model]]
 draws <- fit$draws
 p <- ncol(draws)-1
 for(scenario in c('Scotese','Straume')) {
  area <- if(scenario=='Scotese') 32000 else 114000
  new <- data.frame(x=log(area/10000),group=factor('Other',levels=c('Caribbean','Other')))
  X <- model.matrix(fit$formula,new)
  mu <- as.vector(draws[,1:p] %*% as.vector(X))
  sigma <- exp(draws[,p+1])
  density <- numeric(length(x))
  for(j in seq_along(x)) density[j] <- mean(dlnorm(x[j],mu,sigma))
  curves <- rbind(curves,data.frame(dataset,model,scenario,richness=x,density))
 }
}
write.csv(curves,'area-estimation/results/richness-geological-areas/predictive-densities.csv',row.names=FALSE)
curves$dataset <- factor(curves$dataset,levels=c('27 islands','22 islands'))
s$dataset <- factor(s$dataset,levels=c('27 islands','22 islands'))
curves$model <- factor(curves$model,levels=c('simple_prior1','separate_prior1'),labels=c('Pooled','Separate slopes: non-Caribbean'))
s$model <- factor(s$model,levels=c('simple_prior1','separate_prior1'),labels=levels(curves$model))
colours <- c(Scotese='#0072B2',Straume='#D55E00')
s$annotation <- sprintf('%s: mean %.1f; median %.1f\n50%%: %.1f–%.1f; 95%%: %.1f–%.1f',s$scenario,s$predictive_mean,s$predictive_median,s$lower50,s$upper50,s$lower95,s$upper95)
s$y <- ifelse(s$scenario=='Scotese',.059,.048)
p <- ggplot(curves,aes(richness,density,colour=scenario,fill=scenario))+
 geom_area(alpha=.12,position='identity')+geom_line(linewidth=.8)+
 geom_vline(data=s,aes(xintercept=predictive_median,colour=scenario),linetype=2,linewidth=.45,show.legend=FALSE)+
 geom_text(data=s,aes(x=295,y=y,label=annotation,colour=scenario),inherit.aes=FALSE,hjust=1,vjust=1,size=3.1,show.legend=FALSE)+
 facet_grid(dataset~model)+scale_colour_manual(values=colours)+scale_fill_manual(values=colours)+
 coord_cartesian(xlim=c(0,300),ylim=c(0,.062),expand=FALSE)+theme_classic(base_size=12)+
 theme(legend.position='bottom',strip.background=element_blank(),strip.text=element_text(face='bold'))+
 labs(x='Predicted number of species',y='Predictive density',colour=NULL,fill=NULL,
 caption='10 Ma. Scotese: ~32,000 km²; Straume: ~114,000 km². Dashed lines: medians.\nDisplay cropped at 300 species; statistics use the full distributions. Geological areas treated as effective habitat area.')
ggsave(file.path(plots,'richness-distributions.png'),p,width=12,height=8,dpi=180)
ggsave(file.path(plots,'richness-distributions.pdf'),p,width=12,height=8)
# Interval comparison makes the tails and central 50% easy to compare.
s$scenario <- factor(s$scenario,levels=c('Straume','Scotese'))
p <- ggplot(s,aes(y=scenario,colour=scenario))+
 geom_segment(aes(x=lower95,xend=upper95,yend=scenario),linewidth=.7)+
 geom_segment(aes(x=lower50,xend=upper50,yend=scenario),linewidth=4,alpha=.45)+
 geom_point(aes(x=predictive_median,shape='Median'),size=2.8)+
 geom_point(aes(x=predictive_mean,shape='Mean'),size=3)+
 facet_grid(dataset~model)+scale_colour_manual(values=colours)+scale_shape_manual(values=c(Median=16,Mean=4))+
 theme_classic(base_size=12)+theme(legend.position='bottom',strip.background=element_blank(),strip.text=element_text(face='bold'))+
 labs(x='Predicted number of species',y=NULL,colour=NULL,shape=NULL,
 caption='Thick bars: central 50%; thin bars: 95% predictive intervals.\nScotese: ~32,000 km²; Straume: ~114,000 km². Both geological areas are fixed scenarios at 10 Ma.')
ggsave(file.path(plots,'richness-intervals.png'),p,width=11,height=6,dpi=180)
ggsave(file.path(plots,'richness-intervals.pdf'),p,width=11,height=6)
