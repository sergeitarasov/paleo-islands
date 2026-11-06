# Run from paleo-islands. Only the 32 mapped islands smaller than Timor.
library(ggplot2)
library(ggrepel)
d <- read.csv('data/richness/maps/islands-smaller-than-timor.csv')
d$log_area <- log(d$area_km2)
d$log_richness <- log(d$n_sp)
f <- lm(log_richness ~ log_area,data=d)
r <- cor.test(d$log_area,d$log_richness)
rho <- cor.test(d$area_km2,d$n_sp,method='spearman',exact=FALSE)
curve <- data.frame(log_area=seq(min(d$log_area),max(d$log_area),length.out=200))
ci <- predict(f,curve,interval='confidence')
curve$area <- exp(curve$log_area)
curve$fit <- exp(ci[,1]);curve$lower <- exp(ci[,2]);curve$upper <- exp(ci[,3])
p <- ggplot(d,aes(area_km2,n_sp))+
 geom_ribbon(data=curve,aes(x=area,ymin=lower,ymax=upper),inherit.aes=FALSE,fill='grey88')+
 geom_line(data=curve,aes(area,fit),inherit.aes=FALSE,colour='grey30',linewidth=.8)+
 geom_point(colour='#D73027',size=2.5)+
 geom_text_repel(aes(label=island),colour='#B52720',size=3,seed=27,max.overlaps=Inf,
  min.segment.length=0,segment.colour='grey65',box.padding=.4)+
 scale_x_log10()+scale_y_log10(breaks=c(1,2,5,10,20,50))+
 theme_classic(base_size=12)+labs(x=expression('Area (km'^2*')'),y='Species richness',
 caption='Fit uses only the 32 highlighted islands. Band: 95% confidence interval for fitted median richness.')
ggsave('area-estimation/plots/pilot-regression/small-malay-only-regression.png',p,width=11,height=8,dpi=180)
ggsave('area-estimation/plots/pilot-regression/small-malay-only-regression.pdf',p,width=11,height=8)
capture.output(summary(f),confint(f),r,rho,file='area-estimation/results/pilot-regression/small-malay-only-statistics.txt')
