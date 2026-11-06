# Run from paleo-islands. Requested geographic subset; original workbook unchanged.
library(readxl)
library(ggplot2)
library(ggrepel)
d <- as.data.frame(read_excel('data/richness/island-richness-area.xlsx',sheet='Island richness'))
d <- d[d$region=='Caribbean' | d$island %in% c('Réunion','Mauritius','Madagascar','New Caledonia'),c('island','region','n_sp','area_km2')]
# User-supplied NZ richness; provisional national land area from Australian DFAT.
# https://www.dfat.gov.au/geo/new-zealand/new-zealand-country-brief
# This differs from the Global Islands v3 polygon source used for the other records.
d <- rbind(d,data.frame(island='New Zealand',region='Australasian',n_sp=18,area_km2=263310))
d$group <- factor(ifelse(d$region=='Caribbean','Caribbean','Other selected islands'))
d$x <- log(d$area_km2);d$y <- log(d$n_sp)
out <- 'area-estimation/results/selected-21'
plots.dir <- 'area-estimation/plots/selected-21'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
dir.create(plots.dir,recursive=TRUE,showWarnings=FALSE)
write.csv(d,file.path(out,'input.csv'),row.names=FALSE)
forms <- list(pooled=y~x,group_intercepts=y~x+group,group_slopes=y~x*group)
scores <- data.frame();coefficients <- data.frame()
fits <- list()
for(name in names(forms)) {
 f <- lm(forms[[name]],data=d);fits[[name]] <- f
 pred <- numeric(nrow(d))
 for(i in seq_len(nrow(d))) pred[i] <- predict(lm(forms[[name]],data=d[-i, ]),d[i, ])
 k <- attr(logLik(f),'df')
 scores <- rbind(scores,data.frame(model=name,n=nrow(d),R2=summary(f)$r.squared,
  adjusted_R2=summary(f)$adj.r.squared,AICc=AIC(f)+2*k*(k+1)/(nrow(d)-k-1),LOO_log_RMSE=sqrt(mean((d$y-pred)^2))))
 coefficients <- rbind(coefficients,data.frame(model=name,term=names(coef(f)),estimate=coef(f),confint(f),p=summary(f)$coefficients[,4]))
}
sensitivity <- data.frame()
subsets <- list(all=d,without_Madagascar=d[d$island!='Madagascar', ],Caribbean=d[d$group=='Caribbean', ],other_five=d[d$group!='Caribbean', ])
for(name in names(subsets)) {
 a <- subsets[[name]];f <- lm(y~x,a);ci<-confint(f)
 sensitivity <- rbind(sensitivity,data.frame(subset=name,n=nrow(a),slope=coef(f)[2],lower=ci[2,1],upper=ci[2,2],R2=summary(f)$r.squared,p=summary(f)$coefficients[2,4]))
}
f <- fits$pooled
r <- cor.test(d$area_km2,d$n_sp,method='spearman',exact=FALSE)
write.csv(scores,file.path(out,'models.csv'),row.names=FALSE)
write.csv(coefficients,file.path(out,'coefficients.csv'),row.names=FALSE)
write.csv(sensitivity,file.path(out,'sensitivity.csv'),row.names=FALSE)
capture.output(summary(f),confint(f),r,file=file.path(out,'statistics.txt'))
# Show pooled regression only, with labels for every observation.
grid <- data.frame(x=seq(min(d$x),max(d$x),length.out=200))
ci <- predict(f,grid,interval='confidence')
grid$area<-exp(grid$x);grid$fit<-exp(ci[,1]);grid$lower<-exp(ci[,2]);grid$upper<-exp(ci[,3])
p <- ggplot(d,aes(area_km2,n_sp))+
 geom_ribbon(data=grid,aes(x=area,ymin=lower,ymax=upper),inherit.aes=FALSE,fill='grey88')+
 geom_line(data=grid,aes(area,fit),inherit.aes=FALSE,colour='grey25')+
 geom_point(aes(colour=group),size=2.5)+
 geom_text_repel(aes(label=island,colour=group),seed=27,size=3.2,max.overlaps=Inf,min.segment.length=0)+
 scale_x_log10()+scale_y_log10()+scale_colour_manual(values=c('#0072B2','#D55E00'))+
 theme_classic(base_size=12)+theme(legend.position='bottom')+
 labs(x=expression('Area (km'^2*')'),y='Species richness',colour=NULL,
 caption='Pooled fit to the requested 21 records; band: 95% confidence interval for fitted median richness.')
ggsave(file.path(plots.dir,'regression.png'),p,width=11,height=8,dpi=180)
ggsave(file.path(plots.dir,'regression.pdf'),p,width=11,height=8)
report <- c('# Requested 21-island subset','',
'Included all 16 Caribbean records plus Réunion, Mauritius, Madagascar, New Caledonia and New Zealand. Mauritius retains the workbook value of 6 species. NZ uses the supplied 18 species and 263,310 km2 land area from https://www.dfat.gov.au/geo/new-zealand/new-zealand-country-brief. NZ is a country/island-group observation; area source is not yet harmonized with Global Islands v3. The workbook is unchanged.','',
'Exploratory lognormal regression. AICc comparisons only use models fitted to the same 21 observations. LOO validates islands, not entire geographic regions. The five non-Caribbean observations cover different climates/geological histories and do not define a single homogeneous reference group.','',capture.output(scores),'',capture.output(sensitivity),'',capture.output(r),'',
'A stronger association within this selected subset does not establish that it is the correct calibration for Sunken Islands. Selection and inventory uncertainty are not propagated, and forward association alone does not establish accurate inverse area recovery. No Sunken Islands area estimate was recalculated.')
writeLines(report,'area-estimation/reports/selected-21.md')
print(scores);print(sensitivity);print(r)
