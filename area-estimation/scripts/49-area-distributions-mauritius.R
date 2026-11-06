# Run from paleo-islands. Existing posterior distributions; no refitting.
curves <- read.csv('area-estimation/results/bayesian-selected27-slopes-10ma/area-distributions.csv')
curves$dataset <- '27 islands'
extra <- read.csv('area-estimation/results/bayesian-no-malay-22-10ma/area-distributions.csv')
extra$dataset <- '22 islands'
curves <- rbind(curves, extra)
datasets <- c('27 islands', '22 islands')
islands <- read.csv('area-estimation/results/isolation/regression-input.csv')
mauritius_area <- islands$area_km2[islands$island == 'Mauritius']
stopifnot(length(mauritius_area) == 1)
models <- c('simple_prior1', 'separate_prior1')
titles <- c('Pooled model', 'Separate slopes: non-Caribbean')
colours <- c('#0072B2', '#D55E00')
out <- 'area-estimation/results/two-datasets-mauritius-10ma'
plots <- 'area-estimation/plots/two-datasets-mauritius-10ma'
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(plots, recursive = TRUE, showWarnings = FALSE)
statistics <- data.frame()
for(dataset in datasets) for(i in 1:2) {
 d <- curves[curves$dataset == dataset & curves$model == models[i] & curves$prior == 'log_uniform_10_1million', ]
 mass <- d$density_log10_area * .01
 stopifnot(abs(sum(mass)-1)<1e-8)
 q <- 10^approx(cumsum(mass), log10(d$area_km2), xout = c(.025,.25,.5,.75,.975), ties = 'ordered')$y
 values <- c(q[1], q[2], q[3], sum(d$area_km2*mass), q[4], q[5])
 labels <- c('95% lower', '50% lower', 'Median', 'Mean', '50% upper', '95% upper')
 for(j in 1:6) {
  nearest <- which.min(abs(log(islands$area_km2 / values[j])))
  statistics <- rbind(statistics, data.frame(dataset = dataset, model = models[i], statistic = labels[j], area_km2 = values[j],
   comparison_island = islands$island[nearest], island_area_km2 = islands$area_km2[nearest],
   mauritius_area_km2 = mauritius_area, mauritius_equivalents = values[j]/mauritius_area, size_ratio = islands$area_km2[nearest]/values[j]))
 }
}
write.csv(statistics, file.path(out, 'statistics-mauritius-comparisons.csv'), row.names = FALSE)
write.csv(islands[,c('island','area_km2','source_of_area')], file.path(out,'comparison-island-catalogue.csv'), row.names=FALSE)

# Panels: density above; numerical estimates and island analogues below.
# Distribution is per unit log10 area. Intervals are equal-tailed, not HPD.
for(format in c('png','pdf')) {
 if(format=='png') png(file.path(plots,'area-distributions.png'), width=2600,height=3200,res=200,type='cairo')
 if(format=='pdf') pdf(file.path(plots,'area-distributions.pdf'), width=13,height=16)
 layout(matrix(1:8,4,2,byrow=TRUE), heights=c(3.2,1.8,3.2,1.8))
 for(dataset in datasets) {
 par(mar=c(4.3,4.2,3,1),family='sans')
 for(i in 1:2) {
  d <- curves[curves$dataset==dataset & curves$model==models[i] & curves$prior=='log_uniform_10_1million' & curves$area_km2>=10 & curves$area_km2<=1e6, ]
  z <- log10(d$area_km2)
  density <- d$density_log10_area
  density[c(1,length(density))] <- density[c(1,length(density))]*2
  s <- statistics[statistics$dataset==dataset & statistics$model==models[i], ]
  plot(z,density,type='n',xlim=c(1,6),ylim=c(0,.78),axes=FALSE,xlab=expression('Effective area (km'^2*')'),ylab='Density per log10 area',main=paste(dataset, titles[i], sep=" | "))
  axis(1,at=1:6,labels=expression(10^1,10^2,10^3,10^4,10^5,10^6))
  axis(2,las=1); box(bty='l')
  polygon(c(z,rev(z)),c(density,rep(0,length(z))),col='grey95',border=NA)
  for(bounds in list(c(1,6),c(2,5))) {
   lo <- log10(s$area_km2[bounds[1]]); hi <- log10(s$area_km2[bounds[2]])
   x <- c(lo,z[z>lo & z<hi],hi)
   y <- approx(z,density,xout=x)$y
   alpha <- if(bounds[1]==1) .18 else .42
   polygon(c(x,rev(x)),c(y,rep(0,length(x))),col=adjustcolor(colours[i],alpha.f=alpha),border=NA)
  }
  lines(z,density,col=colours[i],lwd=2)
  abline(v=log10(s$area_km2[3]),lty=2,col=colours[i],lwd=1.5)
  abline(v=log10(s$area_km2[4]),lty=3,col=colours[i],lwd=1.5)
  segments(log10(32000),0,log10(32000),.66,col='#444444',lty=4)
  segments(log10(114000),0,log10(114000),.73,col='#777777',lty=4)
  text(log10(32000),.675,'Scotese ~32,000',cex=.8)
  text(log10(114000),.745,'Straume ~114,000',cex=.8)
  legend('topleft',legend=c('95% interval','50% interval','Median','Mean'),
   fill=c(adjustcolor(colours[i],alpha.f=.18),adjustcolor(colours[i],alpha.f=.42),NA,NA),
   border=NA,lty=c(NA,NA,2,3),col=colours[i],bty='n',cex=.8)
 }
 par(mar=c(1,1,1,1))
 for(i in 1:2) {
  s <- statistics[statistics$dataset==dataset & statistics$model==models[i], ]
  plot.new(); plot.window(xlim=c(0,1),ylim=c(0,1))
  text(.02,.94,'Statistic',adj=0,font=2,cex=.9)
  text(.34,.94,expression('Area (km'^2*')'),adj=1,font=2,cex=.9)
  text(.40,.94,'Multiples of Mauritius area',adj=0,font=2,cex=.85)
  segments(.02,.88,.99,.88,col='grey75')
  for(j in 1:6) {
   y <- .8-(j-1)*.125
   text(.02,y,s$statistic[j],adj=0,cex=.9)
   text(.34,y,format(round(s$area_km2[j]),big.mark=',',trim=TRUE),adj=1,cex=.9)
   label <- sprintf('%.3g x Mauritius', s$mauritius_equivalents[j])
   text(.40,y,label,adj=0,cex=.85)
  }
  text(.02,.04,sprintf('10 Ma; 1 Mauritius = %s km²; equal-tailed intervals.', format(round(mauritius_area),big.mark=',',trim=TRUE)),adj=0,cex=.72,col='grey35')
 }
 }
 dev.off()
}
print(statistics)
