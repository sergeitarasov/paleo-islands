# Same area distributions, shown on linear axes. Run from project root.
a <- read.csv('area-estimation/results/bayesian-selected27-slopes-10ma/area-distributions.csv')
b <- read.csv('area-estimation/results/bayesian-no-malay-22-10ma/area-distributions.csv')
a$dataset <- '27 islands'
b$dataset <- '22 islands'
curves <- rbind(a,b)
statistics <- read.csv('area-estimation/results/two-datasets-mauritius-10ma/statistics-mauritius-comparisons.csv')
models <- c('simple_prior1','separate_prior1')
titles <- c('Pooled','Separate slopes: non-Caribbean')
colours <- c('#0072B2','#D55E00')
plots <- 'area-estimation/plots/two-datasets-mauritius-10ma'
# f_A(a) = f_log10(A)(log10(a)) / (a * ln(10)).
# Correct endpoint quadrature weights before displaying density.
d <- curves[curves$prior=='log_uniform_10_1million' & curves$area_km2>=10 & curves$area_km2<=1e6, ]
d$density_log10_area[d$area_km2 %in% c(10,1e6)] <- 2*d$density_log10_area[d$area_km2 %in% c(10,1e6)]
d$density_area <- d$density_log10_area/(d$area_km2*log(10))
height <- max(d$density_area)*1.12
for(view in c('wide','zoom')) for(format in c('png','pdf')) {
 limit <- if(view=='wide') 200000 else 10000
 filename <- paste0('area-distributions-linear-',view,'.',format)
 if(format=='png') png(file.path(plots,filename),width=2400,height=1800,res=200,type='cairo')
 if(format=='pdf') pdf(file.path(plots,filename),width=12,height=9)
 par(mfrow=c(2,2),mar=c(4.3,5,3,1),oma=c(2,0,0,0))
 for(dataset in c('27 islands','22 islands')) for(i in 1:2) {
  curve <- d[d$dataset==dataset & d$model==models[i], ]
  s <- statistics[statistics$dataset==dataset & statistics$model==models[i], ]
  x <- curve$area_km2
  y <- curve$density_area
  plot(x,y,type='n',xlim=c(0,limit),ylim=c(0,height),axes=FALSE,
   xlab=expression('Effective area (km'^2*')'),ylab=expression('Density per km'^2),main=paste(dataset,titles[i],sep=' | '))
  ticks <- pretty(c(0,limit),n=5)
  axis(1,at=ticks,labels=format(ticks,big.mark=',',scientific=FALSE,trim=TRUE),cex.axis=.8)
  axis(2,las=1,cex.axis=.8); box(bty='l')
  polygon(c(x,rev(x)),c(y,rep(0,length(x))),col='grey94',border=NA)
  for(bounds in list(c(1,6),c(2,5))) {
   lo <- s$area_km2[bounds[1]]; hi <- s$area_km2[bounds[2]]
   xx <- c(lo,x[x>lo & x<hi],hi)
   yy <- approx(x,y,xout=xx)$y
   polygon(c(xx,rev(xx)),c(yy,rep(0,length(xx))),border=NA,
    col=adjustcolor(colours[i],alpha.f=if(bounds[1]==1) .18 else .42))
  }
  lines(x,y,col=colours[i],lwd=2)
  abline(v=s$area_km2[3],lty=2,col=colours[i])
  abline(v=s$area_km2[4],lty=3,col=colours[i])
  if(view=='wide') {
   segments(32000,0,32000,height*.55,lty=4,col='grey40')
   segments(114000,0,114000,height*.55,lty=4,col='grey40')
   text(32000,height*.59,'Scotese\n~32,000',cex=.75)
   text(114000,height*.59,'Straume\n~114,000',cex=.75)
  }
  label <- sprintf('Mean: %s; median: %s\n50%%: %s–%s\n95%%: %s–%s km²',
   format(round(s$area_km2[4]),big.mark=',',trim=TRUE),format(round(s$area_km2[3]),big.mark=',',trim=TRUE),
   format(round(s$area_km2[2]),big.mark=',',trim=TRUE),format(round(s$area_km2[5]),big.mark=',',trim=TRUE),
   format(round(s$area_km2[1]),big.mark=',',trim=TRUE),format(round(s$area_km2[6]),big.mark=',',trim=TRUE))
  text(limit*.97,height*.94,label,adj=c(1,1),cex=.78)
  legend('right',legend=c('95%','50%','Median','Mean'),fill=c(adjustcolor(colours[i],.18),adjustcolor(colours[i],.42),NA,NA),border=NA,lty=c(NA,NA,2,3),col=colours[i],bty='n',cex=.75)
  mass <- curve$density_log10_area*.01
  mass[c(1,length(mass))] <- mass[c(1,length(mass))]/2
  shown <- sum(mass[x<=limit])/sum(mass)
  mtext(sprintf('Displayed range contains approximately %.1f%% of probability',100*shown),side=3,line=.25,cex=.7)
 }
 mtext('10 Ma. Display cropped; distributions and numerical summaries retain the full 10–1,000,000 km² support.',outer=TRUE,side=1,line=.5,cex=.8)
 dev.off()
}
