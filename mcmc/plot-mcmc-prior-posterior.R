# Run from the paleo-islands project root ----
# Same final chains and warmup as analyse-mcmc.R.
draws6 <- NULL
draws7 <- NULL
for (chain in 1:4) {
  x <- readRDS(paste0("mcmc/results/roihu-5000-20261023/gr6-chain-",chain,".rds"))
  draws6 <- rbind(draws6,x$rates[1001:5000,])
  x <- readRDS(paste0("mcmc/results/roihu-5000-20261023/gr7-chain-",chain,".rds"))
  draws7 <- rbind(draws7,x$rates[1001:5000,])
}
stopifnot(all(is.finite(draws6)),all(draws6>0),all(is.finite(draws7)),all(draws7>0))
parameters <- c("lambda010101.1","mu01.1","q0108.1","q0516.1",
  "lambda060606.1","lambda080103.1","lambda070102.1","q0107.1",
  "mu02.1","mu02.1","mu06.1")
labels <- c("[sA]","[xA]","[dAM, dEM]","[dSR, dAS]","[sR, sS]",
  "[scAM, scEM, scSR, scAS]","[scA, scE, sE, sM, sU]","[dAE, dAU, dEA, dEU]",
  "Group 6\n[xE, xM, xU, xR, xS]","Group 7\n[xE, xM, xU]","Group 7\n[xR, xS]")
colours <- c("#0072B2","#D55E00")

# Density in z=log10(rate), with rates labelled in their original units ----
# Transform the prior density too: f_Z(z) = log(10)*10^z*f_rate(10^z).
# This changes only the display; the sampling prior remains Exp(rate=1).
z <- seq(-5,1,length.out=2048)
prior <- log(10)*10^z*dexp(10^z,rate=1)
summary <- data.frame()
dir.create("mcmc/plots",showWarnings=FALSE)
dir.create("mcmc/summaries",showWarnings=FALSE)
for (format in c("pdf","png")) {
  file <- paste0("mcmc/plots/gr6-gr7-prior-posterior.",format)
  if(format=="pdf") pdf(file,width=10.5,height=10,useDingbats=FALSE)
  if(format=="png") png(file,width=10.5,height=10,units="in",res=300)
  par(mfrow=c(4,3),mar=c(3.2,3.4,2.9,.6),mgp=c(1.9,.6,0),
      oma=c(2.5,.5,0,0),bty="l",tcl=-.2,cex=.8)
  for (p in seq_along(parameters)) {
    models <- c(6,7)
    if(p==9) models <- 6
    if(p>=10) models <- 7
    pooled <- numeric()
    for(model in models) {
      v <- if(model==6) draws6[,parameters[p]] else draws7[,parameters[p]]
      pooled <- c(pooled,log10(v))
    }
    bandwidth <- bw.nrd0(pooled)
    curves <- list()
    ymax <- max(prior)
    for(model in models) {
      v <- if(model==6) draws6[,parameters[p]] else draws7[,parameters[p]]
      d <- density(log10(v),bw=bandwidth,from=-5,to=1,n=2048)
      curves[[as.character(model)]] <- d
      ymax <- max(ymax,d$y)
      if(format=="pdf") {
        ci <- quantile(v,c(.025,.975),type=7)
        summary <- rbind(summary,data.frame(group=model,parameter=parameters[p],
          mean=mean(v),median=median(v),lower=ci[1],upper=ci[2],
          prior.median=qexp(.5),prior.lower=qexp(.025),prior.upper=qexp(.975),
          posterior.to.prior.width=(ci[2]-ci[1])/(qexp(.975)-qexp(.025)),
          fraction.outside.display=mean(v<1e-5 | v>10)))
      }
    }
    plot(NA,xlim=c(-5,1),ylim=c(0,ymax*1.08),xaxt="n",
         xlab="",ylab="Density",main=labels[p],cex.main=.95)
    axis(1,at=c(-5,-3,-1,1),labels=expression(10^{-5},10^{-3},10^{-1},10^1))
    polygon(c(z,rev(z)),c(prior,rep(0,length(z))),col="grey93",border=NA)
    for(model in models) {
      d <- curves[[as.character(model)]]
      colour <- colours[model-5]
      polygon(c(d$x,rev(d$x)),c(d$y,rep(0,length(d$y))),
              col=adjustcolor(colour,alpha.f=.25),border=NA)
      lines(d$x,d$y,col=colour,lwd=1.6)
    }
    lines(z,prior,col="grey35",lty=2,lwd=1.5)
  }
  plot.new()
  legend("center",legend=c("Prior: Exp(1)","Posterior: group 6","Posterior: group 7"),
         col=c("grey35",colours),lty=c(2,1,1),lwd=2,bty="n",cex=1,y.intersp=1.4)
  mtext(expression(paste("Rate (",Ma^{-1},"; logarithmic scale)")),side=1,outer=TRUE,line=.7,cex=.95)
  dev.off()
}
rownames(summary) <- NULL
write.csv(summary,"mcmc/summaries/prior-posterior-comparison.csv",row.names=FALSE)
