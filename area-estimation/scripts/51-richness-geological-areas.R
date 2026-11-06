# Richness predicted from fixed approximate geological areas at 10 Ma.
# Same models as the biological area analysis; no refitting.
models27 <- readRDS('area-estimation/results/bayesian-selected27-slopes-10ma/regression-posteriors.rds')
models22 <- readRDS('area-estimation/results/bayesian-no-malay-22-10ma/regression-posteriors.rds')
fits <- list('27 islands'=models27, '22 islands'=models22)
geology <- data.frame(scenario=c('Scotese','Straume'),area_km2=c(32000,114000))
out <- 'area-estimation/results/richness-geological-areas'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
summary <- data.frame()
for(dataset in names(fits)) for(model in c('simple_prior1','separate_prior1')) {
 fit <- fits[[dataset]][[model]]
 draws <- fit$draws
 p <- ncol(draws)-1
 for(g in 1:nrow(geology)) {
  new <- data.frame(x=log(geology$area_km2[g]/10000),group=factor('Other',levels=c('Caribbean','Other')))
  X <- model.matrix(fit$formula,new)
  mu <- as.vector(draws[,1:p] %*% as.vector(X))
  sigma <- exp(draws[,p+1])
  # Invert the mixture CDF of log-richness. Includes coefficient uncertainty
  # and residual island-to-island variation, using all retained MCMC draws.
  probs <- c(.025,.25,.5,.75,.975)
  q <- numeric(5)
  for(j in 1:5) {
   cdf_difference <- function(y) mean(pnorm(y,mu,sigma))-probs[j]
   q[j] <- exp(uniroot(cdf_difference,interval=c(min(mu-12*sigma),max(mu+12*sigma)),tol=1e-9)$root)
  }
  expected <- exp(mu+sigma^2/2)
  summary <- rbind(summary,data.frame(dataset,model,scenario=geology$scenario[g],area_km2=geology$area_km2[g],
    predictive_mean=mean(expected),predictive_median=q[3],lower95=q[1],upper95=q[5],lower50=q[2],upper50=q[4],
    conditional_mean_lower95=unname(quantile(expected,.025)),conditional_mean_upper95=unname(quantile(expected,.975))))
 }
}
write.csv(summary,file.path(out,'predicted-richness.csv'),row.names=FALSE)
write.csv(geology,file.path(out,'geological-inputs.csv'),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,'session-info.txt'))
print(summary)
