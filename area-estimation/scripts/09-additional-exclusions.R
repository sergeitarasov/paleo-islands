# Run from the paleo-islands project root.
library(readxl)
library(ggplot2)
library(ggrepel)
d <- as.data.frame(read_excel("data/richness/island-richness-area.xlsx"))
d$x <- log(d$area_km2)
d$y <- log(d$n_sp)
d$group <- factor(ifelse(d$region == "Caribbean", "Caribbean", "Other islands"))
keep <- !(d$region == "Australasian" & d$n_sp < 6 | d$island == "Réunion")
additional <- c("Tioman", "Langkawi", "Sibuyan", "Lombok", "Bali")
stopifnot(all(additional %in% d$island))
keep <- keep & !d$island %in% additional
sets <- list(Filtered = d[keep, ])
formulas <- list(Pooled = y ~ x, "Two intercepts" = y ~ x + group,
                "Two slopes" = y ~ x * group)
out <- "area-estimation/results/additional-exclusions"
plots.dir <- "area-estimation/plots/additional-exclusions"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(plots.dir, recursive = TRUE, showWarnings = FALSE)
write.csv(d[keep, ], file.path(out,"retained-islands.csv"), row.names=FALSE)
write.csv(d[!keep, ], file.path(out,"excluded-islands.csv"), row.names=FALSE)
scores <- data.frame(); predictions <- data.frame(); curves <- data.frame(); fits <- list()
coefficients <- data.frame()
for (dataset in names(sets)) {
  z <- sets[[dataset]]
  for (name in names(formulas)) {
    fit <- lm(formulas[[name]], data = z)
    fits[[paste(dataset,name)]] <- fit
    k <- attr(logLik(fit),"df"); n <- nrow(z)
    row <- data.frame(dataset=dataset,model=name,n=n,AICc=AIC(fit)+2*k*(k+1)/(n-k-1),
      adjusted_R2=summary(fit)$adj.r.squared,residual_SD=sigma(fit))
    for (i in seq_len(n)) {
      f <- lm(formulas[[name]],data=z[-i, ])
      pr <- as.numeric(predict(f,z[i, ,drop=FALSE]))
      predictions <- rbind(predictions,data.frame(dataset=dataset,model=name,island=z$island[i],group=z$group[i],
        observed=z$n_sp[i],log_prediction=pr,mean_prediction=exp(pr+sigma(f)^2/2)))
    }
    pred <- predictions[predictions$dataset==dataset & predictions$model==name, ]
    row$log_RMSE <- sqrt(mean((log(pred$observed)-pred$log_prediction)^2))
    row$count_RMSE <- sqrt(mean((pred$observed-pred$mean_prediction)^2))
    scores <- rbind(scores,row)
    ct <- cbind(coef(summary(fit)),confint(fit))
    coefficients <- rbind(coefficients,data.frame(dataset=dataset,model=name,term=rownames(ct),ct,row.names=NULL))
    for (group in levels(z$group)) {
      zg <- z[z$group==group, ]
      grid <- data.frame(x=seq(min(zg$x),max(zg$x),length.out=200),group=factor(group,levels=levels(z$group)))
      ci <- predict(fit,grid,interval="confidence")
      grid$area_km2 <- exp(grid$x); grid$richness <- exp(ci[,1])
      grid$lower <- exp(ci[,2]); grid$upper <- exp(ci[,3])
      grid$dataset <- dataset; grid$model <- name
      curves <- rbind(curves,grid)
    }
  }
}
for (dataset in names(sets)) {
  ix <- scores$dataset==dataset
  scores$delta_AICc[ix] <- scores$AICc[ix]-min(scores$AICc[ix])
}
write.csv(scores,file.path(out,"model-comparison.csv"),row.names=FALSE)
write.csv(coefficients,file.path(out,"coefficients.csv"),row.names=FALSE)
write.csv(predictions,file.path(out,"held-out-predictions.csv"),row.names=FALSE)
saveRDS(fits,file.path(out,"models.rds"))

# Pooled and grouped regressions, with all retained island names ----
cols <- c(Caribbean="#D55E00", "Other islands"="#0072B2")
for (model in names(formulas)) {
  line <- curves[curves$model==model, ]
  if (model == "Pooled") {
    grid <- data.frame(x=seq(min(sets$Filtered$x),max(sets$Filtered$x),length.out=300))
    ci <- predict(fits[["Filtered Pooled"]],grid,interval="confidence")
    line <- data.frame(area_km2=exp(grid$x),richness=exp(ci[,1]),lower=exp(ci[,2]),upper=exp(ci[,3]),group="Pooled")
  }
  p <- ggplot(sets$Filtered,aes(area_km2,n_sp,colour=group)) +
    geom_point(size=2.4) +
    geom_text_repel(aes(label=island),seed=26,max.overlaps=Inf,size=3.1,
      max.iter=100000,max.time=10,box.padding=.5,show.legend=FALSE) +
    scale_x_log10(breaks=c(100,1000,10000,100000,1000000),labels=c("100","1k","10k","100k","1m"),expand=expansion(mult=c(.12,.15))) +
    scale_y_log10(expand=expansion(mult=c(.16,.15))) +
    scale_colour_manual(values=cols) +
    labs(x=expression("Island area (km"^2*")"),y="Species richness",colour=NULL,fill=NULL,
      caption=paste(model,"; 39 islands. Bands: 95% confidence intervals for fitted median richness.")) +
    theme_classic(base_size=12) + theme(legend.position="bottom")
  if (model == "Pooled") {
    p <- p + geom_ribbon(data=line,aes(x=area_km2,ymin=lower,ymax=upper),inherit.aes=FALSE,fill="grey75",alpha=.25) +
      geom_line(data=line,aes(area_km2,richness),inherit.aes=FALSE,colour="grey25",linewidth=.8)
  } else {
    p <- p + geom_ribbon(data=line,aes(x=area_km2,ymin=lower,ymax=upper,fill=group),inherit.aes=FALSE,alpha=.12) +
      geom_line(data=line,aes(area_km2,richness,colour=group),inherit.aes=FALSE,linewidth=.8) + scale_fill_manual(values=cols)
  }
  stem <- c(Pooled="pooled-labelled", "Two intercepts"="shared-slope-labelled", "Two slopes"="separate-slopes-labelled")[[model]]
  ggsave(file.path(plots.dir,paste0(stem,".pdf")),p,width=12,height=9)
  ggsave(file.path(plots.dir,paste0(stem,".png")),p,width=12,height=9,dpi=180)
}

report <- c("# Additional island exclusions", "",
  "Starting with the previous exclusions (Australasian richness below six species and Reunion), additionally excluded Tioman, Langkawi, Sibuyan, Lombok and Bali. This leaves 39 records: 16 Caribbean and 23 other islands. The source workbook and earlier results are unchanged.", "",
  "Compared pooled power regression, common slope with group-specific intercepts, and separate slopes. All use Gaussian residuals on natural-log richness. Leave-one-out predictions use identical 39-island test sets for the three candidates.", "",
  "| Model | Adjusted R² | Log CV RMSE | Count CV RMSE | Delta AICc |",
  "|---|---:|---:|---:|---:|")
for (i in seq_len(nrow(scores))) report <- c(report,sprintf("| %s | %.3f | %.3f | %.2f | %.2f |",scores$model[i],scores$adjusted_R2[i],scores$log_RMSE[i],scores$count_RMSE[i],scores$delta_AICc[i]))
report <- c(report,"",
  "Separate slopes are preferred on this selected dataset: adjusted R² 0.893 versus 0.847 for the pooled model; log-scale prediction error is 15.7% lower, and AICc improves by 11.18. The Caribbean area exponent is 0.358, compared with approximately 0.646 for other islands. These findings are conditional on the specified exclusions.", "",
  "Compare models on these same 39 records; do not compare AICc with fits on the earlier 44 or 61 records. Successive exclusions based on inspecting richness and the fitted plot can make apparent fit and predictive scores optimistic for the original population. An independent inventory-quality justification is needed before adopting this filtered calibration as primary. Low-richness Caribbean records remain included. This analysis does not revise Bayesian area estimates.", "",
  "Run area-estimation/scripts/09-additional-exclusions.R from the project root. Outputs: results/additional-exclusions and plots/additional-exclusions.")
writeLines(report,"area-estimation/reports/additional-exclusions.md")
print(scores)
print(coefficients)
