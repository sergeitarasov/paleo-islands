# Run from the paleo-islands project root.
library(readxl)
library(ggplot2)
library(ggrepel)
d <- as.data.frame(read_excel("data/richness/island-richness-area.xlsx"))
d$x <- log(d$area_km2)
d$y <- log(d$n_sp)
d$group <- factor(ifelse(d$region == "Caribbean", "Caribbean", "Other islands"))
keep <- !(d$region == "Australasian" & d$n_sp < 6 | d$island == "Réunion")
sets <- list(All = d, Filtered = d[keep, ])
formulas <- list(Pooled = y ~ x, "Two intercepts" = y ~ x + group,
                "Two slopes" = y ~ x * group)
out <- "area-estimation/results/caribbean-versus-rest"
plots.dir <- "area-estimation/plots/caribbean-versus-rest"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(plots.dir, recursive = TRUE, showWarnings = FALSE)
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

# Show separate-slope fits and shared-slope fits, not just the selected winner.
points <- rbind(data.frame(d,dataset="All"),data.frame(d[keep, ],dataset="Filtered"))
cols <- c(Caribbean="#D55E00","Other islands"="#0072B2")
for (model in c("Two intercepts","Two slopes")) {
  line <- curves[curves$model==model, ]
  p <- ggplot(points,aes(area_km2,n_sp,colour=group)) + geom_point(size=2) +
    geom_ribbon(data=line,aes(x=area_km2,ymin=lower,ymax=upper,fill=group),inherit.aes=FALSE,alpha=.12) +
    geom_line(data=line,aes(area_km2,richness,colour=group),inherit.aes=FALSE,linewidth=.8) +
    facet_wrap(~dataset) + scale_x_log10(breaks=c(100,1000,10000,100000,1000000),labels=c("100","1k","10k","100k","1m")) +
    scale_y_log10() + scale_colour_manual(values=cols) + scale_fill_manual(values=cols) +
    labs(x=expression("Island area (km"^2*")"),y="Species richness",colour=NULL,fill=NULL,
      caption="Lines: fitted median richness; bands: 95% confidence intervals. All: 61 islands; filtered: 44.") +
    theme_classic(base_size=12) + theme(legend.position="bottom")
  stem <- ifelse(model=="Two slopes","separate-slopes","shared-slope")
  ggsave(file.path(plots.dir,paste0(stem,".pdf")),p,width=11,height=6)
  ggsave(file.path(plots.dir,paste0(stem,".png")),p,width=11,height=6,dpi=180)
}
# Labelled filtered plot for inspecting particular islands.
line <- curves[curves$model=="Two slopes" & curves$dataset=="Filtered", ]
p <- ggplot(sets$Filtered,aes(area_km2,n_sp,colour=group)) + geom_point(size=2.3) +
  geom_ribbon(data=line,aes(x=area_km2,ymin=lower,ymax=upper,fill=group),inherit.aes=FALSE,alpha=.12) +
  geom_line(data=line,aes(area_km2,richness,colour=group),inherit.aes=FALSE,linewidth=.8) +
  geom_text_repel(aes(label=island),seed=26,max.overlaps=Inf,size=3,max.time=10,show.legend=FALSE) +
  scale_x_log10() + scale_y_log10() + scale_colour_manual(values=cols) + scale_fill_manual(values=cols) +
  labs(x=expression("Island area (km"^2*")"),y="Species richness",colour=NULL,fill=NULL) +
  theme_classic(base_size=12) + theme(legend.position="bottom")
ggsave(file.path(plots.dir,"filtered-labelled.pdf"),p,width=12,height=9)
ggsave(file.path(plots.dir,"filtered-labelled.png"),p,width=12,height=9,dpi=180)

report <- c("# Caribbean versus other islands","",
"Three lognormal power models were compared: one pooled relationship; two group intercepts with a common slope; separate group intercepts and slopes. Grouping is Caribbean versus Australasian + Indian Ocean. Leave-one-island-out tests use the same observations for all models within each dataset.","",
"| Data | Model | Adjusted R² | Log CV RMSE | Count CV RMSE | Delta AICc |",
"|---|---|---:|---:|---:|---:|")
for(i in seq_len(nrow(scores))) report <- c(report,sprintf("| %s | %s | %.3f | %.3f | %.2f | %.2f |",scores$dataset[i],scores$model[i],scores$adjusted_R2[i],scores$log_RMSE[i],scores$count_RMSE[i],scores$delta_AICc[i]))
report <- c(report,"",
"## Interpretation", "",
"For the filtered data, a common slope with two intercepts is preferred: log-scale prediction error improves by 16.0%, adjusted R² increases from 0.691 to 0.787, and AICc improves by 15.11. Separate slopes are not needed (AICc 2.42 worse than shared slope). However, count-scale RMSE increases from 29.38 to 33.95, so improvement is specific to proportional/log-scale prediction rather than every criterion.", "",
"For all islands, group models improve log-scale RMSE by only 0.7–1.4%, and AICc differences are below one. Thus the split is useful for the filtered sample but is not strongly supported as a universal difference in the complete data.", "",
"Full data: 16 Caribbean and 45 other islands. Filtered data: the same 16 Caribbean and 28 other islands. The filter removes 16 Australasian records below six species and Reunion, leaving 44 records. Compare AICc only within a dataset. All models use Gaussian errors on log richness. Plot bands describe uncertainty in fitted median richness, not new-island prediction intervals.","",
"The filtered comparison deliberately retains low-richness Caribbean records while excluding low-richness Australasian records. This response-based selection can create or exaggerate a group contrast. Differences may also reflect geology, source/inventory differences and isolation. The split is not proof of a causal regional effect; assess the all-data result alongside the filtered one. These exploratory comparisons reuse the dataset to select hypotheses and are not independent validation.","",
"Previous Bayesian area estimates are unchanged. Scripts use the authoritative workbook without editing it. Run area-estimation/scripts/08-caribbean-versus-rest.R from the project root.")
writeLines(report,"area-estimation/reports/caribbean-versus-rest.md")
print(scores)
print(coefficients[coefficients$model=="Two slopes", ])
