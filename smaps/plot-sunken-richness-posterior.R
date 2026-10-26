# Posterior Sunken-island richness ----
# Run from the paleo-islands project root after mapping finishes.
library(divextra)
library(ggplot2)
maps.dir <- "smaps/maps/posterior-aic-saved-asr"
plots.dir <- "smaps/plots/posterior-aic-saved-asr"
summaries.dir <- "smaps/summaries/posterior-aic-saved-asr"
dir.create(plots.dir, recursive = TRUE, showWarnings = FALSE)
dir.create(summaries.dir, recursive = TRUE, showWarnings = FALSE)

# Both observed lineages and extinct side lineages count ----
nesovin <- readRDS(file.path(maps.dir, "smaps-nesovin-full.rds"))
nesosisyphus <- readRDS(file.path(maps.dir, "smaps-nesosisyphus-full.rds"))
stopifnot(length(nesovin) == length(nesosisyphus))
for (i in seq_along(nesovin)) {
  stopifnot(identical(nesovin[[i]]$posterior.draw$posterior.row,
                      nesosisyphus[[i]]$posterior.draw$posterior.row),
            identical(nesovin[[i]]$posterior.draw$model,
                      nesosisyphus[[i]]$posterior.draw$model))
}
height.nesovin <- max(nesovin[[1]]$lineages$start.age)
height.nesosisyphus <- max(nesosisyphus[[1]]$lineages$start.age)
# Use one common grid; do not infer zero richness before a clade's root.
times <- sort(unique(c(seq(0, max(height.nesovin, height.nesosisyphus), by=.1),
                       height.nesovin, height.nesosisyphus, 9)))
ltt.nesovin <- lineages.through.time.simmap(nesovin,
  times=times[times <= height.nesovin], states=c("A.S", "S", "S.R"), quantile.type=7)
ltt.nesosisyphus <- lineages.through.time.simmap(nesosisyphus,
  times=times[times <= height.nesosisyphus], states=c("A.S", "S", "S.R"), quantile.type=7)

# Sum matched histories FIRST, then summarize the combined distribution ----
common.times <- times[times <= min(height.nesovin, height.nesosisyphus)]
counts.combined <- ltt.nesovin$counts[match(common.times,ltt.nesovin$times),] +
  ltt.nesosisyphus$counts[match(common.times,ltt.nesosisyphus$times),]
combined <- data.frame(time=common.times,mean=rowMeans(counts.combined),
  median=NA_real_,lower=NA_real_,upper=NA_real_,n.maps=length(nesovin))
for (i in seq_along(common.times)) {
  combined$median[i] <- median(counts.combined[i,])
  combined$lower[i] <- quantile(counts.combined[i,], .025, type=7)
  combined$upper[i] <- quantile(counts.combined[i,], .975, type=7)
}
summary.nesovin <- ltt.nesovin$summary
summary.nesovin$clade <- "Nesovinsonia"
summary.nesosisyphus <- ltt.nesosisyphus$summary
summary.nesosisyphus$clade <- "Nesosisyphus"
combined$clade <- "Combined"
summary <- rbind(summary.nesovin,summary.nesosisyphus,combined)
summary$clade <- factor(summary$clade, levels=c("Nesovinsonia","Nesosisyphus","Combined"))
write.csv(summary,file.path(summaries.dir,"richness-through-time.csv"),row.names=FALSE)
saveRDS(list(nesovin=ltt.nesovin,nesosisyphus=ltt.nesosisyphus,
  combined=list(counts=counts.combined,summary=combined)),
  file.path(summaries.dir,"richness-through-time.rds"))

# 95% pointwise map-distribution interval, not a CI for the Monte Carlo mean ----
# Lines interpolate the 0.1 Ma grid; no smoothing of counts or uncertainty.
p <- ggplot(summary,aes(x=time)) +
  geom_ribbon(aes(ymin=lower,ymax=upper),fill="#0072B2",alpha=.18) +
  geom_line(aes(y=mean,linetype="Mean",colour="Mean"),linewidth=.7) +
  geom_line(aes(y=median,linetype="Median",colour="Median"),linewidth=.55) +
  scale_linetype_manual(NULL,values=c(Mean="solid",Median="dashed")) +
  scale_colour_manual(NULL,values=c(Mean="#0072B2",Median="black")) +
  scale_x_reverse(expand=expansion(mult=c(.01,.01))) +
  scale_y_continuous(expand=expansion(mult=c(0,.07))) +
  labs(x="Age (Ma)",y="Lineages on Sunken islands",
       caption="Shading: pointwise 95% interval; observed + extinct side lineages") +
  theme_classic(base_size=11) +
  theme(legend.position="top",plot.caption=element_text(size=8),
        strip.background=element_blank(),strip.text=element_text(face="italic"))
for (clade in c("Nesovinsonia","Nesosisyphus","Combined")) {
  plot <- p
  plot$data <- summary[summary$clade == clade,]
  plot <- plot + facet_wrap(~clade)
  filename <- paste0("richness-",tolower(clade))
  ggsave(file.path(plots.dir,paste0(filename,".pdf")),plot,width=6,height=4)
  ggsave(file.path(plots.dir,paste0(filename,".png")),plot,width=6,height=4,dpi=300)
}
plot <- p + facet_wrap(~clade,ncol=1,scales="free")
ggsave(file.path(plots.dir,"richness-all-panels.pdf"),plot,width=6.5,height=9)
ggsave(file.path(plots.dir,"richness-all-panels.png"),plot,width=6.5,height=9,dpi=300)

# Short numerical report ----
report <- c("# Sunken-island lineage richness", "",
 paste("Maps per clade:",length(nesovin)),
 "States pooled: A.S, S, S.R. Each living lineage counts once.",
 "Counts include observed backbone and extinct side lineages.",
 "Models are mixed using AIC weights; paired clades share model and posterior rates.",
 "Fixed ML-ASR root priors are used separately for each clade and model.",
 "This is a two-stage plug-in reconstruction, not a joint full-tree posterior.",
 "Intervals are 2.5% and 97.5% map quantiles (type 7), not CIs for the mean.",
 "The combined curve is restricted to ages represented in BOTH clades.", "")
for (clade in c("Nesovinsonia","Nesosisyphus","Combined")) {
  z <- summary[summary$clade == clade,]
  peak <- z[which.max(z$mean),]
  report <- c(report,sprintf("- %s: peak mean %.2f at %.2f Ma; median %.2f; pointwise 95%% interval %.2f–%.2f.",
    clade,peak$mean,peak$time,peak$median,peak$lower,peak$upper))
}
report <- c(report,"", "Peak-mean ages are grid-based; these are not distributions of each map's maximum.",
 "Full pointwise statistics are in richness-through-time.csv.")
writeLines(report,file.path(summaries.dir,"README.md"))
print(report)
