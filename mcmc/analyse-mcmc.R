# Run from the paleo-islands project root ----
# Install once if needed: install.packages("posterior")
library(posterior)

dir.create("mcmc/plots", showWarnings = FALSE)
dir.create("mcmc/summaries", showWarnings = FALSE)
summary.all <- data.frame()
runtime.all <- data.frame()

for (group in c(6, 7)) {
  # Read four chains; keep all post-warmup sweeps ----
  first <- readRDS(paste0("mcmc/results/roihu-5000-20261023/gr", group, "-chain-1.rds"))
  parameters <- colnames(first$rates)
  rates <- array(NA_real_, c(5000, 4, length(parameters)),
    dimnames = list(NULL, paste0("chain", 1:4), parameters))
  for (chain in 1:4) {
    x <- readRDS(paste0("mcmc/results/roihu-5000-20261023/gr", group, "-chain-", chain, ".rds"))
    # One check prevents incomplete or mismatched chains entering the analysis.
    stopifnot(nrow(x$rates) == 5000, x$warmup == 1000,
      identical(colnames(x$rates), parameters), x$group == group, x$chain == chain,
      all(is.finite(x$rates)), all(x$rates > 0))
    rates[, chain, ] <- x$rates
    runtime.all <- rbind(runtime.all, data.frame(group, chain,
      hours = unname(x$runtime["elapsed"])/3600))
  }
  draws <- rates[1001:5000, , ]

  # Summaries and rank-normalized convergence diagnostics ----
  d <- as.data.frame(summarise_draws(as_draws_array(draws),
    "mean", "median", "sd", "rhat", "ess_bulk", "ess_tail", "mcse_mean"))
  d$group <- group
  d$ML <- unname(first$ML[d$variable])
  for (p in seq_along(parameters)) {
    interval <- quantile(draws[, , p], c(.025, .975), type = 7)
    d$lower[p] <- interval[1]
    d$upper[p] <- interval[2]
    # Compare early and late retained halves; descriptive, not a formal test.
    d$median.first.half[p] <- median(draws[1:2000, , p])
    d$median.second.half[p] <- median(draws[2001:4000, , p])
    d$half.shift.over.sd[p] <- abs(d$median.first.half[p] - d$median.second.half[p])/d$sd[p]
  }
  d$mcse.over.sd <- d$mcse_mean/d$sd
  d$pass <- is.finite(d$rhat) & d$rhat < 1.01 &
    d$ess_bulk >= 400 & d$ess_tail >= 400 & d$mcse.over.sd < .05
  write.csv(d, paste0("mcmc/summaries/gr", group, "-parameters.csv"), row.names = FALSE)
  summary.all <- rbind(summary.all, d)
  print(d[, c("variable", "median", "lower", "upper", "rhat", "ess_bulk", "ess_tail", "pass")])

  # Posterior distributions and ML: vector PDF + 300 dpi PNG ----
  colour <- if (group == 6) "#0072B2" else "#D55E00"
  for (format in c("pdf", "png")) {
    filename <- paste0("mcmc/plots/gr", group, "-posterior.", format)
    if (format == "pdf") pdf(filename, width = 7.2, height = 10, useDingbats = FALSE)
    if (format == "png") png(filename, width = 7.2, height = 10, units = "in", res = 300)
    par(mfrow = c(5, 2), mar = c(3.4, 3.7, 2, .7), mgp = c(2, .6, 0),
      oma = c(3, 0, 0, 0), bty = "l", tcl = -.2, cex = .8)
    for (p in seq_along(parameters)) {
      values <- as.vector(draws[, , p])
      xmax <- max(quantile(values, .995), d$upper[p], d$ML[p])*1.05
      # Boundary-corrected KDE. Quantiles above use actual draws, not the KDE.
      density <- density(c(-values, values), bw = bw.nrd0(values), from = 0, to = xmax, n = 2048)
      density$y <- 2*density$y
      plot(density$x, density$y, type = "n", xlab = "Rate", ylab = "Density",
        main = parameters[p], ylim = c(-.07, 1.05)*max(density$y))
      polygon(c(density$x, rev(density$x)), c(density$y, rep(0, length(density$y))),
        col = adjustcolor(colour, alpha.f = .18), border = NA)
      inside <- density$x >= d$lower[p] & density$x <= d$upper[p]
      polygon(c(density$x[inside], rev(density$x[inside])),
        c(density$y[inside], rep(0, sum(inside))),
        col = adjustcolor(colour, alpha.f = .35), border = NA)
      lines(density$x, density$y, col = colour, lwd = 1.3)
      abline(v = d$ML[p], lty = 2, lwd = 1, col = "grey25")
      height <- -.035*max(density$y)
      segments(d$lower[p], height, d$upper[p], height, col = colour, lwd = 2)
      points(d$median[p], height, pch = 16, col = colour, cex = .65)
      points(d$mean[p], height, pch = 18, col = "black", cex = .65)
    }
    mtext("Shading/bar: 95% credible interval; circle: median; diamond: mean; dashed: ML",
      side = 1, outer = TRUE, line = .7, cex = .65)
    dev.off()
  }

  # Trace diagnostics: four chains, log rates, after warmup ----
  for (format in c("pdf", "png")) {
    filename <- paste0("mcmc/plots/gr", group, "-traces.", format)
    if (format == "pdf") pdf(filename, width = 7.2, height = 10, useDingbats = FALSE)
    if (format == "png") png(filename, width = 7.2, height = 10, units = "in", res = 300)
    par(mfrow = c(5, 2), mar = c(3.4, 3.7, 2, .7), mgp = c(2, .6, 0),
      oma = c(3, 0, 0, 0), bty = "l", tcl = -.2, cex = .8)
    for (p in seq_along(parameters)) {
      matplot(1001:5000, log(draws[, , p]), type = "l", lty = 1, lwd = .4,
        col = c("#0072B2", "#D55E00", "#009E73", "#CC79A7"),
        xlab = "Sweep", ylab = "Log rate", main = parameters[p])
    }
    mtext("Chains 1-4: blue, orange, green, purple; first 1,000 sweeps omitted",
      side = 1, outer = TRUE, line = .7, cex = .65)
    dev.off()
  }
}

# Combined tables and short report ----
write.csv(summary.all, "mcmc/summaries/all-parameters.csv", row.names = FALSE)
write.csv(runtime.all, "mcmc/summaries/runtime.csv", row.names = FALSE)
report <- c("FINAL EMPIRICAL MCMC: GROUPS 6 AND 7", "",
  "Four chains per model; 5,000 sweeps each, first 1,000 discarded.",
  "16,000 retained joint draws per model, without thinning.",
  "Independent exponential mean-1 priors on free rates; fixed log-rate transformation.", "")
for (group in c(6, 7)) {
  d <- summary.all[summary.all$group == group, ]
  report <- c(report, paste0("Group ", group),
    sprintf("Diagnostic screen: %d/%d parameters pass.", sum(d$pass), nrow(d)),
    sprintf("Maximum R-hat %.5f; minimum bulk ESS %.0f; minimum tail ESS %.0f.",
      max(d$rhat), min(d$ess_bulk), min(d$ess_tail)),
    sprintf("Largest mean MCSE/posterior SD: %.3f; largest half-to-half median shift/SD: %.3f.",
      max(d$mcse.over.sd), max(d$half.shift.over.sd)), "")
}
report <- c(report,
  "Pass: rank-normalized split R-hat <1.01, bulk/tail ESS >=400, mean MCSE/SD <0.05.",
  "Prefer ESS >=1,000 for final tail summaries. Inspect traces; thresholds do not prove convergence.",
  "Intervals are pointwise equal-tailed 95% posterior credible intervals, not ML confidence intervals.",
  "Density curves are zero-reflected KDEs; plots display approximately the first 99.5% of draws.",
  "The displayed curves are not renormalized after cropping. All draws enter summaries.",
  "Inference is conditional on each model, its priors, sampling assumptions and fixed tree.",
  "Near-zero ML outside a credible interval is not evidence against an exact-zero boundary model.")
writeLines(report, "mcmc/summaries/convergence-report.txt")
writeLines(capture.output(packageVersion("posterior")), "mcmc/summaries/posterior-version.txt")
