# Run from the paleo-islands project root, after analyse-mcmc.R ----
summary <- read.csv("mcmc/summaries/all-parameters.csv")
draws6 <- NULL
draws7 <- NULL
for (chain in 1:4) {
  x <- readRDS(paste0("mcmc/results/roihu-5000-20261023/gr6-chain-", chain, ".rds"))
  draws6 <- rbind(draws6, x$rates[1001:5000, ])
  x <- readRDS(paste0("mcmc/results/roihu-5000-20261023/gr7-chain-", chain, ".rds"))
  draws7 <- rbind(draws7, x$rates[1001:5000, ])
}

# Eight comparable rates, followed by three model-specific contractions ----
parameters <- c("lambda010101.1", "mu01.1", "q0108.1", "q0516.1",
  "lambda060606.1", "lambda080103.1", "lambda070102.1", "q0107.1",
  "mu02.1", "mu02.1", "mu06.1")
labels <- c("[sA]", "[xA]", "[dAM, dEM]", "[dSR, dAS]", "[sR, sS]",
  "[scAM, scEM, scSR, scAS]", "[scA, scE, sE, sM, sU]", "[dAE, dAU, dEA, dEU]",
  "Group 6\n[xE, xM, xU, xR, xS]", "Group 7\n[xE, xM, xU]", "Group 7\n[xR, xS]")
colours <- c("#0072B2", "#D55E00")
# Identical x limits for the three model-specific contraction panels.
contraction.limit <- max(quantile(draws6[, "mu02.1"], .995),
  quantile(draws7[, "mu02.1"], .995), quantile(draws7[, "mu06.1"], .995),
  summary$ML[summary$variable %in% c("mu02.1", "mu06.1")])*1.05

dir.create("mcmc/plots", showWarnings = FALSE)
for (format in c("pdf", "png")) {
  file <- paste0("mcmc/plots/gr6-gr7-posterior.", format)
  if (format == "pdf") pdf(file, width = 10.5, height = 10, useDingbats = FALSE)
  if (format == "png") png(file, width = 10.5, height = 10, units = "in", res = 300)
  par(mfrow = c(4, 3), mar = c(3.4, 3.7, 3, .7), mgp = c(2, .6, 0),
    oma = c(2.2, 0, 0, 0), bty = "l", tcl = -.2, cex = .8)
  for (p in seq_along(parameters)) {
    models <- c(6, 7)
    if (p == 9) models <- 6
    if (p >= 10) models <- 7
    values <- list()
    pooled <- numeric()
    xmax <- 0
    for (model in models) {
      v <- if (model == 6) draws6[, parameters[p]] else draws7[, parameters[p]]
      values[[as.character(model)]] <- v
      pooled <- c(pooled, v)
      s <- summary[summary$group == model & summary$variable == parameters[p], ]
      xmax <- max(xmax, quantile(v, .995), s$ML)
    }
    xmax <- xmax*1.05
    if (p >= 9) xmax <- contraction.limit
    bandwidth <- bw.nrd0(pooled)
    curves <- list()
    ymax <- 0
    for (model in models) {
      v <- values[[as.character(model)]]
      d <- density(c(-v, v), bw = bandwidth, from = 0, to = xmax, n = 2048)
      d$y <- 2*d$y
      curves[[as.character(model)]] <- d
      ymax <- max(ymax, d$y)
    }
    plot(NA, xlim = c(0, xmax), ylim = c(-.22, 1.42)*ymax,
      main = labels[p], xlab = "Rate", ylab = "Density", yaxt = "n", cex.main = .95)
    ticks <- pretty(c(0, ymax))
    axis(2, at = ticks[ticks >= 0 & ticks <= ymax])
    for (model in models) {
      d <- curves[[as.character(model)]]
      colour <- colours[model - 5]
      polygon(c(d$x, rev(d$x)), c(d$y, rep(0, length(d$y))),
        col = adjustcolor(colour, alpha.f = .18), border = NA)
      lines(d$x, d$y, col = colour, lwd = 1.3)
      s <- summary[summary$group == model & summary$variable == parameters[p], ]
      height <- -.07*match(model, models)*ymax
      segments(s$lower, height, s$upper, height, col = colour, lwd = 2)
      points(s$median, height, pch = 16, col = colour, cex = .8)
      points(s$ML, height, pch = 4, col = colour, cex = .8, lwd = 1.2)
      ml.label <- paste0("G", model, " ML = ", formatC(s$ML, format = "g", digits = 3))
      text(.98*xmax, (1.34 - .16*(match(model, models)-1))*ymax,
        labels = ml.label, adj = c(1, .5), col = colour, cex = .85)
    }
  }
  plot.new()
  legend("center", legend = c("Group 6", "Group 7"), fill = colours,
    border = NA, bty = "n", cex = 1.1, y.intersp = 1.5)
  # CI here denotes equal-tailed posterior credible intervals.
  mtext("Point: median; Cross: ML estimate;  Bar: 95% CI",
    side = 1, outer = TRUE, line = .4, cex = .8)
  dev.off()
}
