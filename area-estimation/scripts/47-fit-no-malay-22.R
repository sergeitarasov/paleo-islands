# Run from the paleo-islands project root.
library(MASS)
library(posterior)
library(ggplot2)
library(readxl)
dir.create("area-estimation/results/bayesian-no-malay-22-10ma", recursive = TRUE, showWarnings = FALSE)
dir.create("area-estimation/plots/bayesian-no-malay-22-10ma", recursive = TRUE, showWarnings = FALSE)
dir.create("area-estimation/logs/bayesian-no-malay-22-10ma", recursive = TRUE, showWarnings = FALSE)
# Compare pooled and region-specific relationships on the same 22 islands after the stated exclusions.
d <- read.csv("area-estimation/results/selected-27/input.csv")
d <- d[!d$island %in% c("Borneo", "Java", "Sumatra", "Sulawesi", "New Guinea"), ]
stopifnot(nrow(d) == 22)
d$x <- log(d$area_km2 / 10000)
write.csv(d, "area-estimation/results/bayesian-no-malay-22-10ma/calibration-data.csv", row.names = FALSE)
d$group <- factor(ifelse(d$region == "Caribbean", "Caribbean", "Other"), levels = c("Caribbean", "Other"))
write.csv(d, "area-estimation/results/bayesian-no-malay-22-10ma/calibration-data.csv", row.names = FALSE)
formulas <- list(simple = ~ x, separate = ~ x * group)
results <- list()
diagnostics <- data.frame()
acceptance <- data.frame()
traces <- data.frame()
started <- proc.time()[3]

# Four chains; proposal covariance fixed at posterior-mode curvature.
# Only proposal scale is tuned during warmup, then held fixed.
for (model in names(formulas)) {
  X <- model.matrix(formulas[[model]], d)
  p <- ncol(X)
  for (prior_width in 1) {
    id <- paste0(model, "_prior", prior_width)
    prior_mean <- rep(0, p)
    prior_mean[1] <- log(10)
    prior_mean[2] <- .4
    prior_sd <- rep(1.5, p)
    prior_sd[1] <- 2
    prior_sd[grepl("x", colnames(X), fixed = TRUE)] <- .5
    prior_sd <- prior_sd * prior_width
    # log(sigma) or log(NB size) is sampled with the appropriate prior.
    log_target <- function(theta) {
      beta <- theta[seq_len(p)]
      eta <- as.vector(X %*% beta)
      scale <- exp(theta[p+1])
      if (!is.finite(scale) || scale < 1e-10 || scale > 1e10) return(-Inf)
      lp <- sum(dnorm(beta, prior_mean, prior_sd, log = TRUE))
      if (model %in% c("simple", "separate")) {
        ll <- sum(dnorm(log(d$n_sp), eta, scale, log = TRUE))
        lp <- lp + dnorm(scale, 0, 1.5*prior_width, log = TRUE) + log(scale)
      } else {
        ll <- sum(dnbinom(d$n_sp, mu = exp(eta), size = scale, log = TRUE))
        lp <- lp + dnorm(log(scale), log(2), 1.5*prior_width, log = TRUE)
      }
      if (!is.finite(ll + lp)) return(-Inf)
      ll + lp
    }
    start <- c(prior_mean, 0)
    mode <- optim(start, function(z) -log_target(z), method = "BFGS", hessian = TRUE,
                  control = list(maxit = 2000, reltol = 1e-10))
    stopifnot(mode$convergence == 0, min(eigen(mode$hessian, symmetric = TRUE)$values) > 0)
    L <- t(chol(solve(mode$hessian)))
    chains <- array(NA_real_, c(25000, 4, p+1),
                    dimnames = list(NULL, paste0("chain",1:4), c(colnames(X), "log_scale")))
    for (chain in 1:4) {
      set.seed(261000 + 100*match(model, names(formulas)) + 10*prior_width + chain)
      theta <- as.vector(mode$par + L %*% rnorm(p+1, sd = 2))
      value <- log_target(theta)
      step <- 2.38/sqrt(p+1)
      accepted <- 0
      block <- 0
      for (iteration in 1:30000) {
        proposal <- as.vector(theta + step * L %*% rnorm(p+1))
        candidate <- log_target(proposal)
        if (log(runif(1)) < candidate - value) {
          theta <- proposal; value <- candidate
          block <- block + 1
          if (iteration > 5000) accepted <- accepted + 1
        }
        if (iteration <= 5000 && iteration %% 100 == 0) {
          step <- step * exp(.3*(block/100 - .234))
          block <- 0
        }
        if (iteration > 5000) chains[iteration-5000, chain, ] <- theta
      }
      acceptance <- rbind(acceptance, data.frame(model = id, chain = chain, acceptance = accepted/25000))
    }
    summary <- as.data.frame(summarise_draws(as_draws_array(chains)))
    summary$model <- id
    diagnostics <- rbind(diagnostics, summary)
    draws <- rbind(chains[,1,], chains[,2,], chains[,3,], chains[,4,])
    results[[id]] <- list(draws = draws, chains = chains, formula = formulas[[model]],
      prior_mean = prior_mean, prior_sd = prior_sd, mode = mode$par)
    for (j in seq_len(p+1)) {
      for (chain in 1:4) {
        traces <- rbind(traces, data.frame(model = id, parameter = dimnames(chains)[[3]][j],
          chain = factor(chain), iteration = seq(1,25000,by=50), value = chains[seq(1,25000,by=50),chain,j]))
      }
    }
    cat(id, "max Rhat", max(summary$rhat), "min ESS", min(summary$ess_bulk, summary$ess_tail), "\n")
    saveRDS(results, "area-estimation/results/bayesian-no-malay-22-10ma/regression-posteriors.rds")
  }
}
write.csv(diagnostics, "area-estimation/results/bayesian-no-malay-22-10ma/convergence.csv", row.names = FALSE)
write.csv(acceptance, "area-estimation/results/bayesian-no-malay-22-10ma/acceptance.csv", row.names = FALSE)
for (model in names(results)) {
  z <- traces[traces$model == model, ]
  plot <- ggplot(z, aes(iteration, value, colour = chain)) + geom_line(linewidth = .2, alpha = .65) +
    facet_wrap(~parameter, scales = "free_y", ncol = 3) + theme_bw(base_size = 9) + labs(y = NULL)
  ggsave(paste0("area-estimation/plots/bayesian-no-malay-22-10ma/traces-", model, ".pdf"), plot, width = 10, height = 7)
}
cat("Elapsed seconds:", proc.time()[3]-started, "\n")
capture.output(sessionInfo(), file = "area-estimation/logs/bayesian-no-malay-22-10ma/session-info.txt")
stopifnot(all(diagnostics$rhat < 1.01), all(diagnostics$ess_bulk > 400), all(diagnostics$ess_tail > 400))
