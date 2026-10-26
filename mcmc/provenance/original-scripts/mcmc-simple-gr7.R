# Run from the divextra-sunken project root ----
library(divextra)
library(diversitree)

# Data ----
phy <- readRDS("data/phy-271.rds")
geo.regions <- readRDS("data/geo_regions-271.rds")
states <- c(Afr = 1L, OP = 2L, Mada = 3L, Aus = 4L, Maur = 6L)
states <- states[as.character(geo.regions)]
names(states) <- names(geo.regions)

# Select the model ----
group <- 7
fit <- readRDS("runs/rds/AfrSunken16_gr-7_r-1.rds")
transform <- readRDS("asr/bayesian-gr7/tuned/transform.rds")

# Likelihood: same model as the saved ML analysis ----
pars <- read_yaml_pars_td("yml/backbone-Afr-Sunken-16.yml")
groups <- read_par_groups_yaml("yml/groups-Afr-Sunken-16.yml")
pars <- regroup_parameters(pars, groups[[group]])

lik <- make.classe.td(
  phy, states, k = 16, n.epoch = 2, strict = FALSE,
  sampling.f = c(.039, .027, .446, .008, 1e-6, .833, rep(1e-6, 10)),
  control = list(backend = "gslode")
)
lik <- constrain(lik, formulae = make_constraints_sse_td(pars))

# Posterior ----
# u is a transformed log-rate vector. The fixed pilot transformation
# helps sampling correlated rates; it does not change the posterior.
logposterior <- function(u) {
  log.rates <- as.vector(transform$center + transform$L %*% u)
  rates <- exp(log.rates)
  names(rates) <- names(fit$par)
  if (any(!is.finite(rates)) || any(rates <= 0)) return(-Inf)

  lnL <- lik(rates, root = ROOT.GIVEN,
    root.p = c(1, rep(0, 15)), condition.surv = TRUE)
  prior <- sum(dexp(rates, rate = 1, log = TRUE))
  # Jacobian: exponential priors are on rates, not log rates.
  lnL + prior + sum(log.rates)
}

# Four fresh chains ----
# Each chain has 5,000 sweeps; discard the first 1,000 in analysis.
# By default run all four sequentially; MCMC_CHAIN selects one on Roihu.
# Each chain is saved immediately after completion.
# Rerunning this section replaces files with the same names.
chains <- as.integer(Sys.getenv("MCMC_CHAIN", "0"))
if (chains == 0) chains <- 1:4
output <- Sys.getenv("MCMC_OUTPUT", "asr/mcmc-results")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
for (chain in chains) {
  set.seed(20261027 + group * 100 + chain)
  start <- rnorm(length(fit$par), sd = 1.5)
  names(start) <- paste0("u", seq_along(start))

  runtime <- system.time(
    samples <- mcmc(
      logposterior, x.init = start, nsteps = 5000,
      w = rep(1, length(fit$par)), print.every = 100,
      keep.func = FALSE, fail.value = NULL
    )
  )

  # Convert sampler coordinates back into biological rates.
  log.rates <- as.matrix(coef(samples)) %*% t(transform$L)
  log.rates <- sweep(log.rates, 2, transform$center, "+")
  rates <- exp(log.rates)
  colnames(rates) <- names(fit$par)

  saveRDS(
    list(rates = rates, warmup = 1000, samples = samples, runtime = runtime,
      group = group, chain = chain, prior.rate = 1, transform = transform,
      ML = fit$par, seed = 20261027 + group * 100 + chain,
      session = sessionInfo()),
    file.path(output, paste0("gr", group, "-chain-", chain, ".rds"))
  )
}
