# Generated-script template for the Oriental16 16-state hypothesis.
library(ape)
library(phytools)
library(readxl)
library(geiger)
library(plyr)
library(dplyr)
library(phangorn)
library(tidyverse)
library(jsonlite)
library(rstudioapi)
library(fs)
library(divextra)

# Values inside double braces are replaced by make-scripts.R.
PARAM_GR <- 6
RUN <- 5
BASE_NAME <- "Oriental16"
YAML_BACKBONE <- "../../yml/backbone-Oriental-16.yml"
YAML_GROUPS <- "../../yml/groups-Oriental-16.yml"
K <- 16
N_EPOCH <- 2
EPSs <- seq(0.1, 0.8, length.out = 5)

EPS <- EPSs[RUN]
OUTPUT_BASE <- paste0(BASE_NAME, "_gr-", PARAM_GR)
print(OUTPUT_BASE)
print(RUN)

# Data and observed tip states.
phy <- readRDS("../../data/phy-271.rds")
geo_regions_abc <- readRDS("../../data/geo_regions-271.rds")
geo_regions <- mapvalues(
  geo_regions_abc,
  from = c("Afr", "OP", "Mada", "Aus", "Maur"),
  to = c(1:4, 6)
) %>% as.numeric()
names(geo_regions) <- names(geo_regions_abc)

# State order:
# A, E, M, U, S, R, A.E, A.M, A.U, A.S, A.R, E.M, E.U, E.S, E.R, S.R
single_transitions.f <- c(0.039, 0.027, 0.446, 0.008, 1e-06, 0.833)
dual_transitions.f <- rep(1e-06, 10)
sampling.f <- c(single_transitions.f, dual_transitions.f)
names(sampling.f) <- seq_len(K)

# Read the backbone and select one equality-constraint model.
par.categories.td <- read_yaml_pars_td(YAML_BACKBONE)
par_groups_yaml <- read_par_groups_yaml(YAML_GROUPS)
if (PARAM_GR < 1 || PARAM_GR > length(par_groups_yaml)) {
  stop("PARAM_GR is outside the available grouping range")
}
focal_group <- par_groups_yaml[[PARAM_GR]]
new_pars <- regroup_parameters(par.categories.td, focal_group)

create_parameter_html(
  par.categories.td = new_pars,
  model_name = OUTPUT_BASE,
  output_file = path("../", "html", paste0(OUTPUT_BASE, ".html"))
)

# Maximum-likelihood fit.
lik.td <- make.classe.td(
  phy,
  geo_regions,
  k = K,
  n.epoch = N_EPOCH,
  control = list(backend = "gslode"),
  strict = FALSE,
  sampling.f = sampling.f
)
formula.td <- make_constraints_sse_td(new_pars)
lik.const.td <- constrain(lik.td, formulae = formula.td)
starting.point <- init.pars.classe_td(
  lik.const.td,
  phy,
  k = K,
  n.epoch = N_EPOCH,
  eps = EPS
)
print(starting.point)

mle.td <- tryCatch(
  find.mle(
    lik.const.td,
    starting.point,
    condition.surv = TRUE,
    keep.func = FALSE,
    root = ROOT.GIVEN,
    root.p = c(1, rep(0, K - 1)),
    control = list(maxit = 20000)
  ),
  error = function(e) {
    message(sprintf("Failed at EPS=%s: %s", EPS, conditionMessage(e)))
    NULL
  }
)
if (is.null(mle.td)) {
  quit(save = "no", status = 1)
}

# Validate and save the fitted parameterization.
mle.pars.groups <- pars2groups(mle.td$par.full)
check_parametrization_cons <- check_parametrization_consistency(
  mle.pars.groups,
  new_pars$pars
)

mle.td <- append(list(AIC = AIC(mle.td)), mle.td)
mle.td <- append(list(NPAR = length(mle.td$par)), mle.td)
mle.td <- append(list(PARAM_CONS = check_parametrization_cons), mle.td)
mle.td <- append(list(EPS = EPS), mle.td)
mle.td <- append(list(RUN = RUN), mle.td)
mle.td <- append(list(PARAM_GR = PARAM_GR), mle.td)
mle.td <- append(list(BASE_NAME = BASE_NAME), mle.td)
mle.td$par_json <- as.list(mle.td$par)
mle.td$par.full_json <- as.list(mle.td$par.full)

write_json(
  mle.td,
  path = path("../", "json", paste0(OUTPUT_BASE, "_r-", RUN, ".json")),
  pretty = TRUE,
  auto_unbox = FALSE,
  null = "null",
  digits = 6,
  force = TRUE
)
saveRDS(
  mle.td,
  file = path("../", "rds", paste0(OUTPUT_BASE, "_r-", RUN, ".rds"))
)


