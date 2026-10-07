##
## Tuberculosis: Household, Regular, and Casual Contacts Over a Multilayer Network
## EpiModel Gallery (https://github.com/EpiModel/EpiModel-Gallery)
##
## Author: Samuel M. Jenness (Emory University)
## Date: September 2026
##

# Load EpiModel
suppressMessages(library(EpiModel))

# Standard Gallery unit test lines
rm(list = ls())
eval(parse(text = print(commandArgs(TRUE)[1])))

# Run settings. The full settings are used when the script is run
# interactively (for example sourced in RStudio). Rscript is not interactive,
# so it uses the small CI settings unless the first command-line argument
# defines run_full, which the unit test line above evaluates as R code:
#   Rscript examples/tuberculosis/model.R "run_full <- TRUE"
# The full run takes about half an hour on five cores; CI mode runs in about
# a minute and its results are not meant to be interpreted.
if (interactive() || exists("run_full")) {
  N <- 20000
  nsims <- 10
  ncores <- 5
  burnin <- 600
  horizon <- 120
} else {
  N <- 1000
  nsims <- 1
  ncores <- 1
  burnin <- 24
  horizon <- 12
}

# model.R is written to be run from the repository root; the fallback lets
# the downloaded script run from a directory that holds both files.
source(if (file.exists("examples/tuberculosis/module-fx.R")) {
  "examples/tuberculosis/module-fx.R"
} else {
  "module-fx.R"
})


# 1. Population and Households ----------------------------------------------

# Background mortality by age band (annual rates starting at the ages in
# mort.breaks). A stationary population with these rates has a life
# expectancy at birth of about 60 years and about a quarter of its members
# under 15; births at a per-capita rate of 1 / e0 keep its size stable.
mort.breaks <- c(0, 1, 5, 15, 30, 45, 60, 70, 80)
mort.rates <- c(0.025, 0.002, 0.0012, 0.004, 0.009, 0.020, 0.040, 0.080, 0.180)
age_grid <- seq(0, 100 - 1 / 12, by = 1 / 12)
mort_grid <- mort.rates[findInterval(age_grid, mort.breaks)]
surv <- exp(-cumsum(c(0, head(mort_grid, -1) / 12)))
e0 <- sum(surv) / 12

# Ages are drawn from the stationary distribution, then people are grouped
# into households of a target size distribution with every child placed in a
# household that has an adult (assign_groups).
set.seed(12)
age <- sample(age_grid, N, replace = TRUE, prob = surv) + runif(N, 0, 1 / 12)
agegrp <- ifelse(age < 15, "child", "adult")
hh_size_dist <- c("1" = 0.22, "2" = 0.18, "3" = 0.17, "4" = 0.16, "5" = 0.11,
                  "6" = 0.08, "7" = 0.05, "8" = 0.03)
hh_id <- assign_groups(hh_size_dist, role = agegrp, anchor = "adult",
                       dependent = "child")

cat(sprintf("N = %d in %d households; life expectancy %.1f years; %.1f%% under 15\n",
            N, max(hh_id), e0, 100 * mean(age < 15)))


# 2. Contact Layers ---------------------------------------------------------

nw <- network_initialize(N)
nw <- set_vertex_attribute(nw, "age", age)
nw <- set_vertex_attribute(nw, "agegrp", agegrp)
nw <- set_vertex_attribute(nw, "hh_id", hh_id)

# Layer 1, household: every pair of co-residents is an edge for as long as
# they live together. Newborns join the household of a random adult aged 15
# to 49 (newborn_household in module-fx.R).
est_hh <- netclique(nw, group.attr = "hh_id", arrivals = "join",
                    arrivals.FUN = newborn_household)
print(est_hh, by = "agegrp")

# Layers 2 and 3, regular and casual non-household contacts: TERGMs with a
# fully targeted age-mixing matrix. Each profile entry "a.b" is the mean
# number of b-partners per a-node.
counts <- table(agegrp)
cells <- sub("^mix\\.agegrp\\.", "",
             names(summary(nw ~ nodemix("agegrp", levels2 = TRUE))))
mix_targets <- function(profile) {
  out <- setNames(numeric(length(cells)), cells)
  for (nm in names(profile)) {
    ab <- strsplit(nm, ".", fixed = TRUE)[[1]]
    n_a <- as.numeric(counts[ab[1]])
    e <- if (ab[1] == ab[2]) n_a * profile[[nm]] / 2 else n_a * profile[[nm]]
    cell <- paste(sort(ab), collapse = ".")
    out[cell] <- out[cell] + e
  }
  pmax(round(out), 1)
}
reg_profile <- c(adult.adult = 4.6, child.adult = 1.0, child.child = 5.0)
cas_profile <- c(adult.adult = 5.5, child.adult = 1.5, child.child = 1.5)
reg_targets <- mix_targets(reg_profile)
cas_targets <- mix_targets(cas_profile)

formation <- ~edges + nodemix("agegrp", levels2 = -1)
san_ctrl <- control.ergm(SAN = control.san(SAN.maxit = 20, SAN.nsteps = 2^21))
d.rate <- 1 / e0 / 12

est_reg <- netest(nw, formation,
                  target.stats = c(sum(reg_targets), as.numeric(reg_targets[-1])),
                  coef.diss = dissolution_coefs(~offset(edges), duration = 12,
                                                d.rate = d.rate),
                  set.control.ergm = san_ctrl, verbose = FALSE)
est_cas <- netest(nw, formation,
                  target.stats = c(sum(cas_targets), as.numeric(cas_targets[-1])),
                  coef.diss = dissolution_coefs(~offset(edges), duration = 1),
                  set.control.ergm = san_ctrl, verbose = FALSE)

# Per-step Markov chain lengths for simulating the TERGM layers. tergm's
# default chain is far too short for these layers: a layer whose ties last
# one step must be redrawn completely every step (with the default chain
# about half of the casual ties carried over from one month to the next),
# and the regular layer's dense mixing among children fell about 20% below
# its target. Both need a chain that scales with the number of edges, hence
# with N.
mcmc <- function(steps) {
  control.simulate.formula.tergm(MCMC.burnin.min = steps, MCMC.burnin.max = steps)
}
tergm_persist <- mcmc(20 * N)
tergm_redraw <- mcmc(25 * N)

# Diagnostics: realized mean degree by age group from simulations of the
# fitted models. The regular layer is simulated dynamically, which also
# checks its mean tie duration; the casual layer's ties last one step, so
# its cross-sectional model is checked with static simulations.
dx_reg <- netdx(est_reg, nsims = 2, nsteps = 60, set.control.tergm = tergm_persist,
                nwstats.formula = ~edges + nodefactor("agegrp", levels = TRUE),
                verbose = FALSE)
dx_cas <- netdx(est_cas, nsims = 10, dynamic = FALSE,
                nwstats.formula = ~edges + nodefactor("agegrp", levels = TRUE),
                verbose = FALSE)
print(dx_reg$stats.table.duration)
realized_degree <- function(dx) {
  st <- colMeans(do.call(rbind, lapply(dx$stats, as.matrix)))
  st[paste0("nodefactor.agegrp.", c("adult", "child"))] /
    as.numeric(counts[c("adult", "child")])
}

# Mean degree by age group and layer: the household layer's is read from the
# cliques, and the community layers' targets are set beside the values from
# the diagnostic simulations
hh_el <- as.edgelist(est_hh$newnetwork)
hh_deg <- tabulate(c(hh_el), nbins = N)
deg_tbl <- rbind(
  household = tapply(hh_deg, agegrp, mean),
  regular = c(adult = 2 * reg_targets[["adult.adult"]] + reg_targets[["adult.child"]],
              child = reg_targets[["adult.child"]] + 2 * reg_targets[["child.child"]]) /
    as.numeric(counts[c("adult", "child")]),
  casual = c(adult = 2 * cas_targets[["adult.adult"]] + cas_targets[["adult.child"]],
             child = cas_targets[["adult.child"]] + 2 * cas_targets[["child.child"]]) /
    as.numeric(counts[c("adult", "child")])
)
deg_show <- cbind(deg_tbl, rbind(household = deg_tbl["household", ],
                                 regular = realized_degree(dx_reg),
                                 casual = realized_degree(dx_cas)))
colnames(deg_show) <- c("adult, target", "child, target",
                        "adult, simulated", "child, simulated")
print(round(deg_show, 2))

# Contact time. McCreesh and White (2018) report, for adults in the Western
# Cape, about 17.0 person-hours a day with household members, 6.0 with
# repeated non-household contacts, and 24.2 with non-repeated contacts. The
# per-edge monthly hours below divide each budget by the adult mean degree
# on the layer, so an adult with the average number of contacts has the
# reported daily contact time on each layer.
hours_day <- c(household = 17.0, regular = 6.0, casual = 24.2)
layer_hours <- hours_day / deg_tbl[, "adult"] * 365.25 / 12
print(round(layer_hours, 1))


# 3. Parameters ---------------------------------------------------------------

init <- init.net(status.vector = rep("s", N))

param_base <- param.net(
  # transmission
  beta = 0.0046,
  layer.hours = unname(layer_hours),
  layer.names = names(layer_hours),
  inf.age.breaks = c(0, 10, 15),
  inf.age.mult = c(0, 0.5, 1),
  inf.k = 0.15,
  sus.latent = 0.21,
  sus.recovered = 1,
  # natural history, annual rates for ages 0-4, 5-14, 15+
  prog.early = c(0.773, 0.231, 0.219),
  stab.early = c(4.38, 4.38, 1.97),
  react.late = c(0, 0.00234, 0.00121),
  clear.late = 0.028,
  diag.rate = 1.0,
  tb.mort = c(0.12, 0.03, 0.10),
  self.cure = 0.15,
  tx.duration = 6,
  tx.mort = 0.0059,
  tx.success = 0.93,
  relapse.rate = 0.024,
  rec.stab = 0.5,
  # programs, all off in the base parameter set
  hhci = 0,
  hhci.reach = 0.65,
  hhci.eval = 0.85,
  screen.sens = c(0.7, 0.8),
  tpt.age.max = Inf,
  tpt.start = 0.85,
  tpt.complete = 0.8,
  tpt.eff = 0.9,
  acf = 0,
  acf.interval = 12,
  acf.coverage = 0.42,
  acf.sens = 0.85,
  acf.min.age = 15,
  random.rate = 0,
  # demography
  birth.rate = 1 / e0,
  mort.rates = mort.rates,
  mort.breaks = mort.breaks,
  leave.rate = 0.0015,
  leave.ages = c(18, 35),
  # initial conditions
  init.arti = 0.025,
  init.early.frac = 0.05,
  init.prev.adult = 0.004,
  init.rec.adult = 0.02
)

make_control <- function(nsteps, tergm, start = 1, save.run = FALSE) {
  control.net(
    type = NULL, nsims = nsims, ncores = ncores, nsteps = nsteps,
    start = start, save.run = save.run,
    tergmLite = TRUE, resimulate.network = TRUE,
    initTB.FUN = init_tb, aging.FUN = aging, infection.FUN = infect,
    progress.FUN = progress, screen.FUN = screen, departures.FUN = deaths,
    arrivals.FUN = births, households.FUN = households, tally.FUN = tally,
    module.order = c("resim_nets.FUN", "summary_nets.FUN", "initTB.FUN",
                     "aging.FUN", "infection.FUN", "progress.FUN",
                     "screen.FUN", "departures.FUN", "arrivals.FUN",
                     "nwupdate.FUN", "households.FUN", "tally.FUN",
                     "prevalence.FUN"),
    set.control.tergm = tergm, save.nwstats = FALSE, verbose = FALSE
  )
}


# 4. Part 1: Three Contact Structures ---------------------------------------

# The same people, the same contact time on every layer, and the same
# per-contact-hour risk, arranged three ways:
#   baseline:          household cliques, regular ties lasting a year on
#                      average, casual ties redrawn every month
#   regular_redrawn:   regular ties redrawn every month (same mean degree)
#   household_redrawn: household contact time moved to a layer of random
#                      ties with the household layer's mean degree and age
#                      mixing, redrawn every month; the clique layer is kept
#                      for bookkeeping but carries no transmission
# Each structure runs its own burn-in from the same initial conditions.
est_reg1 <- netest(nw, formation,
                   target.stats = c(sum(reg_targets), as.numeric(reg_targets[-1])),
                   coef.diss = dissolution_coefs(~offset(edges), duration = 1),
                   set.control.ergm = san_ctrl, verbose = FALSE)
hh_pair <- paste(pmin(agegrp[hh_el[, 1]], agegrp[hh_el[, 2]]),
                 pmax(agegrp[hh_el[, 1]], agegrp[hh_el[, 2]]), sep = ".")
hh_targets <- table(factor(hh_pair, levels = cells))
est_hhr <- netest(nw, formation,
                  target.stats = c(sum(hh_targets), as.numeric(hh_targets[-1])),
                  coef.diss = dissolution_coefs(~offset(edges), duration = 1),
                  set.control.ergm = san_ctrl, verbose = FALSE)

tergm_none <- control.simulate.formula.tergm()
structures <- list(
  baseline = list(
    nets = list(est_hh, est_reg, est_cas),
    hours = layer_hours,
    tergm = multilayer(tergm_none, tergm_persist, tergm_redraw)),
  regular_redrawn = list(
    nets = list(est_hh, est_reg1, est_cas),
    hours = layer_hours,
    tergm = multilayer(tergm_none, tergm_redraw, tergm_redraw)),
  household_redrawn = list(
    nets = list(est_hh, est_reg, est_cas, est_hhr),
    hours = c(household = 0, layer_hours[c("regular", "casual")],
              household_redrawn = layer_hours[["household"]]),
    tergm = multilayer(tergm_none, tergm_persist, tergm_redraw, tergm_redraw))
)

burn <- list()
for (s in names(structures)) {
  cat("Burn-in:", s, "\n")
  st <- structures[[s]]
  param <- param_base
  param$layer.hours <- unname(st$hours)
  param$layer.names <- names(st$hours)
  set.seed(101)
  burn[[s]] <- netsim(st$nets, param, init,
                      make_control(burnin, st$tergm, save.run = (s == "baseline")))
}


# 5. Analysis Helpers and Part 1 Results -------------------------------------

# Every summary is computed per simulation, so that Monte Carlo intervals can
# be built from the between-simulation variation. Rows of sim$epi are time
# steps; a resumed simulation carries the burn-in history in its first
# `burnin` rows. Rows are selected by position.
epi_sum <- function(sim, vars, rows) {
  colSums(Reduce(`+`, lapply(vars, function(v) as.matrix(sim$epi[[v]])[rows, , drop = FALSE])),
          na.rm = TRUE)
}
epi_mean <- function(sim, vars, rows) epi_sum(sim, vars, rows) / length(rows)
srcs <- c("recent", "reinf", "react", "relapse")
inc_vars <- function(grp) paste0("inc.", srcs, ".", grp)
per100k_yr <- function(events, pop, rows) 1e5 * events / pop / (length(rows) / 12)

# Mean with a 95% Monte Carlo interval across simulations (t quantile), and a
# formatter; with one simulation the interval is not defined.
mc <- function(x) {
  m <- mean(x)
  if (length(x) < 2) return(c(est = m, lo = NA, hi = NA))
  h <- qt(0.975, length(x) - 1) * sd(x) / sqrt(length(x))
  c(est = m, lo = m - h, hi = m + h)
}
fmt <- function(v, digits = 0) {
  if (is.na(v["lo"])) return(formatC(v["est"], format = "f", digits = digits))
  sprintf("%s (%s, %s)", formatC(v["est"], format = "f", digits = digits),
          formatC(v["lo"], format = "f", digits = digits),
          formatC(v["hi"], format = "f", digits = digits))
}

## 5a. The baseline against published targets (last 20 years of burn-in) ----

rows_b <- max(2, burnin - 239):burnin
calib <- function(sim, rows) {
  pop_c <- epi_mean(sim, "num.child", rows)
  pop_a <- epi_mean(sim, "num.adult", rows)
  inc_c <- epi_sum(sim, inc_vars("child"), rows)
  inc_a <- epi_sum(sim, inc_vars("adult"), rows)
  inf <- function(l, g) epi_sum(sim, paste0("inf.", l, ".", g), rows)
  layers <- sim$param$layer.names
  hh_share <- function(g) {
    inf(layers[1], g) / Reduce(`+`, lapply(layers, inf, g = g))
  }
  bands <- c("u5", "5to14", "adult")
  hhc <- function(v, b = bands) epi_sum(sim, paste0("hhc.", v, ".", b), rows)
  data.frame(
    inc = per100k_yr(inc_c + inc_a, pop_c + pop_a, rows),
    inc_child = per100k_yr(inc_c, pop_c, rows),
    inc_adult = per100k_yr(inc_a, pop_a, rows),
    child_share = 100 * inc_c / (inc_c + inc_a),
    prev = 1e5 * epi_mean(sim, c("i.num.child", "i.num.adult"), rows) / (pop_c + pop_a),
    deaths = per100k_yr(epi_sum(sim, c("tbdeath.child", "tbdeath.adult"), rows),
                        pop_c + pop_a, rows),
    cfr = 100 * epi_sum(sim, c("tbdeath.child", "tbdeath.adult"), rows) /
      (inc_c + inc_a),
    detected = 100 * epi_sum(sim, c("dx.passive.child", "dx.passive.adult"), rows) /
      (inc_c + inc_a),
    recent = 100 * epi_sum(sim, c(inc_vars("child")[1:2], inc_vars("adult")[1:2]), rows) /
      (inc_c + inc_a),
    hh_adult = 100 * hh_share("adult"),
    hh_child = 100 * hh_share("child"),
    hhc_active = 100 * hhc("active") / hhc("n"),
    hhc_infected = 100 * hhc("infected") / hhc("n"),
    hhc_active_u5 = 100 * hhc("active", "u5") / hhc("n", "u5"),
    hhc_infected_u5 = 100 * hhc("infected", "u5") / hhc("n", "u5"),
    hhc_risk = 100 * hhc("inc") / hhc("n"),
    hhc_tb_u5 = 100 * (hhc("active", "u5") + hhc("inc", "u5")) / hhc("n", "u5"),
    infected_15to34 = 100 * epi_sum(sim, "infected.15to34", rows) / epi_sum(sim, "num.15to34", rows),
    hh_size = (pop_c + pop_a) / epi_mean(sim, "hh.num", rows),
    hh_alone = 100 * epi_sum(sim, "hh.alone", rows) / epi_sum(sim, "hh.num", rows),
    n_active = epi_mean(sim, c("i.num.child", "i.num.adult"), rows)
  )
}
cal <- calib(burn$baseline, rows_b)

# Concentration of transmission: the share of the window's infections caused
# by the 20% of cases who caused the most. Every transmission record carries
# the infector's unique id; cases who infected no one are counted from the
# number of incident cases in the window.
top20_share <- function(sim, rows) {
  sapply(seq_len(sim$control$nsims), function(s) {
    tm <- as.data.frame(get_transmat(sim, sim = s))
    tm <- tm[tm$at %in% rows, ]
    n_cases <- sum(sapply(c(inc_vars("child"), inc_vars("adult")),
                          function(v) sum(sim$epi[[v]][rows, s], na.rm = TRUE)))
    off <- sort(c(as.numeric(table(tm$infUid)),
                  rep(0, max(0, n_cases - length(unique(tm$infUid))))),
                decreasing = TRUE)
    100 * sum(off[seq_len(ceiling(0.2 * length(off)))]) / sum(off)
  })
}
cal$top20 <- top20_share(burn$baseline, rows_b)

targets <- data.frame(
  measure = c("TB incidence per 100,000 per year, all ages",
              "Share of incident TB in children under 15 (%)",
              "Share of incident TB diagnosed and treated (%)",
              "TB deaths as a share of incident TB (%)",
              "Infections of adults acquired in the household (%)",
              "Infections of children acquired in the household (%)",
              "Infections caused by the most infectious 20% of cases (%)",
              "Household contacts with active TB at index diagnosis (%)",
              "Household contacts ever infected at index diagnosis (%)",
              "Contacts under 5 with active TB at index diagnosis (%)",
              "Contacts under 5 ever infected at index diagnosis (%)",
              "Household contacts with TB onset within 2 years (%)",
              "Contacts under 5 with TB at investigation or within 2 years (%)",
              "Ever infected, ages 15 to 34 (%)",
              "Mean household size",
              "Single-person households (%)"),
  var = c("inc", "child_share", "detected", "cfr", "hh_adult", "hh_child",
          "top20", "hhc_active", "hhc_infected", "hhc_active_u5",
          "hhc_infected_u5", "hhc_risk", "hhc_tb_u5", "infected_15to34",
          "hh_size", "hh_alone"),
  target = c("about 400 (WHO 2024: South Africa 389, Indonesia 382)",
             "4 to 19 across high-burden countries (WHO 2024: South Africa 9, Indonesia 17)",
             "74 to 92 (WHO 2024 treatment coverage)",
             "6 to 12 in low-HIV high-burden countries (WHO 2024)",
             "8 to 19 of disease (molecular studies); 13 of disease (McCreesh and White)",
             "10 to 30 (Martinez 2019)",
             "about 80 to 90 (McCreesh and White: 93% of disease)",
             "3.1 (2.1 to 4.5), Fox 2013",
             "45 to 52 (65 in adults), Fox 2013; Velen 2021",
             "10.0 (5.0 to 18.9), Fox 2013",
             "35.5, Fox 2013",
             "about 2.5 (2.0 in the first year), Velen 2021",
             "7.6, prevalent plus incident, Martinez 2020",
             "50 to 80 in high-burden settings",
             "3.5 (South Africa, 2022 census)",
             "26 (South Africa, 2023)")
)
targets$model <- sapply(targets$var, function(v) fmt(mc(cal[[v]]), 1))
cat("\n=== Baseline (last 20 years of burn-in) against published targets ===\n")
print(targets[, c("measure", "model", "target")], right = FALSE)

## 5b. Part 1: the three contact structures -------------------------------------

struct_summary <- function(sim, rows) {
  layers <- sim$param$layer.names
  inf_l <- sapply(layers, function(l) {
    epi_sum(sim, paste0("inf.", l, c(".child", ".adult")), rows)
  })
  inf_l <- matrix(inf_l, ncol = length(layers), dimnames = list(NULL, layers))
  c_ <- calib(sim, rows)
  half <- split(rows, seq_along(rows) > length(rows) / 2)
  out <- list(
    inc = c_$inc, inc_child = c_$inc_child, inc_adult = c_$inc_adult,
    drift = 100 * (calib(sim, half[[2]])$inc / calib(sim, half[[1]])$inc - 1),
    reinf = 100 * epi_sum(sim, "inf.reinf", rows) /
      epi_sum(sim, c("inf.reinf", "inf.first"), rows),
    recent = c_$recent
  )
  for (l in layers) out[[paste0("share_", l)]] <- 100 * inf_l[, l] / rowSums(inf_l)
  inf_a <- matrix(sapply(layers, function(l) epi_sum(sim, paste0("inf.", l, ".adult"), rows)),
                  ncol = length(layers), dimnames = list(NULL, layers))
  for (l in layers) out[[paste0("adult_share_", l)]] <- 100 * inf_a[, l] / rowSums(inf_a)
  out
}
st <- lapply(burn, struct_summary, rows = rows_b)

# Percent change in incidence from the baseline structure, with a 95% Monte
# Carlo interval. The structures are separate burn-ins, so the simulations
# are not paired and the interval uses both groups' variances (it treats the
# baseline mean in the denominator as fixed).
rel_change <- function(x, x0) {
  d <- mean(x) - mean(x0)
  h <- if (length(x) > 1) {
    qt(0.975, 2 * length(x) - 2) * sqrt(var(x) / length(x) + var(x0) / length(x0))
  } else NA
  100 * c(est = d, lo = d - h, hi = d + h) / mean(x0)
}
struct_tbl <- data.frame(
  structure = names(st),
  incidence = sapply(st, function(x) fmt(mc(x$inc))),
  change = sapply(names(st), function(s) {
    if (s == "baseline") "reference" else fmt(rel_change(st[[s]]$inc, st$baseline$inc), 1)
  }),
  inc_child = sapply(st, function(x) fmt(mc(x$inc_child))),
  inc_adult = sapply(st, function(x) fmt(mc(x$inc_adult))),
  household = sapply(st, function(x) fmt(mc(if (is.null(x$share_household_redrawn))
    x$share_household else x$share_household_redrawn), 1)),
  regular = sapply(st, function(x) fmt(mc(x$share_regular), 1)),
  casual = sapply(st, function(x) fmt(mc(x$share_casual), 1)),
  reinfection = sapply(st, function(x) fmt(mc(x$reinf), 1)),
  row.names = NULL
)
cat("\n=== Part 1: incidence per 100,000 per year and share of infections by layer (%) ===\n")
print(struct_tbl, right = FALSE)

# Between-simulation variability of incidence in each structure
cv <- sapply(st, function(x) sd(x$inc) / mean(x$inc))
print(round(cv, 2))

# 6. Part 2: Interventions Resumed from the Baseline Burn-in -------------------

scenarios.df <- data.frame(
  .scenario.id = c("none", "hhci_screen", "hhci_u5", "hhci_all", "acf"),
  .at = 1,
  hhci = c(0, 1, 1, 1, 0),
  tpt.age.max = c(Inf, 0, 5, Inf, Inf),
  acf = c(0, 0, 0, 0, 1)
)
scenarios.list <- create_scenario_list(scenarios.df)

resume <- make_control(burnin + horizon, structures$baseline$tergm,
                       start = burnin + 1)
sims <- list()
for (scn in scenarios.list) {
  cat("Scenario:", scn$id, "\n")
  set.seed(202)
  sims[[scn$id]] <- netsim(burn$baseline, use_scenario(param_base, scn), init,
                           resume)
}

# Equal-effort comparator: screen people chosen at random, with the same
# package, at the average monthly rate at which household investigation
# screened contacts in the hhci_all scenario.
rows_h <- burnin + seq_len(horizon)
hh_screened <- mean(epi_sum(sims$hhci_all, "hh.screened", rows_h)) / horizon
pop_h <- mean(epi_sum(sims$hhci_all, c("num.child", "num.adult"), rows_h)) / horizon
random.rate <- hh_screened / pop_h
rnd <- create_scenario_list(data.frame(.scenario.id = "random", .at = 1,
                                       random.rate = random.rate,
                                       tpt.age.max = Inf))
set.seed(202)
sims$random <- netsim(burn$baseline, use_scenario(param_base, rnd[[1]]), init,
                      resume)


## 6a. Part 2 results -------------------------------------

scn_ids <- names(sims)
labels <- c(none = "Passive case finding only",
            hhci_screen = "Household investigation, no TPT",
            hhci_u5 = "Household investigation, TPT under 5",
            hhci_all = "Household investigation, TPT all ages",
            acf = "Community-wide screening of adults",
            random = "Random screening, equal effort")
bands <- c("u5", "5to14", "adult")

# Per-simulation outcomes over the horizon. The contact outcomes use years 3
# to 10: TB onsets within 24 months of being identified as a household
# contact of a newly diagnosed case, and the contacts identified, in the
# same months. In a steady state their ratio is the two-year risk.
rows_c <- burnin + (min(24, horizon - 1) + 1):horizon
rows_y4 <- burnin + (min(36, horizon - 1) + 1):min(48, horizon)  # year 4
outcomes <- function(sim) {
  pop <- epi_mean(sim, c("num.child", "num.adult"), rows_h)
  list(
    pop = pop,
    inc = epi_sum(sim, c(inc_vars("child"), inc_vars("adult")), rows_h),
    inc_child = epi_sum(sim, inc_vars("child"), rows_h),
    deaths = epi_sum(sim, c("tbdeath.child", "tbdeath.adult"), rows_h),
    screened = epi_sum(sim, c("hh.screened", "rnd.screened", "acf.screened"), rows_h),
    found = epi_sum(sim, c("hh.found", "rnd.found", "acf.found"), rows_h),
    tpt = epi_sum(sim, "tpt.completed", rows_h),
    prev_adult_y4 = epi_mean(sim, "i.num.adult", rows_y4),
    contacts_u5 = epi_sum(sim, "hhc.n.u5", rows_c),
    contacts_5p = epi_sum(sim, c("hhc.n.5to14", "hhc.n.adult"), rows_c),
    contact_tb_u5 = epi_sum(sim, "hhc.inc.u5", rows_c),
    contact_tb_5p = epi_sum(sim, c("hhc.inc.5to14", "hhc.inc.adult"), rows_c),
    tpt_c = epi_sum(sim, "tpt.completed", rows_c)
  )
}
oc <- lapply(sims, outcomes)
per100k <- function(x, pop) 1e5 * x / pop

# Population impact: differences from passive case finding, paired by
# simulation (every scenario resumes simulation s from the same burn-in
# state s), with 95% Monte Carlo intervals.
impact_tbl <- data.frame(
  scenario = labels[scn_ids],
  incidence = sapply(scn_ids, function(s)
    fmt(mc(per100k_yr(oc[[s]]$inc, oc[[s]]$pop, rows_h)))),
  tb_averted = sapply(scn_ids, function(s)
    fmt(mc(per100k(oc$none$inc - oc[[s]]$inc, oc$none$pop)))),
  pct_averted = sapply(scn_ids, function(s)
    fmt(mc(100 * (oc$none$inc - oc[[s]]$inc) / mean(oc$none$inc)), 1)),
  child_tb_averted = sapply(scn_ids, function(s)
    fmt(mc(per100k(oc$none$inc_child - oc[[s]]$inc_child, oc$none$pop)))),
  deaths_averted = sapply(scn_ids, function(s)
    fmt(mc(per100k(oc$none$deaths - oc[[s]]$deaths, oc$none$pop)))),
  row.names = NULL
)
cat("\n=== Part 2: population impact over", horizon / 12,
    "years (incidence per 100,000 per year; averted per 100,000 population) ===\n")
print(impact_tbl, right = FALSE)

# Program yield. Screenings count every screen, so a person screened in
# several community rounds counts several times. TB averted per case found
# divides the mean cases averted by the mean cases found by the program.
program_tbl <- data.frame(
  scenario = labels[scn_ids],
  screenings = sapply(scn_ids, function(s) round(mean(per100k(oc[[s]]$screened, oc[[s]]$pop)))),
  found = sapply(scn_ids, function(s) round(mean(per100k(oc[[s]]$found, oc[[s]]$pop)))),
  found_per_1000 = sapply(scn_ids, function(s)
    if (sum(oc[[s]]$screened) == 0) NA else
      round(1000 * sum(oc[[s]]$found) / sum(oc[[s]]$screened), 1)),
  tpt_completed = sapply(scn_ids, function(s) round(mean(per100k(oc[[s]]$tpt, oc[[s]]$pop)))),
  averted_per_found = sapply(scn_ids, function(s)
    if (sum(oc[[s]]$found) == 0) NA else
      round(sum(oc$none$inc - oc[[s]]$inc) / sum(oc[[s]]$found), 1)),
  row.names = NULL
)
cat("\n=== Part 2: program yield over", horizon / 12, "years (per 100,000 population) ===\n")
print(program_tbl, right = FALSE)

# The direct effect of TPT on household contacts. The three household arms
# screen the same way and differ only in who is offered TPT, so contacts
# under 5 are compared between no TPT and TPT under 5, and contacts 5 and
# older between TPT under 5 and TPT for all ages. Risks are TB onsets within
# 24 months per 1,000 contacts identified, pooled over simulations.
risk <- function(o, grp) 1000 * sum(o[[paste0("contact_tb_", grp)]]) /
  sum(o[[paste0("contacts_", grp)]])
tpt_effect <- function(without, with, grp, courses) {
  r0 <- risk(oc[[without]], grp)
  r1 <- risk(oc[[with]], grp)
  averted <- (r0 - r1) / 1000 * sum(oc[[with]][[paste0("contacts_", grp)]])
  c(risk_without = round(r0, 1), risk_with = round(r1, 1),
    pct_reduction = round(100 * (r0 - r1) / r0),
    courses_per_case_averted = if (isTRUE(averted > 0)) round(courses / averted) else NA)
}
tpt_tbl <- rbind(
  "Contacts under 5" = tpt_effect("hhci_screen", "hhci_u5", "u5",
                                  sum(oc$hhci_u5$tpt_c)),
  "Contacts 5 and older" = tpt_effect("hhci_u5", "hhci_all", "5p",
                                      sum(oc$hhci_all$tpt_c) - sum(oc$hhci_u5$tpt_c))
)
cat("\n=== Part 2: TB within 2 years per 1,000 household contacts, without and with TPT ===\n")
print(tpt_tbl)

# The most a household program could prevent directly: the share of all
# incident TB, under passive case finding, that occurs in people identified
# as a household contact in the previous 24 months
contact_share <- 100 * (oc$none$contact_tb_u5 + oc$none$contact_tb_5p) /
  epi_sum(sims$none, c(inc_vars("child"), inc_vars("adult")), rows_c)
cat(sprintf("\nIncident TB in recent household contacts, passive case finding: %s%% of all incident TB\n",
            fmt(mc(contact_share), 1)))

# Community-wide screening against ACT3: adult prevalence of active TB in
# year 4 (after three annual rounds), relative to passive case finding
acf_prev_ratio <- mc(oc$acf$prev_adult_y4 / oc$none$prev_adult_y4)
cat(sprintf("Adult TB prevalence in year 4, community screening vs passive: %s\n",
            fmt(acf_prev_ratio, 2)))


# 7. Plots ---------------------------------------------------------------------

cols_st <- c(baseline = "black", regular_redrawn = "#e67e22",
             household_redrawn = "#c0392b")
annual <- function(sim, vars, rows) {
  yr <- ceiling(seq_along(rows) / 12)
  ev <- Reduce(`+`, lapply(vars, function(v) as.matrix(sim$epi[[v]])[rows, , drop = FALSE]))
  pop <- Reduce(`+`, lapply(c("num.child", "num.adult"),
                            function(v) as.matrix(sim$epi[[v]])[rows, , drop = FALSE]))
  sapply(split(seq_along(rows), yr), function(i) {
    mean(1e5 * colSums(ev[i, , drop = FALSE]) / colMeans(pop[i, , drop = FALSE]))
  })
}
inc_all <- c(inc_vars("child"), inc_vars("adult"))

## Plot 1: incidence over the burn-in by contact structure
inc_st <- lapply(burn, annual, vars = inc_all, rows = 2:burnin)
par(mfrow = c(1, 1), mar = c(4, 4.5, 2, 1))
plot(NA, xlim = c(1, length(inc_st$baseline)),
     ylim = c(0, max(unlist(inc_st)) * 1.05),
     xlab = "Year of burn-in", ylab = "TB incidence per 100,000 per year",
     main = "Incidence by Contact Structure")
for (s in names(inc_st)) lines(inc_st[[s]], col = cols_st[s], lwd = 2)
legend("topright", legend = names(inc_st), col = cols_st, lwd = 2, bty = "n")

## Plot 2: contact time versus transmission by layer (baseline)
time_share <- 100 * hours_day / sum(hours_day)
inf_share <- c(household = mean(st$baseline$adult_share_household),
               regular = mean(st$baseline$adult_share_regular),
               casual = mean(st$baseline$adult_share_casual))
barplot(rbind(time_share, inf_share), beside = TRUE,
        col = c("gray70", "#2c7fb8"), ylim = c(0, 100),
        ylab = "Percent", main = "Contact Time and Transmission by Layer")
legend("topleft", fill = c("gray70", "#2c7fb8"), bty = "n",
       legend = c("Share of an adult's contact time",
                  "Share of adults' infections acquired"))

## Plot 3: incidence over the horizon by intervention scenario
cols_scn <- c(none = "black", hhci_screen = "#a6bddb", hhci_u5 = "#74a9cf",
              hhci_all = "#0570b0", acf = "#d95f0e", random = "#756bb1")
rows_ctx <- max(2, burnin - 59):(burnin + horizon)
inc_scn <- lapply(sims, annual, vars = inc_all, rows = rows_ctx)
yrs <- seq_along(inc_scn$none) - ceiling((burnin - min(rows_ctx) + 1) / 12)
plot(NA, xlim = range(yrs), ylim = c(0, max(unlist(inc_scn)) * 1.05),
     xlab = "Years since programs started", ylab = "TB incidence per 100,000 per year",
     main = "Incidence by Intervention")
abline(v = 0.5, lty = 3)
for (s in names(inc_scn)) lines(yrs, inc_scn[[s]], col = cols_scn[s], lwd = 2)
legend("bottomleft", legend = labels[names(inc_scn)], col = cols_scn[names(inc_scn)],
       lwd = 2, bty = "n", cex = 0.8)

## Plot 4: population impact by program, one point per simulation
scn_int <- setdiff(scn_ids, "none")
short <- c(hhci_screen = "Household, no TPT", hhci_u5 = "Household, TPT <5",
           hhci_all = "Household, TPT all", acf = "Community-wide",
           random = "Random, equal effort")
diffs <- matrix(sapply(scn_int, function(s) per100k(oc$none$inc - oc[[s]]$inc, oc$none$pop)),
                ncol = length(scn_int), dimnames = list(NULL, scn_int))
par(mfrow = c(1, 1), mar = c(9, 5, 3, 1))
plot(NA, xlim = c(0.5, length(scn_int) + 0.5), ylim = range(c(0, diffs)) * 1.1,
     xaxt = "n", xlab = "", ylab = "TB averted per 100,000 over 10 years",
     main = "Population Impact")
abline(h = 0, lty = 3)
for (j in seq_along(scn_int)) {
  points(j + runif(nrow(diffs), -0.15, 0.15), diffs[, j], pch = 19,
         col = adjustcolor(cols_scn[scn_int[j]], 0.6))
  ci <- mc(diffs[, j])
  segments(j - 0.25, ci["est"], j + 0.25, ci["est"], lwd = 3)
  if (!is.na(ci["lo"])) arrows(j, ci["lo"], j, ci["hi"], angle = 90, code = 3, length = 0.05)
}
axis(1, at = seq_along(scn_int), labels = short[scn_int], las = 2, cex.axis = 0.85)

## Plot 5: yield of screening by program
yield <- program_tbl$found_per_1000[match(scn_int, scn_ids)]
bp2 <- barplot(yield, names.arg = short[scn_int], las = 2, col = cols_scn[scn_int],
               cex.names = 0.85, ylim = c(0, max(yield, na.rm = TRUE) * 1.2),
               ylab = "Active TB found per 1,000 screened", main = "Yield of Screening")
text(bp2, yield + max(yield, na.rm = TRUE) * 0.05, yield, cex = 0.9)
