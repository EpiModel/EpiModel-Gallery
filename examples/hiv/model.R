
##
## HIV Transmission Model with Care Cascade and PrEP
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##

# Load EpiModel
suppressMessages(library(EpiModel))

# Standard Gallery unit test lines
rm(list = ls())
eval(parse(text = print(commandArgs(TRUE)[1])))

if (interactive()) {
  N <- 1500
  nsims <- 5
  ncores <- 5
  nsteps <- 1040          # 20 years on a weekly timestep
  dx_sims <- 20           # network diagnostics: 20 runs of 2,000 weeks
  dx_steps <- 2000
} else {
  N <- 300
  nsims <- 2
  ncores <- 2
  nsteps <- 100
  dx_sims <- 2
  dx_steps <- 100
}

set.seed(2026)


# 1. Network Setup -----------------------------------------------------------

# Two persistent partnership layers represent the HIV-specific partnership
# pattern that an edges-only Bernoulli graph cannot:
#
#   - Main:   low-degree, long-duration, high per-week act count. Most of
#             an individual's infectious time is spent under prolonged
#             exposure to a single main partner. Capped at degree 2 by
#             degrange because "5 simultaneous main partners" is not a
#             plausible main-partnership structure.
#   - Casual: lower mean degree, shorter-duration, fewer acts per week,
#             but degree distribution is left UNCAPPED. The right tail
#             (a small number of nodes with many concurrent casual ties)
#             is precisely the high-risk subgroup that PrEP targeting is
#             designed for; truncating it would erase the structural
#             heterogeneity the intervention is meant to address.
#
# Both layers include a `concurrent` ERGM term to encode overlapping
# partnerships, the structural feature most strongly implicated in
# generalized HIV epidemics (Morris & Kretzschmar 1997).

departure_rate <- 0.0005      # ~38-year mean tenure in the sexually-active pop

# Partnerships also end when a partner leaves the population, and
# dissolution_coefs() adjusts each layer's dissolution model for that. The
# rate it needs is the total departure rate the network sees: background
# exits plus AIDS deaths, about 0.00065 per week in the untreated epidemic
# (the same value as the arrival rate below). With the background rate alone,
# the main layer ran about 3% below its mean degree target in netsim().
network_exit_rate <- 0.00065

# A tighter convergence tolerance for the ERGM fits. With the default, the
# casual fit came out about 1% below its targets; this costs a few seconds.
fit_ctrl <- control.ergm(MCMLE.MCMC.precision = 0.01)

nw <- network_initialize(N)

# --- Main partnership layer ---
mean_deg_main <- 0.5
concurrent_main <- round(0.04 * N)    # ~4% of nodes have 2+ main partners
formation_main <- ~edges + concurrent + degrange(from = 3)
target_main <- c(mean_deg_main * N / 2, concurrent_main, 0)
diss_main <- dissolution_coefs(~offset(edges), duration = 200,
                               d.rate = network_exit_rate)

# --- Casual partnership layer ---
# No degrange cap: the right tail of casual degree is the high-activity
# subgroup that drives transmission and motivates targeted PrEP.
# Concurrent target 10% is ~2.7x the Poisson baseline for mean degree 0.3
# (which would give 1 - exp(-0.3) * 1.3 ~= 3.7% concurrent), so the term
# meaningfully shifts the network above the Bernoulli-graph default
# without over-constraining the degree distribution.
mean_deg_cas <- 0.3
concurrent_cas <- round(0.10 * N)
formation_cas <- ~edges + concurrent
target_cas <- c(mean_deg_cas * N / 2, concurrent_cas)
diss_cas <- dissolution_coefs(~offset(edges), duration = 26,
                              d.rate = network_exit_rate)

cat("Fitting main partnership layer ERGM...\n")
est_main <- netest(nw, formation_main, target_main, diss_main,
                   set.control.ergm = fit_ctrl, verbose = FALSE)
cat("Fitting casual partnership layer ERGM...\n")
est_cas <- netest(nw, formation_cas, target_cas, diss_cas,
                  set.control.ergm = fit_ctrl, verbose = FALSE)


# 2. Network Diagnostics ------------------------------------------------------

# tergm draws each week's network with a Markov chain that starts from the
# previous week's network and, by default, may stop after 1,000 proposals.
# That is too short for the casual layer, whose share of people with 2+
# partners settled about 5% below target. Fixing the chain at 10,000
# proposals a week gives the same diagnostics as 100,000. The same control
# goes to netdx() here and to netsim() below.
chain_ctrl <- control.simulate.formula.tergm(MCMC.burnin.min = 1e4,
                                             MCMC.burnin.max = 1e4)

# Dynamic diagnostics simulate each layer on its own, in a closed
# population, and compare the simulated statistics with the targets.
cat("Running network diagnostics...\n")
dx_main <- netdx(est_main, nsims = dx_sims, nsteps = dx_steps,
                 ncores = ncores, set.control.tergm = chain_ctrl,
                 verbose = FALSE)
dx_cas <- netdx(est_cas, nsims = dx_sims, nsteps = dx_steps,
                ncores = ncores, set.control.tergm = chain_ctrl,
                verbose = FALSE)
print(dx_main)
print(dx_cas)

par(mfrow = c(1, 2), mar = c(4, 4, 2.5, 1), mgp = c(2.5, 0.8, 0))
plot(dx_main, legend = TRUE)
title("Main layer")
plot(dx_cas, legend = TRUE)
title("Casual layer")


# 3. Disease and Intervention Parameters -------------------------------------

source("examples/hiv/module-fx.R")

# Per-act transmission probability is a single biological parameter; per-
# step per-partnership probability is computed in infect() from the per-
# layer act rate. Stage and ART-status multipliers reflect viral load:
#
#   Acute     (~12 wk)       :  5x  -- elevated VL in the seroconversion window
#                                       (Hollingsworth et al., 2008)
#   Chronic   (~10 yr)       :  1x  -- reference
#   AIDS      (~2 yr untx)   :  2x  -- some viremic rebound, but in practice
#                                       counterbalanced by reduced sexual
#                                       activity in late-stage disease
#   ART, not yet suppressed  : 0.3x -- 2-3 mo before VL becomes undetectable
#   ART, virally suppressed  : 0.01x-- U=U: HPTN 052 and PARTNER established
#                                       that durable suppression eliminates
#                                       sexual transmission risk
#
# PrEP efficacy (95%) reflects the high end of demonstrated effectiveness
# with consistent adherence (oral TDF/FTC and injectable cabotegravir).

# Base parameter set holds the defaults for every parameter the modules
# read. Per-scenario overrides are applied via the EpiModel scenarios
# API (create_scenario_list + use_scenario) further below. Cascade and
# PrEP rates default to zero so that the "none" scenario inherits an
# all-off baseline directly from this base set.
param_base <- param.net(
  # --- Transmission biology ---
  inf.prob.act = 0.0025,
  rel.inf.acute = 5,
  rel.inf.aids = 2,
  rel.inf.art.unsupp = 0.30,
  rel.inf.art.supp = 0.01,
  prep.efficacy = 0.95,
  acts.main = 3,
  acts.casual = 1,
  # --- Disease progression (per week) ---
  acute.to.chronic.rate = 1 / 12,    # 12-week acute window
  chronic.to.aids.rate = 1 / 520,    # ~10-year chronic stage
  aids.depart.rate = 1 / 104,        # ~2-year AIDS survival untreated
  art.prog.mult = 0.5,               # suppressive ART halves progression
  art.aids.surv.mult = 0.1,          # ART extends AIDS survival 10x
  # --- Care cascade (off by default) ---
  test.rate = 0,
  aids.dx.rate = 0,
  linkage.rate = 0,
  art.reinit.rate = 0,
  suppression.rate = 0,
  art.disc.rate = 0,
  # --- PrEP (off by default; indication threshold set up front) ---
  prep.init.cov = 0,
  prep.start.rate = 0,
  prep.stop.rate = 0,
  prep.indic.deg = 2,
  # --- Vital dynamics ---
  departure.rate = departure_rate,
  arrival.rate = 0.00065             # offsets background + AIDS mortality
)

# Initial conditions: 8% seroprevalence, with the seed cohort distributed
# across stages by mean stage duration (handled in init_attrs).
init <- init.net(i.num = round(0.08 * N))

control <- control.net(
  type = NULL,
  nsims = nsims,
  ncores = ncores,
  nsteps = nsteps,
  tergmLite = TRUE,
  resimulate.network = TRUE,
  infection.FUN = infect,
  progress.FUN = progress,
  cascade.FUN = cascade,
  prep.FUN = prep,
  departures.FUN = dfunc,
  arrivals.FUN = afunc,
  set.control.tergm = chain_ctrl,
  verbose = FALSE
)


# 4. Scenarios --------------------------------------------------------------

# Scenarios are defined as a flat data frame (one row per scenario) and
# converted to a list via EpiModel's create_scenario_list(). At sim time,
# use_scenario(param_base, scn) returns a parameter object with the
# scenario's row applied to the base defaults. .scenario.id labels the
# scenario; .at = 0 means the override applies before step 1.
#
# Cascade rates are calibrated so the cascade scenario lands near the
# UNAIDS 95-95-95 attainment (95% of PLHIV diagnosed, 95% of those on
# ART, 95% of those on ART suppressed). PrEP scenarios use risk-based
# eligibility set in param_base (prep.indic.deg = 2) with equal start and
# stop rates, which would give 50% coverage in a fixed indicated group.
# Indication changes as partnerships form and end, so realized coverage
# among indicated susceptibles is lower (about a third).
scenarios.df <- data.frame(
  .scenario.id     = c("none", "cascade", "prep",  "both"),
  .at              = 0,
  test.rate        = c(0,       0.015,     0,       0.015),
  aids.dx.rate     = c(0,       0.050,     0,       0.050),
  linkage.rate     = c(0,       0.100,     0,       0.100),
  art.reinit.rate  = c(0,       0.030,     0,       0.030),
  suppression.rate = c(0,       1 / 12,    0,       1 / 12),
  art.disc.rate    = c(0,       0.002,     0,       0.002),
  prep.init.cov    = c(0,       0,         0.50,    0.50),
  prep.start.rate  = c(0,       0,         0.015,   0.015),
  prep.stop.rate   = c(0,       0,         0.015,   0.015)
)
scenarios.list <- create_scenario_list(scenarios.df)

labels <- c(none = "No intervention",
            cascade = "Cascade (95-95-95)",
            prep = "PrEP (50% coverage)",
            both = "Cascade + PrEP")

sims <- list()
for (scn in scenarios.list) {
  cat(sprintf("Scenario: %s\n", scn$id))
  sims[[scn$id]] <- netsim(list(est_main, est_cas),
                           use_scenario(param_base, scn),
                           init, control)
}


# 5. Network Targets During the Epidemic ------------------------------------

# netdx() checks each layer in a closed population. The epidemic model adds
# arrivals, background departures, and AIDS deaths, so check the targets
# during netsim() as well. netsim() records each layer's formation statistics
# every week (save.nwstats = TRUE by default); dividing by the population
# size each week gives per-person quantities that stay comparable as the
# population changes.
net_track <- function(sim) {
  num <- as.data.frame(sim)[, c("time", "sim", "num")]
  out <- NULL
  for (k in 1:2) {
    ns <- merge(get_nwstats(sim, network = k), num, by = c("time", "sim"))
    out <- rbind(out, data.frame(layer = c("main", "casual")[k],
                                 time = ns$time, sim = ns$sim,
                                 deg = 2 * ns$edges / ns$num,
                                 conc = ns$concurrent / ns$num))
  }
  out
}
nw_track <- lapply(sims, net_track)

net_targets <- data.frame(
  layer = c("main", "main", "casual", "casual"),
  stat = c("deg", "conc", "deg", "conc"),
  label = c("Main: mean degree", "Main: share with 2+ partners",
            "Casual: mean degree", "Casual: share with 2+ partners"),
  target = c(mean_deg_main, concurrent_main / N,
             mean_deg_cas, concurrent_cas / N)
)

# Percent difference from target, averaged over all weeks and simulations.
pct_off <- sapply(nw_track, function(d) {
  sapply(seq_len(nrow(net_targets)), function(i) {
    x <- d[d$layer == net_targets$layer[i], net_targets$stat[i]]
    100 * (mean(x) / net_targets$target[i] - 1)
  })
})
nw_tbl <- data.frame(Statistic = net_targets$label,
                     Target = net_targets$target,
                     round(pct_off, 1), check.names = FALSE)
names(nw_tbl)[-(1:2)] <- labels[colnames(pct_off)]
cat("\n=== Network statistics in netsim(): % difference from target ===\n")
print(nw_tbl, row.names = FALSE)

cols_scn <- c(none = "firebrick", cascade = "steelblue",
              prep = "darkorange", both = "darkgreen")

par(mfrow = c(2, 2), mar = c(4, 4, 2.5, 1), mgp = c(2.5, 0.8, 0),
    oma = c(3, 0, 0, 0))
for (i in seq_len(nrow(net_targets))) {
  yearly <- lapply(nw_track, function(d) {
    d <- d[d$layer == net_targets$layer[i], ]
    tapply(d[[net_targets$stat[i]]], ceiling(d$time / 52), mean)
  })
  tg <- net_targets$target[i]
  plot(NULL, xlim = c(1, ceiling(nsteps / 52)), ylim = tg * c(0.8, 1.2),
       xlab = "Year", ylab = "Simulated value", main = net_targets$label[i])
  abline(h = tg, lty = 2, col = "gray40")
  for (s in names(yearly)) {
    lines(as.numeric(names(yearly[[s]])), yearly[[s]],
          col = cols_scn[s], lwd = 2)
  }
}
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0),
    new = TRUE)
plot.new()
legend("bottom", legend = c(labels, "Target"), ncol = 3,
       col = c(cols_scn, "gray40"), lwd = c(2, 2, 2, 2, 1),
       lty = c(1, 1, 1, 1, 2), bty = "n", cex = 0.8)


# 6. Analysis ---------------------------------------------------------------

# Derived measures: prevalence, incidence rate per 100 person-years, and
# UNAIDS cascade attainment (% diagnosed, % on ART, % suppressed).
# NB: mutate_epi rejects ifelse, so we do plain division and let NAs
# propagate where i.num == 0 (handled later via na.rm).
for (s in names(sims)) {
  sims[[s]] <- mutate_epi(sims[[s]],
    prev = i.num / num,
    ir100py = 52 * 100 * si.flow / s.num,
    pct.dx = dx.num / i.num,
    pct.art = art.num / i.num,
    pct.supp = supp.num / i.num,
    prep.cov = prep.num / s.num
  )
}

# Headline summary: cumulative infections + final-state metrics.
df_summary <- function(s) {
  df <- as.data.frame(sims[[s]])
  last_t <- max(df$time)
  late <- df$time >= 0.5 * last_t   # average over the back half (post-burn-in)
  data.frame(
    Scenario = labels[s],
    Cum_infections = round(sum(df$si.flow, na.rm = TRUE) / nsims),
    Final_prev = round(mean(df$prev[df$time == last_t], na.rm = TRUE), 3),
    Mean_inc_100py = round(mean(df$ir100py[late], na.rm = TRUE), 2),
    Pct_diagnosed = round(100 * mean(df$pct.dx[late], na.rm = TRUE), 1),
    Pct_on_ART = round(100 * mean(df$pct.art[late], na.rm = TRUE), 1),
    Pct_suppressed = round(100 * mean(df$pct.supp[late], na.rm = TRUE), 1),
    stringsAsFactors = FALSE
  )
}
summary_tbl <- do.call(rbind, lapply(names(sims), df_summary))
rownames(summary_tbl) <- NULL
cat("\n=== Scenario summary ===\n")
print(summary_tbl, row.names = FALSE)

# Infections averted versus the no-intervention baseline.
none_cum <- summary_tbl$Cum_infections[1]
summary_tbl$Infections_averted <- none_cum - summary_tbl$Cum_infections
cat("\n=== Infections averted vs. no intervention ===\n")
print(summary_tbl[, c("Scenario", "Cum_infections", "Infections_averted")],
      row.names = FALSE)


# 7. Plots ------------------------------------------------------------------

## --- Plot 1: HIV prevalence over time (the headline) ---
par(mfrow = c(1, 1), mar = c(4, 4, 2.5, 1), mgp = c(2.5, 0.8, 0))
prev_max <- max(sapply(sims, function(s) {
  df <- as.data.frame(s); max(df$prev, na.rm = TRUE)
}))
plot(sims[[1]], y = "prev",
     main = "HIV Prevalence by Scenario",
     ylab = "Prevalence", xlab = "Week",
     mean.col = cols_scn[1], mean.lwd = 2, mean.smooth = TRUE,
     qnts = FALSE, legend = FALSE,
     ylim = c(0, 1.05 * prev_max))
for (i in 2:length(sims)) {
  plot(sims[[i]], y = "prev", add = TRUE,
       mean.col = cols_scn[i], mean.lwd = 2, mean.smooth = TRUE,
       qnts = FALSE, legend = FALSE)
}
legend("topleft", legend = labels, col = cols_scn, lwd = 2,
       bty = "n", cex = 0.9)


## --- Plot 2: HIV incidence per 100 person-years ---
inc_max <- max(sapply(sims, function(s) {
  df <- as.data.frame(s); max(df$ir100py, na.rm = TRUE)
}))
plot(sims[[1]], y = "ir100py",
     main = "HIV Incidence Rate by Scenario",
     ylab = "New infections per 100 person-years",
     xlab = "Week",
     mean.col = cols_scn[1], mean.lwd = 2, mean.smooth = TRUE,
     qnts = FALSE, legend = FALSE,
     ylim = c(0, 1.05 * inc_max))           # lower bound at 0 prevents
for (i in 2:length(sims)) {                  # the smoothed cascade line
  plot(sims[[i]], y = "ir100py", add = TRUE, # from being clipped below 0
       mean.col = cols_scn[i], mean.lwd = 2, mean.smooth = TRUE,
       qnts = FALSE, legend = FALSE)
}
legend("topright", legend = labels, col = cols_scn, lwd = 2,
       bty = "n", cex = 0.9)


## --- Plot 3: Cascade attainment (cascade scenario only) ---
# Single panel because the cascade attainment in the "both" scenario
# is identical by construction (same cascade rates). The reference
# lines mark the UNAIDS 95-95-95 cumulative targets: 95% of PLHIV
# diagnosed, 90.25% on ART, 85.7% virally suppressed.
par(mfrow = c(1, 1), mar = c(4, 4, 2.5, 1), mgp = c(2.5, 0.8, 0))
plot(sims[["cascade"]], y = "pct.dx",
     main = "Care Cascade Attainment (cascade scenario)",
     ylab = "Fraction of PLHIV", xlab = "Week",
     mean.col = "steelblue", mean.lwd = 2, mean.smooth = TRUE,
     qnts = FALSE, legend = FALSE, ylim = c(0, 1))
plot(sims[["cascade"]], y = "pct.art", add = TRUE,
     mean.col = "darkorange", mean.lwd = 2, mean.smooth = TRUE,
     qnts = FALSE, legend = FALSE)
plot(sims[["cascade"]], y = "pct.supp", add = TRUE,
     mean.col = "darkgreen", mean.lwd = 2, mean.smooth = TRUE,
     qnts = FALSE, legend = FALSE)
abline(h = c(0.95, 0.95 ^ 2, 0.95 ^ 3), lty = 3, col = "gray50")
legend("bottomright",
       legend = c("% Diagnosed", "% on ART", "% Suppressed",
                  "UNAIDS 95-95-95 targets"),
       col = c("steelblue", "darkorange", "darkgreen", "gray50"),
       lwd = c(2, 2, 2, 1), lty = c(1, 1, 1, 3),
       bty = "n", cex = 0.9)


## --- Plot 4: PrEP coverage among indicated vs. all susceptibles ---
# The interesting comparison is between total susceptibles (PrEP fraction
# is small because most susceptibles don't meet indication criteria) and
# the indicated subgroup (PrEP fraction is high because that's where the
# intervention concentrates).
for (s in c("prep", "both")) {
  sims[[s]] <- mutate_epi(sims[[s]],
    prep.cov.indic = prep.num.indic / prep.indic.num
  )
}
plot(sims[["prep"]], y = "prep.cov.indic",
     main = "PrEP Coverage Among Indicated vs. All Susceptibles",
     ylab = "Fraction on PrEP", xlab = "Week",
     mean.col = "darkorange", mean.lwd = 2, mean.smooth = TRUE,
     qnts = FALSE, legend = FALSE, ylim = c(0, 1))
plot(sims[["prep"]], y = "prep.cov", add = TRUE,
     mean.col = "darkorange", mean.lwd = 2, mean.smooth = TRUE,
     mean.lty = 2, qnts = FALSE, legend = FALSE)
legend("topright",
       legend = c("Indicated susceptibles (target population)",
                  "All susceptibles (population coverage)"),
       col = c("darkorange", "darkorange"),
       lwd = 2, lty = c(1, 2), bty = "n", cex = 0.9)
