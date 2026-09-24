
##
## Partner Notification for an Endemic Bacterial STI
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
#   Rscript examples/sis-partner-notification/model.R "run_full <- TRUE"
# The full run takes about five minutes on five cores; CI mode runs in well
# under a minute and its results are not meant to be interpreted.
if (interactive() || exists("run_full")) {
  N <- 10000
  nsims <- 10
  ncores <- 5
  burnin <- 520
  horizon <- 260
} else {
  N <- 1000
  nsims <- 1
  ncores <- 1
  burnin <- 52
  horizon <- 26
}


# 1. Network Model Estimation ----------------------------------------------

# A heterosexual population of young adults, one time step per week, with
# two partnership layers. Attributes belong to the nodes, and the layers share
# one node set, so every attribute is set on a single network object that both
# layers are estimated from, whether or not a layer's model uses it.
#   sex   F or M; both layers allow only female-male partnerships
#   risk  H for the 10% of people with high casual-partner activity
set.seed(2026)
nw <- network_initialize(N)
nw <- set_vertex_attribute(nw, "sex", rep(c("F", "M"), length.out = N))
risk <- ifelse(seq_len(N) <= 0.1 * N, "H", "L")
nw <- set_vertex_attribute(nw, "risk", risk)

# People leave the population at 1/520 per week (a mean of 10 years in the 15
# to 24 age band), and the dissolution models are adjusted for that exit.
departure_rate <- 1 / 520

# Main partnerships: half the population is in one at any time, nobody has
# two at once, and they last 78 weeks on average. The offset terms with
# coefficients of -Inf are hard constraints.
formation_main <- ~edges + offset(nodematch("sex")) + offset(concurrent)
target_main <- 0.5 * N / 2
diss_main <- dissolution_coefs(~offset(edges), duration = 78,
                               d.rate = departure_rate)
est_main <- netest(nw, formation_main, target_main, diss_main,
                   coef.form = c(-Inf, -Inf), verbose = FALSE)

# Casual partnerships last 8 weeks on average. High-activity people have a
# mean of 0.5 casual partners at any time and everyone else 0.04, and half of
# the casual partnerships of high-activity people are with each other.
n_high <- sum(risk == "H")
deg_high <- 0.5
deg_low <- 0.04
formation_cas <- ~edges + offset(nodematch("sex")) +
  nodefactor("risk", levels = "H") + nodematch("risk", diff = TRUE, levels = "H")
target_cas <- round(c((n_high * deg_high + (N - n_high) * deg_low) / 2,
                      n_high * deg_high,
                      0.5 * n_high * deg_high / 2))
diss_cas <- dissolution_coefs(~offset(edges), duration = 8,
                              d.rate = departure_rate)
est_cas <- netest(nw, formation_cas, target_cas, diss_cas,
                  coef.form = -Inf, verbose = FALSE)

# Dynamic diagnostics. tergm draws each week's network with a short Markov
# chain, and its default length (MCMC.burnin.min = 1000) is too short for a
# layer with 8-week partnerships concentrated in a small group: at N = 10,000
# the casual layer settles about a fifth below its edge target and far below
# its high-activity targets. A longer chain fixes this at little cost, and
# the same setting goes to netsim() below.
tergm_ctrl <- control.simulate.formula.tergm(MCMC.burnin.min = 1e5)
dx_main <- netdx(est_main, nsims = 5, ncores = ncores, nsteps = 260,
                 nwstats.formula = ~edges + concurrent,
                 set.control.tergm = tergm_ctrl, verbose = FALSE)
dx_cas <- netdx(est_cas, nsims = 5, ncores = ncores, nsteps = 260,
                nwstats.formula = ~edges + nodefactor("risk", levels = TRUE) +
                  nodematch("risk", diff = TRUE, levels = "H"),
                set.control.tergm = tergm_ctrl, verbose = FALSE)
print(dx_main)
print(dx_cas)

# Realized mean casual degree by group, and the share of the casual
# partnerships of high-activity people that are with each other
nws <- colMeans(get_nwstats(dx_cas)[, c("nodefactor.risk.H", "nodefactor.risk.L",
                                        "nodematch.risk.H")])
deg_tbl <- data.frame(
  Quantity = c("Mean casual degree, high activity",
               "Mean casual degree, everyone else",
               "High-activity casual partnerships with each other (%)"),
  Target = as.character(c(deg_high, deg_low, 50)),
  Simulated = as.character(c(round(nws[["nodefactor.risk.H"]] / n_high, 3),
                             round(nws[["nodefactor.risk.L"]] / (N - n_high), 3),
                             round(100 * 2 * nws[["nodematch.risk.H"]] /
                                     nws[["nodefactor.risk.H"]], 1)))
)
print(deg_tbl)

# The lesson in numbers: the same casual-layer diagnostic with tergm's default
# chain length, against the targets and the longer chain
dx_cas_default <- netdx(est_cas, nsims = 5, ncores = ncores, nsteps = 260,
                        nwstats.formula = ~edges +
                          nodefactor("risk", levels = "H") +
                          nodematch("risk", diff = TRUE, levels = "H"),
                        verbose = FALSE)
stat_names <- c("edges", "nodefactor.risk.H", "nodematch.risk.H")
chain_mean <- function(dx) colMeans(get_nwstats(dx)[, stat_names])
pct_off <- function(m) sprintf("%.0f (%+.0f%%)", m, 100 * (m / target_cas - 1))
chain_tbl <- data.frame(
  Statistic = c("Casual partnerships", "Casual partnership ends, high activity",
                "Casual partnerships, high with high"),
  Target = target_cas,
  `Default chain` = pct_off(chain_mean(dx_cas_default)),
  `MCMC.burnin.min = 1e5` = pct_off(chain_mean(dx_cas)),
  check.names = FALSE, row.names = NULL
)
print(chain_tbl)


# 2. Parameters, Initial Conditions, and Controls --------------------------

if (file.exists("examples/sis-partner-notification/module-fx.R")) {
  source("examples/sis-partner-notification/module-fx.R")
} else {
  source("module-fx.R")
}

# Natural history and testing (each step is one week):
#   inf.prob        0.11   weekly transmission probability per partnership
#                          (one act per week at a per-act probability in the
#                          published range), tuned to about 3% prevalence
#   rec.rate        1/70   untreated infections clear after 70 weeks on average
#   symp.prob.f/m   0.1, 0.2   share of infections that are symptomatic
#   symp.test.rate  0.25   symptomatic infections come to care after a mean of
#                          4 weeks
#   screen.rate.f/m        weekly screening probabilities equal to 40% of women
#                          and 8% of men tested per year
#   tx.prob         0.95   probability that treatment cures
#   reinf.window    13     repeat infections are counted within 13 weeks
#
# Partner services (the burn-in runs patient referral, the US standard of care):
#   pn.arm          "none", "PR" (patient referral), or "EPT" (expedited
#                   partner therapy)
#   pn.lookback     partners from this many weeks before diagnosis (9 weeks
#                   is the CDC 60-day window)
#   reach.main      probability of reaching an ongoing main partner
#   reach.cas       probability of reaching an ongoing casual partner
#   reach.ended     probability of reaching a partner whose partnership ended
#                   within the window
param_base <- param.net(
  inf.prob = 0.11,
  rec.rate = 1 / 70,
  symp.prob.f = 0.1,
  symp.prob.m = 0.2,
  symp.test.rate = 0.25,
  screen.rate.f = 1 - (1 - 0.4)^(1 / 52),
  screen.rate.m = 1 - (1 - 0.08)^(1 / 52),
  tx.prob = 0.95,
  reinf.window = 13,
  pn.arm = "PR",
  pn.lookback = 9,
  reach.main = 0.45,
  reach.cas = 0.20,
  reach.ended = 0.10,
  departure.rate = departure_rate,
  arrival.rate = departure_rate
)

# 4% of the population starts infected; the burn-in forgets this.
init <- init.net(i.num = round(0.04 * N))

# The same controls serve the burn-in and the scenarios, which differ only in
# the steps they run and in what is saved:
#   save.run = TRUE keeps the full state of every burn-in simulation
#   (attributes, networks, cumulative edgelist) so it can be resumed, and
#   save.cumulative.edgelist = TRUE returns the partnership history.
# truncate.el.cuml = 52 keeps partnerships that ended within the last year:
# the longest lookback is 26 weeks, and the analysis counts partners over 52.
# attr.rules gives arrivals a sex and risk group drawn from the proportions
# at the start (step 1), and set.control.tergm passes the longer Markov chain
# used in the diagnostics above.
make_control <- function(start, nsteps, save) {
  control.net(
    type = NULL,
    nsims = nsims,
    ncores = ncores,
    start = start,
    nsteps = nsteps,
    tergmLite = TRUE,
    resimulate.network = TRUE,
    cumulative.edgelist = TRUE,
    truncate.el.cuml = 52,
    initAttr.FUN = init_attrs,
    infection.FUN = infect,
    test.FUN = test_treat,
    notify.FUN = notify,
    recovery.FUN = clear,
    departures.FUN = depart,
    arrivals.FUN = arrive,
    module.order = c("resim_nets.FUN", "summary_nets.FUN", "initAttr.FUN",
                     "infection.FUN", "test.FUN", "notify.FUN",
                     "recovery.FUN", "departures.FUN", "arrivals.FUN",
                     "nwupdate.FUN", "prevalence.FUN"),
    attr.rules = list(sex = "t1", risk = "t1"),
    epi.by = "sex",
    set.control.tergm = tergm_ctrl,
    save.run = save,
    save.cumulative.edgelist = save,
    verbose = FALSE
  )
}


# 3. Burn-in to Endemic Equilibrium ----------------------------------------

set.seed(2027)
sim_burnin <- netsim(list(est_main, est_cas), param_base, init,
                     make_control(start = 1, nsteps = burnin, save = TRUE))

# as.data.frame() numbers the steps of a resumed simulation from its start
# step although the output also holds the burn-in history (EpiModel 2.6.2),
# so steps are numbered here by position, which is correct either way.
epi_df <- function(sim) {
  df <- as.data.frame(sim)
  df$time <- ave(df$time, df$sim, FUN = seq_along)
  df
}

# Calibration checks over the last two years of the burn-in, pooled over
# simulations
bdf <- epi_df(sim_burnin)
bdf <- bdf[bdf$time > burnin - 104, ]
tot <- function(x) sum(x, na.rm = TRUE)

# Partners in the past year for people present all year, from the saved
# cumulative edgelist (one data frame per simulation, both layers, unique ids)
partners_past_year <- sapply(seq_len(nsims), function(s) {
  el <- sim_burnin$cumulative.edgelist[[s]]
  el <- el[is.na(el$stop) | el$stop > burnin - 52, ]
  attr <- sim_burnin$run[[s]]$attr
  ids <- attr$unique_id[attr$entrTime <= burnin - 51]
  n <- table(factor(c(el$head, el$tail), levels = ids))
  c(mean = mean(n), two_plus = 100 * mean(n >= 2))
})

found <- tot(bdf$pn.found.main.flow + bdf$pn.found.cas.flow +
               bdf$pn.found.ended.flow)
reached_cur <- tot(bdf$pn.reach.main.flow + bdf$pn.reach.cas.flow)
dx_f <- tot(bdf$dx.flow.f + bdf$pn.dx.flow.f)
dx_m <- tot(bdf$dx.flow + bdf$pn.dx.flow) - dx_f
calib_tbl <- data.frame(
  Quantity = c("Prevalence, women (%)", "Prevalence, men (%)",
               "New infections diagnosed, women (%)",
               "New infections diagnosed, men (%)",
               "Diagnoses, women per man",
               "Partners in the past year, mean",
               "Two or more partners in the past year (%)",
               "Partners in the 60-day window per index",
               "Reached ongoing partners who are infected (%)",
               "Repeat infection within 13 weeks of cure (%)"),
  Model = sprintf(c("%.1f", "%.1f", "%.0f", "%.0f", "%.1f", "%.2f", "%.0f",
                    "%.2f", "%.0f", "%.0f"),
                  c(100 * mean(bdf$i.num.sexF / bdf$num.sexF),
                    100 * mean(bdf$i.num.sexM / bdf$num.sexM),
                    100 * dx_f / tot(bdf$si.flow.f),
                    100 * dx_m / tot(bdf$si.flow - bdf$si.flow.f),
                    dx_f / dx_m,
                    mean(partners_past_year["mean", ]),
                    mean(partners_past_year["two_plus", ]),
                    found / tot(bdf$pn.index.flow),
                    100 * tot(bdf$pn.inf.main.flow + bdf$pn.inf.cas.flow) /
                      reached_cur,
                    100 * tot(bdf$reinf.flow) / tot(bdf$tx.flow))),
  Target = c("4.7", "2.7", "35 to 50", "25 to 40", "about 2.5", "1.70",
             "20 to 35", "1.18 to 1.48 named in 6 months", "50 to 75",
             "10 to 17")
)
print(calib_tbl)


# 4. Scenarios Resumed from the Burn-in ------------------------------------

# Every scenario resumes each of the burn-in simulations at week burnin + 1
# and runs for another five years, so the scenarios start from identical
# endemic states and are compared within simulation. Reach probabilities for
# EPT are 20, 15, and 5 percentage points above patient referral.
#   none     no partner services
#   pr       patient referral, 60-day lookback (continues the burn-in)
#   pr_6mo   patient referral, 6-month lookback (UK guidance)
#   ept      expedited partner therapy, 60-day lookback
#   ept_cur  expedited partner therapy for ongoing partners only
scenarios.df <- data.frame(
  .scenario.id = c("none", "pr", "pr_6mo", "ept", "ept_cur"),
  .at          = 0,
  pn.arm       = c("none", "PR", "PR", "EPT", "EPT"),
  pn.lookback  = c(9, 9, 26, 9, 0),
  reach.main   = c(0.45, 0.45, 0.45, 0.65, 0.65),
  reach.cas    = c(0.20, 0.20, 0.20, 0.35, 0.35),
  reach.ended  = c(0.10, 0.10, 0.10, 0.15, 0.15)
)
scenarios.list <- create_scenario_list(scenarios.df)

labels <- c(none = "No partner services",
            pr = "Patient referral, 60 days",
            pr_6mo = "Patient referral, 6 months",
            ept = "EPT, 60 days",
            ept_cur = "EPT, ongoing partners only")
cols <- c(none = "gray40", pr = "goldenrod", pr_6mo = "mediumpurple",
          ept = "seagreen", ept_cur = "steelblue")

control_scn <- make_control(start = burnin + 1, nsteps = burnin + horizon,
                            save = FALSE)
sims <- list()
for (scn in scenarios.list) {
  cat("\n--- Running scenario:", scn$id, "---\n")
  set.seed(2028)
  sims[[scn$id]] <- netsim(sim_burnin, use_scenario(param_base, scn), init,
                           control_scn)
}


# 5. Analysis --------------------------------------------------------------

# Per-simulation totals over the five scenario years, one value per
# simulation, so that the Monte Carlo intervals below can be built from the
# within-simulation differences between scenarios.
types <- c("main", "cas", "ended")
summarize_scenario <- function(sim) {
  df <- epi_df(sim)
  df <- df[df$time > burnin, ]
  by_sim <- function(x) as.numeric(tapply(x, df$sim, sum, na.rm = TRUE))
  by_type <- function(prefix) {
    sapply(types, function(k) by_sim(df[[paste0(prefix, k, ".flow")]]))
  }
  end_prev <- function(i, n) {
    100 * as.numeric(tapply(i / n, df$sim, function(v) mean(tail(v, 52))))
  }
  list(inf = by_sim(df$si.flow),
       prev = end_prev(df$i.num, df$num),
       tx = by_sim(df$tx.flow),
       index = by_sim(df$pn.index.flow),
       found = matrix(by_type("pn.found."), ncol = 3),
       reach = matrix(by_type("pn.reach."), ncol = 3),
       reach_inf = matrix(by_type("pn.inf."), ncol = 3),
       part_dx = by_sim(df$pn.dx.flow),
       reinf = by_sim(df$reinf.flow),
       reinf_prior = by_sim(df$reinf.prior.flow))
}
res <- lapply(sims, summarize_scenario)

# Monte Carlo intervals, with a t quantile because there are only ten
# simulations. The scenarios resume the same burn-in simulations, so a
# difference between two scenarios is taken within each simulation before it
# is averaged.
mc_mean <- function(x) {
  m <- mean(x)
  if (length(x) < 2) return(c(est = m, lo = NA, hi = NA))
  half <- qt(0.975, length(x) - 1) * sd(x) / sqrt(length(x))
  c(est = m, lo = m - half, hi = m + half)
}
mc_paired <- function(x0, x1) mc_mean(x0 - x1)
fmt_ci <- function(v, digits = 1) {
  if (is.na(v["lo"])) return(sprintf("%.*f", digits, v["est"]))
  sprintf("%.*f (%.*f, %.*f)", digits, v["est"], digits, v["lo"], digits, v["hi"])
}


## --- Plot 1: prevalence over the last three burn-in years and five scenario years

par(mfrow = c(1, 1), mar = c(4, 4.2, 2.5, 1), mgp = c(2.4, 0.8, 0))
prev_t <- lapply(sims, function(sim) {
  df <- epi_df(sim)
  100 * tapply(df$i.num / df$num, df$time, mean)
})
show <- max(1, burnin - 155):(burnin + horizon)
plot(NA, xlim = range(show - burnin) / 52,
     ylim = c(0, max(sapply(prev_t, max)) * 1.1),
     xlab = "Years since the change in partner services",
     ylab = "Prevalence (%, mean of simulations)",
     main = "Prevalence by Partner Services Strategy")
abline(v = 0, lty = 2, col = "gray50")
for (s in names(sims)) {
  lines((show - burnin) / 52, prev_t[[s]][show], lwd = 2, col = cols[s])
}
legend("topleft", legend = labels, col = cols, lwd = 2, bty = "n", cex = 0.75)


## --- Infections averted relative to patient referral ---------------------

alt <- c("none", "pr_6mo", "ept", "ept_cur")
averted <- lapply(alt, function(s) mc_paired(res$pr$inf, res[[s]]$inf))
names(averted) <- alt
added_none <- mc_paired(res$none$inf, res$pr$inf)
ept_window <- mc_paired(res$ept_cur$inf, res$ept$inf)
treated <- sapply(res, function(r) mean(rowSums(r$reach)))
out_tbl <- data.frame(
  Scenario = labels[names(res)],
  `Prevalence after 5 years (%)` = sapply(res, function(r) fmt_ci(mc_mean(r$prev))),
  `New infections` = sapply(res, function(r) round(mean(r$inf))),
  `Infections averted vs patient referral` =
    c(fmt_ci(averted$none, 0), "", sapply(averted[-1], fmt_ci, digits = 0)),
  `Percent change in infections` = sapply(res, function(r)
    round(100 * (mean(r$inf) - mean(res$pr$inf)) / mean(res$pr$inf), 1)),
  `Partners treated` = round(treated),
  `Infections averted per partner treated, vs none` =
    c("", sapply(names(res)[-1], function(s)
      sprintf("%.1f", (mean(res$none$inf) - mean(res[[s]]$inf)) / treated[[s]]))),
  check.names = FALSE, row.names = NULL
)
print(out_tbl)
cat("EPT, 60 days vs ongoing partners only:", fmt_ci(ept_window, 0), "infections\n")


## --- Plot 2: infections averted relative to patient referral ------------

par(mfrow = c(1, 1), mar = c(4, 4.5, 3, 1), mgp = c(2.6, 0.8, 0))
short <- c(none = "No partner services", pr_6mo = "PR, 6 months",
           ept = "EPT, 60 days", ept_cur = "EPT, ongoing only")
av <- sapply(averted, function(v) 1000 * v / N)
bp <- barplot(av["est", ], names.arg = short[alt], col = cols[alt],
              ylab = "Infections averted per 1,000 people",
              main = "Five-Year Infections Averted Relative to Patient Referral",
              ylim = range(c(0, av), na.rm = TRUE) * 1.15, cex.names = 0.85)
abline(h = 0)
if (!all(is.na(av["lo", ]))) {
  arrows(bp, av["lo", ], bp, av["hi", ], angle = 90, code = 3, length = 0.04)
}


## --- Repeat infection of cured cases -------------------------------------

reinf_tbl <- data.frame(
  Scenario = labels[names(res)],
  `Cured cases` = sapply(res, function(r) round(mean(r$tx))),
  `Repeat infection within 13 weeks (%)` = sapply(res, function(r)
    round(100 * sum(r$reinf) / sum(r$tx), 1)),
  `From a partner they already had (%)` = sapply(res, function(r)
    round(100 * sum(r$reinf_prior) / sum(r$reinf))),
  `Relative to patient referral` = sapply(res, function(r)
    round((sum(r$reinf) / sum(r$tx)) / (sum(res$pr$reinf) / sum(res$pr$tx)), 2)),
  check.names = FALSE, row.names = NULL
)
print(reinf_tbl)

par(mfrow = c(1, 1), mar = c(4, 4.5, 3, 1), mgp = c(2.6, 0.8, 0))
reinf <- sapply(res, function(r) {
  100 * c(prior = sum(r$reinf_prior), new = sum(r$reinf - r$reinf_prior)) / sum(r$tx)
})
barplot(reinf, names.arg = c("None", "PR, 60 days", "PR, 6 months",
                             "EPT, 60 days", "EPT, ongoing"),
        col = c("firebrick", "gray75"), cex.names = 0.85,
        ylab = "Reinfected within 13 weeks of cure (%)",
        main = "Repeat Infection of Cured Cases",
        ylim = c(0, max(colSums(reinf)) * 1.3))
legend("topright", fill = c("firebrick", "gray75"), bty = "n", cex = 0.8,
       legend = c("By a partner they had when cured", "By a new partner"))


## --- What partner services finds ----------------------------------------

pn_scn <- c("pr", "pr_6mo", "ept", "ept_cur")
per_index <- function(s, x) sum(x) / sum(res[[s]]$index)
yield_tbl <- data.frame(
  Scenario = labels[pn_scn],
  `Indices per year` = sapply(pn_scn, function(s)
    round(mean(res[[s]]$index) / (horizon / 52))),
  `Partners found per index` = sapply(pn_scn, function(s)
    round(per_index(s, res[[s]]$found), 2)),
  `Ongoing main, casual, ended (%)` = sapply(pn_scn, function(s) {
    paste(round(100 * colSums(res[[s]]$found) / sum(res[[s]]$found)),
          collapse = ", ")
  }),
  `Reached per index` = sapply(pn_scn, function(s)
    round(per_index(s, res[[s]]$reach), 2)),
  `Reached who are infected (%)` = sapply(pn_scn, function(s)
    round(100 * sum(res[[s]]$reach_inf) / sum(res[[s]]$reach))),
  `Infected by type: main, casual, ended (%)` = sapply(pn_scn, function(s) {
    p <- 100 * colSums(res[[s]]$reach_inf) / colSums(res[[s]]$reach)
    paste(ifelse(is.nan(p), "-", round(p)), collapse = ", ")
  }),
  `Partners diagnosed per index` = sapply(pn_scn, function(s)
    round(per_index(s, res[[s]]$part_dx), 2)),
  check.names = FALSE, row.names = NULL
)
print(yield_tbl)

par(mfrow = c(1, 1), mar = c(4, 4.5, 3, 1), mgp = c(2.6, 0.8, 0))
yield <- sapply(pn_scn, function(s) {
  r <- res[[s]]
  c(colSums(r$reach_inf), colSums(r$reach) - colSums(r$reach_inf)) / sum(r$index)
})
ycols <- c("#8e44ad", "#e67e22", "#566573", "#d7bde2", "#f5cba7", "#d5d8dc")
barplot(yield, names.arg = c("PR, 60 days", "PR, 6 months", "EPT, 60 days",
                             "EPT, ongoing only"),
        col = ycols, cex.names = 0.85, ylab = "Partners reached per index",
        main = "Who Partner Services Reaches",
        ylim = c(0, max(colSums(yield)) * 1.45))
legend("top", ncol = 3, bty = "n", fill = ycols[c(1, 4, 2, 5, 3, 6)], cex = 0.75,
       legend = c("Main, infected", "Main, uninfected", "Casual, infected",
                  "Casual, uninfected", "Ended, infected", "Ended, uninfected"))


## --- Between-simulation variability ---------------------------------------

cv_tbl <- data.frame(
  Scenario = labels[names(res)],
  `Mean prevalence after 5 years (%)` = sapply(res, function(r) round(mean(r$prev), 1)),
  `Range across simulations (%)` = sapply(res, function(r)
    sprintf("%.1f to %.1f", min(r$prev), max(r$prev))),
  `CV of five-year infections` = sapply(res, function(r)
    round(sd(r$inf) / mean(r$inf), 2)),
  check.names = FALSE, row.names = NULL
)
print(cv_tbl)
