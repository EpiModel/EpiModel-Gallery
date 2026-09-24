##
## SEIR with Contact Tracing for an Acute, Immunizing Infection
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
#   Rscript examples/seir-contact-tracing/model.R "run_full <- TRUE"
# The full run takes two to three minutes on five cores; CI mode runs in
# well under a minute and its results are not meant to be interpreted.
if (interactive() || exists("run_full")) {
  N <- 5000
  nsims <- 10
  ncores <- 5
  nsteps <- 250
} else {
  N <- 1000
  nsims <- 1
  ncores <- 1
  nsteps <- 50
}


# 1. Network Model Estimation ----------------------------------------------

# A single dynamic network of close contacts. Each node has 6 contacts on
# average and each contact relationship lasts a week on average, so over an
# infectious period of about 8 days a case has roughly 12 distinct close
# contacts and about half of the contacts that a tracer looks for have
# already ended. Each time step is one day. The seed makes the whole run,
# including the network fit, reproducible.
set.seed(2026)
mean_degree <- 6
duration <- 7

nw <- network_initialize(N)
formation <- ~edges
target.stats <- round(mean_degree * N / 2)
coef.diss <- dissolution_coefs(~offset(edges), duration = duration)
est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)

dx <- netdx(est, nsims = 5, ncores = ncores, nsteps = nsteps,
            nwstats.formula = ~edges + degree(0:4, by = NULL) + meandeg,
            verbose = FALSE)
print(dx)
if (interactive()) plot(dx)


# 2. Parameters, Initial Conditions, and Controls --------------------------

if (file.exists("examples/seir-contact-tracing/module-fx.R")) {
  source("examples/seir-contact-tracing/module-fx.R")
} else {
  source("module-fx.R")
}

# Natural history (each step is one day):
#   ei.rate      1/3    mean latent period 3 days; with the 2.5-day
#                       presymptomatic stage this gives a mean incubation
#                       period of 5.5 days
#   ips.rate     0.4    mean presymptomatic infectious period 2.5 days
#   isr.rate     1/6    mean symptomatic infectious period 6 days
#   iar.rate     1/8    mean asymptomatic infectious period 8 days
#   asymp.prob   0.3    share of infections that never develop symptoms
#   is.inf.mult  0.5    per-contact infectiousness in the symptomatic stage
#                       relative to the presymptomatic stage, so that about
#                       45% of transmission from symptomatic infections
#                       occurs before symptom onset
#   ia.inf.mult  0.35   infectiousness of asymptomatic infections relative
#                       to the presymptomatic stage
#   inf.prob     0.065  per-contact daily transmission probability in the
#                       presymptomatic stage, tuned so that the reproduction
#                       number without interventions is about 1.7 on this
#                       network (a partially mitigated epidemic)
#
# Case-based interventions:
#   dx.prob      0.5    share of symptomatic cases that are ever diagnosed
#   dx.delay     3      mean days from symptom onset to diagnosis
#   iso.duration 10     days of isolation after diagnosis
#   iso.mult     0.2    contact multiplier for an isolated index (80% cut)
#   trace.reach.prob    share of identified contacts reached and quarantined
#   trace.delay         days from the index's diagnosis to contact reach
#   trace.window 2      contacts are elicited from this many days before the
#                       index's symptom onset up to the day of diagnosis
#   quar.duration 10    days of quarantine for a reached contact
#   quar.mult    0.3    contact multiplier for a quarantined contact (70% cut)
#
# Tracing is off in the base parameter set; the scenarios turn it on.
param_base <- param.net(
  inf.prob = 0.065,
  is.inf.mult = 0.5,
  ia.inf.mult = 0.35,
  ei.rate = 1 / 3,
  ips.rate = 0.4,
  isr.rate = 1 / 6,
  iar.rate = 1 / 8,
  asymp.prob = 0.3,
  dx.prob = 0.5,
  dx.delay = 3,
  iso.duration = 10,
  iso.mult = 0.2,
  trace.reach.prob = 0,
  trace.delay = 1,
  trace.window = 2,
  quar.duration = 10,
  quar.mult = 0.3
)

# 0.5% of the population starts infectious, so that the epidemic takes off
# without a long stochastic lag and few simulations die out.
init <- init.net(i.num = round(0.005 * N))

# The cumulative edgelist is switched on here. truncate.el.cuml drops
# partnerships that ended more than 30 days ago, longer than the elicitation
# window the trace module asks for (2 days before onset, the diagnosis
# delay, and the tracing delay) for all but a vanishing share of indices.
# Its default of 0 does not mean "keep everything": it means dissolved
# partnerships are never recorded, so a tracer could only find current
# partners. Module order is set explicitly so that the cumulative edgelist
# is updated (in resim_nets) before infection and tracing read the network,
# and so that tracing runs after the day's diagnoses.
control <- control.net(
  type = NULL,
  nsims = nsims,
  ncores = ncores,
  nsteps = nsteps,
  tergmLite = TRUE,
  resimulate.network = TRUE,
  cumulative.edgelist = TRUE,
  truncate.el.cuml = 30,
  initAttr.FUN = init_attrs,
  infection.FUN = infect,
  progress.FUN = progress,
  trace.FUN = trace,
  module.order = c("resim_nets.FUN", "summary_nets.FUN", "initAttr.FUN",
                   "infection.FUN", "progress.FUN", "trace.FUN",
                   "nwupdate.FUN", "prevalence.FUN"),
  verbose = FALSE
)


# 3. Scenarios -------------------------------------------------------------

# Five scenarios on the same network, natural history, and seeds:
#   none       no case-based intervention
#   iso        symptom-based diagnosis with isolation of diagnosed cases
#   fast_high  isolation plus tracing: contacts reached 1 day after the
#              index's diagnosis, 80% of contacts reached
#   slow_high  contacts reached 4 days after diagnosis, 80% reached
#   fast_low   contacts reached 1 day after diagnosis, 30% reached
scenarios.df <- data.frame(
  .scenario.id     = c("none", "iso", "fast_high", "slow_high", "fast_low"),
  .at              = 0,
  dx.prob          = c(0,   0.5, 0.5, 0.5, 0.5),
  trace.reach.prob = c(0,   0,   0.8, 0.8, 0.3),
  trace.delay      = c(1,   1,   1,   4,   1)
)
scenarios.list <- create_scenario_list(scenarios.df)

labels <- c(none = "No intervention",
            iso = "Isolation only",
            fast_high = "Tracing: fast (1 d), 80% reached",
            slow_high = "Tracing: slow (4 d), 80% reached",
            fast_low = "Tracing: fast (1 d), 30% reached")
cols <- c(none = "gray40", iso = "goldenrod", fast_high = "seagreen",
          slow_high = "firebrick", fast_low = "steelblue")

sims <- list()
for (scn in scenarios.list) {
  cat("\n--- Running scenario:", scn$id, "---\n")
  sims[[scn$id]] <- netsim(est, use_scenario(param_base, scn), init, control)
}
print(sims$fast_high)


# 4. Analysis --------------------------------------------------------------

# Per-simulation outcomes, kept one row per simulation so that the Monte
# Carlo intervals below can be built from the between-simulation variance.
# The counters are NA on step 1, before any module has run, so that step is
# dropped; seeds carry infTime = 1 and never appear as incident infections.
summarize_scenario <- function(sim) {
  df <- as.data.frame(sim)
  df <- df[df$time > 1, ]                          # counters are NA on step 1
  by_sim <- function(x) as.numeric(tapply(x, df$sim, sum, na.rm = TRUE))
  n_sim <- length(unique(df$sim))

  cum_inf <- by_sim(df$se.flow)
  symp <- by_sim(df$ips.flow)
  dxs <- by_sim(df$dx.flow)
  index <- by_sim(df$trace.index.flow)
  part <- by_sim(df$trace.part.flow)
  part_ended <- by_sim(df$trace.part.ended.flow)
  reach <- by_sim(df$trace.reach.flow)
  quar_start <- by_sim(df$quar.start.flow)
  quar_days <- by_sim(df$quar.num)
  iso_days <- by_sim(df$iso.num)
  reach_state <- cbind(s = by_sim(df$reach.s.flow), e = by_sim(df$reach.e.flow),
                       i = by_sim(df$reach.i.flow), r = by_sim(df$reach.r.flow))
  quar_pct <- 100 * tapply(df$quar.num, df$time, mean) / N
  active_end <- with(df[df$time == max(df$time), ],
                     mean(e.num + ip.num + is.num + ia.num))

  # Transmission record, pooled over simulations
  tm <- do.call(rbind, lapply(seq_len(n_sim), function(k) {
    as.data.frame(get_transmat(sim, sim = k))
  }))
  r_seed <- sum(tm$infTime == 1) / (init$i.num * n_sim)
  gen_time <- mean(tm$at - tm$infTime)
  stage_share <- prop.table(table(factor(tm$infStage, levels = c("ip", "is", "ia"))))
  restricted_share <- mean(tm$infIsolated == 1 | tm$anyQuarantined == 1)

  list(cum_inf = cum_inf, attack = 100 * cum_inf / N,
       symp = symp, dx = dxs, index = index, part = part,
       part_ended = part_ended, reach = reach, quar_start = quar_start,
       quar_days = quar_days, iso_days = iso_days, reach_state = reach_state,
       quar_pct = quar_pct, active_end = active_end,
       r_seed = r_seed, gen_time = gen_time, stage_share = stage_share,
       restricted_share = restricted_share)
}

res <- lapply(sims, summarize_scenario)

# Monte Carlo interval helpers, as in the other Gallery examples
mc_mean <- function(x) {
  m <- mean(x)
  if (length(x) < 2) return(c(est = m, lo = NA, hi = NA))
  se <- sd(x) / sqrt(length(x))
  c(est = m, lo = m - 1.96 * se, hi = m + 1.96 * se)
}
mc_diff <- function(x0, x1) {
  d <- mean(x0) - mean(x1)
  if (length(x0) < 2 || length(x1) < 2) return(c(est = d, lo = NA, hi = NA))
  se <- sqrt(var(x0) / length(x0) + var(x1) / length(x1))
  c(est = d, lo = d - 1.96 * se, hi = d + 1.96 * se)
}
fmt_ci <- function(v, digits = 1) {
  if (is.na(v["lo"])) return(sprintf("%.*f", digits, v["est"]))
  sprintf("%.*f (%.*f, %.*f)", digits, v["est"], digits, v["lo"], digits, v["hi"])
}


## --- Epidemic size and the natural history check -------------------------

epi_tbl <- data.frame(
  Scenario = labels[names(res)],
  `Attack rate (%)` = sapply(res, function(r) fmt_ci(mc_mean(r$attack))),
  `Range` = sapply(res, function(r) sprintf("%.1f to %.1f", min(r$attack), max(r$attack))),
  `Still infected at end` = sapply(res, function(r) round(r$active_end)),
  `Symptomatic cases diagnosed (%)` = sapply(res, function(r)
    ifelse(sum(r$symp) > 0, round(100 * sum(r$dx) / sum(r$symp)), 0)),
  check.names = FALSE, row.names = NULL
)
print(epi_tbl)

# Where transmission comes from, by scenario: the share of transmissions
# from each substage of the infector and the share that occurred while the
# infector was isolated or either partner was quarantined.
source_tbl <- t(sapply(res, function(r) {
  c(round(100 * as.numeric(r$stage_share)),
    round(100 * r$restricted_share),
    round(r$r_seed, 2), round(r$gen_time, 1))
}))
colnames(source_tbl) <- c("Presymptomatic (%)", "Symptomatic (%)",
                          "Asymptomatic (%)", "Under isolation or quarantine (%)",
                          "Secondary infections per seed", "Generation time (days)")
print(source_tbl)


## --- Plot 1: daily and cumulative incidence -------------------------------

smooth_ma <- function(x, k = 7) {
  out <- as.numeric(stats::filter(x, rep(1 / k, k), sides = 2))
  names(out) <- names(x)
  out
}
inc <- lapply(sims, function(sim) {
  df <- as.data.frame(sim)
  df$se.flow[is.na(df$se.flow)] <- 0
  tapply(df$se.flow, df$time, mean)
})
cum <- lapply(inc, cumsum)

par(mfrow = c(1, 2), mar = c(4, 4.2, 2.5, 1), mgp = c(2.4, 0.8, 0))
ymax <- max(sapply(inc, function(v) max(smooth_ma(v), na.rm = TRUE)))
plot(NA, xlim = c(1, nsteps), ylim = c(0, ymax * 1.05),
     xlab = "Day", ylab = "New infections per day (mean, 7-day smoothed)",
     main = "Daily Incidence")
for (s in names(sims)) {
  lines(as.numeric(names(inc[[s]])), smooth_ma(inc[[s]]), lwd = 2, col = cols[s])
}
legend("topright", legend = labels, col = cols, lwd = 2, bty = "n", cex = 0.75)
plot(NA, xlim = c(1, nsteps), ylim = c(0, max(sapply(cum, max)) / N * 105),
     xlab = "Day", ylab = "Cumulative attack rate (%)",
     main = "Cumulative Incidence")
for (s in names(sims)) {
  lines(as.numeric(names(cum[[s]])), 100 * cum[[s]] / N, lwd = 2, col = cols[s])
}


## --- Infections averted and the cost of averting them --------------------

# Tracing scenarios are compared with isolation only, the standard of care
# they add to, and every difference carries a 95% Monte Carlo interval.
trace_scn <- c("fast_high", "slow_high", "fast_low")
averted <- lapply(trace_scn, function(s) mc_diff(res$iso$cum_inf, res[[s]]$cum_inf))
names(averted) <- trace_scn
speed_diff <- mc_diff(res$slow_high$cum_inf, res$fast_high$cum_inf)
coverage_diff <- mc_diff(res$fast_low$cum_inf, res$fast_high$cum_inf)

int_tbl <- data.frame(
  Scenario = labels[trace_scn],
  `Infections averted` = sapply(averted, fmt_ci, digits = 0),
  `Percent of isolation-only infections` = sapply(trace_scn, function(s)
    round(100 * averted[[s]]["est"] / mean(res$iso$cum_inf), 1)),
  `Contacts per index` = sapply(trace_scn, function(s)
    round(sum(res[[s]]$part) / sum(res[[s]]$index), 1)),
  `Ended partnerships (%)` = sapply(trace_scn, function(s)
    round(100 * sum(res[[s]]$part_ended) / sum(res[[s]]$part))),
  `Quarantine episodes` = sapply(trace_scn, function(s) round(mean(res[[s]]$quar_start))),
  `Quarantine days per infection averted` = sapply(trace_scn, function(s)
    round(mean(res[[s]]$quar_days) / averted[[s]]["est"])),
  `Peak share of population in quarantine (%)` = sapply(trace_scn, function(s)
    round(max(res[[s]]$quar_pct), 1)),
  check.names = FALSE, row.names = NULL
)
print(int_tbl)
cat("Fast vs slow at 80% reach:", fmt_ci(speed_diff, 0), "infections\n")
cat("80% vs 30% reach at 1 day:", fmt_ci(coverage_diff, 0), "infections\n")


## --- Plot 2: people under isolation or quarantine over time --------------

par(mfrow = c(1, 1), mar = c(4, 4.2, 2.5, 1), mgp = c(2.4, 0.8, 0))
restricted <- lapply(sims, function(sim) {
  df <- as.data.frame(sim)
  100 * tapply(df$iso.num + df$quar.num, df$time, mean) / N
})
plot(NA, xlim = c(1, nsteps), ylim = c(0, max(sapply(restricted, max, na.rm = TRUE)) * 1.1),
     xlab = "Day", ylab = "Percent of population isolated or quarantined",
     main = "Population Under Movement Restriction")
for (s in names(sims)) {
  lines(as.numeric(names(restricted[[s]])), restricted[[s]], lwd = 2, col = cols[s])
}
legend("topright", legend = labels, col = cols, lwd = 2, bty = "n", cex = 0.75)


## --- Plot 3: infections averted and quarantine per infection averted -----

par(mfrow = c(1, 2), mar = c(8, 4.5, 3, 1), mgp = c(3, 0.8, 0))
short <- c(fast_high = "Fast, 80%", slow_high = "Slow, 80%", fast_low = "Fast, 30%")
av_est <- sapply(averted, function(v) 100 * v["est"] / N)
av_lo <- sapply(averted, function(v) 100 * v["lo"] / N)
av_hi <- sapply(averted, function(v) 100 * v["hi"] / N)
bp <- barplot(av_est, names.arg = short[trace_scn], col = cols[trace_scn], las = 2,
              ylab = "Infections averted per 100 population",
              main = "Averted vs Isolation Only",
              ylim = c(min(0, av_lo, av_est, na.rm = TRUE),
                       max(0, av_est, av_hi, na.rm = TRUE) * 1.2))
if (!all(is.na(av_lo))) arrows(bp, av_lo, bp, av_hi, angle = 90, code = 3, length = 0.04)
text(bp, pmax(av_hi, av_est, na.rm = TRUE) + max(c(av_est, av_hi), na.rm = TRUE) * 0.05,
     sprintf("%.1f", av_est), cex = 0.85, font = 2)
qd <- sapply(trace_scn, function(s) mean(res[[s]]$quar_days) / averted[[s]]["est"])
bp2 <- barplot(qd, names.arg = short[trace_scn], col = cols[trace_scn], las = 2,
               ylab = "Quarantine person-days per infection averted",
               main = "Cost of Averting One Infection",
               ylim = c(min(0, qd), max(0, qd) * 1.2))
text(bp2, qd + max(qd) * 0.05, sprintf("%.0f", qd), cex = 0.85, font = 2)


## --- Tracing yield: what state were the reached contacts in? --------------

yield <- sapply(trace_scn, function(s) {
  m <- colSums(res[[s]]$reach_state)
  100 * m / sum(m)
})
yield_tbl <- data.frame(
  Scenario = labels[trace_scn],
  `Contacts reached` = sapply(trace_scn, function(s) round(mean(res[[s]]$reach))),
  `Susceptible (%)` = round(yield["s", ], 1),
  `Latent (%)` = round(yield["e", ], 1),
  `Infectious (%)` = round(yield["i", ], 1),
  `Recovered (%)` = round(yield["r", ], 1),
  check.names = FALSE, row.names = NULL
)
print(yield_tbl)

par(mfrow = c(1, 1), mar = c(4, 4.5, 3, 1), mgp = c(2.6, 0.8, 0))
ycols <- c(s = "#3498db", e = "#8e44ad", i = "#e74c3c", r = "#27ae60")
bp3 <- barplot(yield, names.arg = short[trace_scn], col = ycols, las = 1,
               ylab = "Percent of reached contacts", ylim = c(0, 118),
               main = "State of Contacts When Reached")
legend("top", horiz = TRUE, bty = "n", fill = ycols, cex = 0.85,
       legend = c("Susceptible", "Latent", "Infectious", "Recovered"))
text(bp3, 100 - yield["s", ] / 2, sprintf("%.0f%%", yield["s", ]), col = "white")


## --- Between-simulation variability ---------------------------------------

cv_tbl <- t(sapply(res, function(r) {
  c(`Mean attack rate (%)` = round(mean(r$attack), 1),
    CV = round(sd(r$attack) / mean(r$attack), 2))
}))
print(cv_tbl)
