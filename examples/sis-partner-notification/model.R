
##
## Partner Notification for an Endemic STI (SIS)
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##
## Author: Samuel M. Jenness (Emory University)
## Date: May 2026
##

# Load EpiModel
suppressMessages(library(EpiModel))

# Standard Gallery unit test lines
rm(list = ls())
eval(parse(text = print(commandArgs(TRUE)[1])))

if (interactive()) {
  nsims     <- 5
  ncores    <- 5
  nsteps    <- 600
  n         <- 1000
  pn.start  <- 300
} else {
  nsims     <- 2
  ncores    <- 2
  nsteps    <- 80
  n         <- 250
  pn.start  <- 40
}


# 1. Network Model Estimation ------------------------------------------------

# A sparse sexual contact network with some concurrency. The network
# model is the same across scenarios so the only thing that varies in the
# comparison is the partner-notification configuration.
nw <- network_initialize(n)

# Mean degree 0.7 with ~7% concurrent; partnerships average 100 weeks.
# These are in the empirically observed range for sexual partnership
# networks.
formation    <- ~edges + concurrent
target.stats <- c(round(0.7 * n / 2), round(0.07 * n))
coef.diss    <- dissolution_coefs(~offset(edges), duration = 100)

est <- netest(nw, formation, target.stats, coef.diss)

dx <- netdx(est, nsims = nsims, ncores = ncores, nsteps = nsteps,
            nwstats.formula = ~edges + concurrent + degree(0:3))
print(dx)
plot(dx)


# 2. Epidemic Model Parameters -----------------------------------------------

source("examples/sis-partner-notification/module-fx.R")

# Disease parameters frame an endemic, chlamydia-like STI:
#
#   inf.prob   = 0.35    per-act transmission probability (elevated to
#                        sustain endemic dynamics on the sparse network;
#                        illustrative for teaching, not calibrated to data)
#   act.rate   = 1       1 act per partnership per week
#   rec.rate   = 0.01    slow natural clearance (mean ~100 weeks)
#   screen.rate= 0.015   weekly screening of infecteds (mean ~67 wk
#                        interval between screens for those in I)
#   tx.efficacy = 0.95   probability a treated index actually clears
#   ept.efficacy= 0.85   probability a partner takes EPT meds and clears
#   pn.test.prob= 0.85   sensitivity of returning-partner test in PR
#
# Pedagogical note: pn.lookback = 60 weeks is much longer than the real
# CDC guidance (60 days for chlamydia/gonorrhea). The partnership
# durations and act rates here are tuned for clean endemic-equilibrium
# teaching, not for calibrating to U.S. chlamydia data. Treat all numbers
# as illustrative.

init <- init.net(i.num = round(0.10 * n))

# control.net with the cumulative-edgelist machinery turned on:
#
#   cumulative.edgelist      = TRUE   -- engine collects edge history
#   truncate.el.cuml         = ...    -- destructive lookback truncation;
#                                        edges older than this many steps
#                                        are dropped from storage. We pass
#                                        the maximum lookback across the
#                                        scenarios so the edge history
#                                        covers every scenario's window.
#   save.cumulative.edgelist = TRUE   -- attach the edgelist to the
#                                        returned sim object for post-hoc
#                                        inspection.
max.lookback <- 120

make_control <- function(steps = nsteps) {
  control.net(
    type = NULL,
    nsims  = nsims,
    ncores = ncores,
    nsteps = steps,
    infection.FUN        = infect,
    screen.FUN           = screen,
    partner_services.FUN = partner_services,
    treat.FUN            = treat,
    recovery.FUN         = recov,
    cumulative.edgelist      = TRUE,
    truncate.el.cuml         = max.lookback,
    save.cumulative.edgelist = TRUE,
    verbose = FALSE
  )
}

control <- make_control()

# Helper for parameter sets. PN is off during burn-in (steps 1 .. pn.start)
# and on for the remainder.
make_param <- function(pn.arm = "none",
                       pn.trace.prob = 0,
                       pn.lookback = 60) {
  param.net(
    inf.prob      = 0.35,
    act.rate      = 1,
    rec.rate      = 0.01,
    screen.rate   = 0.015,
    tx.efficacy   = 0.95,
    ept.efficacy  = 0.85,
    pn.test.prob  = 0.85,
    pn.arm        = pn.arm,
    pn.trace.prob = pn.trace.prob,
    pn.lookback   = pn.lookback,
    pn.start      = pn.start
  )
}


# 3. Scenarios ---------------------------------------------------------------

# All five scenarios share the same network, the same disease parameters,
# and the same burn-in (steps 1..pn.start with PN off). They differ only
# in what happens after pn.start.

# Scenario 1: Baseline (screening only, no PN)
sim_none <- netsim(est,
                   make_param(pn.arm = "none", pn.trace.prob = 0),
                   init, control)
print(sim_none)

# Scenario 2: Patient Referral, 50% trace, 60-week lookback
sim_pr <- netsim(est,
                 make_param(pn.arm = "PR",
                            pn.trace.prob = 0.5,
                            pn.lookback = 60),
                 init, control)
print(sim_pr)

# Scenario 3: EPT, 50% trace, 60-week lookback
sim_ept <- netsim(est,
                  make_param(pn.arm = "EPT",
                             pn.trace.prob = 0.5,
                             pn.lookback = 60),
                  init, control)
print(sim_ept)

# Scenario 4: EPT, 50% trace, longer 120-week lookback
sim_ept_long <- netsim(est,
                       make_param(pn.arm = "EPT",
                                  pn.trace.prob = 0.5,
                                  pn.lookback = 120),
                       init, control)
print(sim_ept_long)

# Scenario 5: EPT, high trace (0.8) + longer lookback (120)
sim_ept_max <- netsim(est,
                      make_param(pn.arm = "EPT",
                                 pn.trace.prob = 0.8,
                                 pn.lookback = 120),
                      init, control)
print(sim_ept_max)


# 4. Analysis ----------------------------------------------------------------

sims <- list(
  none      = sim_none,
  pr        = sim_pr,
  ept       = sim_ept,
  ept_long  = sim_ept_long,
  ept_max   = sim_ept_max
)
labels <- c(none     = "Screening only",
            pr       = "PR, 50%, lookback 60",
            ept      = "EPT, 50%, lookback 60",
            ept_long = "EPT, 50%, lookback 120",
            ept_max  = "EPT, 80%, lookback 120")
cols <- c(none     = "gray40",
          pr       = "#8e44ad",
          ept      = "steelblue",
          ept_long = "#16a085",
          ept_max  = "firebrick")

sims <- lapply(sims, function(s) mutate_epi(s, prev = i.num / num))


## --- Plot 1: Prevalence over time, all scenarios overlaid ---
# Shared endemic burn-in (steps 1 .. pn.start) followed by post-PN
# bifurcation as each scenario settles to a new lower equilibrium.
par(mfrow = c(1, 1), mar = c(3, 3, 2, 1), mgp = c(2, 1, 0))
plot(sims$none, y = "prev",
     main = "Prevalence by Partner Notification Scenario",
     ylab = "Prevalence (I / N)", xlab = "Week",
     mean.col = cols["none"], mean.lwd = 2, mean.smooth = TRUE,
     qnts = FALSE, legend = FALSE,
     ylim = c(0, 1))
for (k in c("pr", "ept", "ept_long", "ept_max")) {
  plot(sims[[k]], y = "prev",
       mean.col = cols[k], mean.lwd = 2, mean.smooth = TRUE,
       qnts = FALSE, legend = FALSE, add = TRUE)
}
abline(v = pn.start, lty = 2, col = "gray40")
legend("topright", legend = labels, col = cols, lwd = 2, bty = "n",
       cex = 0.75)


## --- Cascade and reinfection helpers ---

# Per-scenario summary of the second-half period (post-PN). Uses the
# scenario data frames directly because the cascade fields are recorded
# per-step by the modules.
post_window <- function(df) df$time > pn.start

cascade_summary <- function(s) {
  df <- as.data.frame(s)
  df <- df[post_window(df), ]
  c(indices        = mean(tapply(df$n.indices,             df$sim, sum, na.rm = TRUE)),
    notified       = mean(tapply(df$n.partners.reached,    df$sim, sum, na.rm = TRUE)),
    treated        = mean(tapply(df$n.partner.tx,          df$sim, sum, na.rm = TRUE)),
    cleared_pn     = mean(tapply(df$n.partner.cleared.pn,  df$sim, sum, na.rm = TRUE)),
    wasted         = mean(tapply(df$n.partner.tx.wasted,   df$sim, sum, na.rm = TRUE)),
    nat_clear      = mean(tapply(df$is.flow.natural,       df$sim, sum, na.rm = TRUE)),
    index_clear    = mean(tapply(df$n.index.tx,            df$sim, sum, na.rm = TRUE)),
    eq_prev        = mean(df$prev[df$time > nsteps - 100], na.rm = TRUE))
}

casc <- sapply(sims, cascade_summary)


## --- Plot 2: Reinfection-per-index distribution ---
# For each ever-diagnosed individual, how many subsequent infections
# (reinfections) did they accumulate? We use the per-node `infections`
# counter, which is incremented every time the infection module flips a
# node S -> I. We then look at second-half infections only by reading
# si.flow per step, but the canonical "count of subsequent infections
# given ever diagnosed" requires the per-node attribute, which is only
# present in the saved attr.history. As a robust alternative that does
# not depend on attribute history, we use the running counts in the
# epi dataframe directly: among indices in the post-PN window, the
# expected reinfection rate is (subsequent si.flow per index). Plot a
# scenario boxplot of per-sim, per-window reinfection rate.
reinf_rate <- function(s) {
  df <- as.data.frame(s)
  df <- df[post_window(df), ]
  # Reinfections per ever-diagnosed person in this window (lower bound on
  # individual-level reinfection burden; if no one is diagnosed, NA).
  n_subsequent <- tapply(df$si.flow,       df$sim, sum, na.rm = TRUE)
  n_indices    <- tapply(df$n.indices,     df$sim, sum, na.rm = TRUE)
  ifelse(n_indices > 0, n_subsequent / n_indices, NA)
}
ri_by_scen <- lapply(sims, reinf_rate)

ri_mean <- sapply(ri_by_scen, mean, na.rm = TRUE)
xrange  <- range(unlist(ri_by_scen), na.rm = TRUE)

par(mfrow = c(1, 1), mar = c(4, 12, 2, 1), mgp = c(2.5, 1, 0))
plot(NA, xlim = c(0, xrange[2] * 1.05),
     ylim = c(0.5, length(ri_by_scen) + 0.5),
     yaxt = "n", xlab = "Subsequent infections per index (post-PN)",
     ylab = "",
     main = "Reinfection Pressure by Scenario")
axis(2, at = seq_along(ri_by_scen),
     labels = labels[names(ri_by_scen)],
     las = 1, cex.axis = 0.85)
abline(h = seq_along(ri_by_scen), col = "gray90", lty = 1)
for (i in seq_along(ri_by_scen)) {
  vals <- ri_by_scen[[i]]
  points(vals, rep(i, length(vals)),
         pch = 21, bg = adjustcolor(cols[names(ri_by_scen)[i]], 0.6),
         col = cols[names(ri_by_scen)[i]], cex = 1.2)
  segments(ri_mean[i], i - 0.25, ri_mean[i], i + 0.25,
           lwd = 3, col = cols[names(ri_by_scen)[i]])
}


## --- Plot 3: Tracing cascade bar chart ---
# Four bars per scenario: indices, partners notified, partners treated,
# partners cleared by PN. The point of the picture is the falloff at
# each step ("cascade") and how each scenario reshapes that funnel.
cascade_mat <- casc[c("indices", "notified", "treated", "cleared_pn"), ]
cascade_labels <- c(none     = "Screening\nonly",
                    pr       = "PR 50%\nlb 60",
                    ept      = "EPT 50%\nlb 60",
                    ept_long = "EPT 50%\nlb 120",
                    ept_max  = "EPT 80%\nlb 120")
colnames(cascade_mat) <- cascade_labels[colnames(cascade_mat)]

par(mfrow = c(1, 1), mar = c(4, 5, 4, 1), mgp = c(3.5, 1, 0))
bp <- barplot(cascade_mat, beside = TRUE,
              col = c("#7f8c8d", "#3498db", "#f39c12", "#27ae60"),
              names.arg = colnames(cascade_mat),
              las = 1, cex.names = 0.85,
              ylab = "Mean count over post-PN window",
              main = "Partner-Notification Cascade",
              ylim = c(0, max(cascade_mat, na.rm = TRUE) * 1.25))
legend("top",
       legend = c("Indices", "Notified", "Treated", "Cleared"),
       fill = c("#7f8c8d", "#3498db", "#f39c12", "#27ae60"),
       bty = "n", cex = 0.85, horiz = TRUE, inset = c(0, -0.08),
       xpd = TRUE)


## --- Summary table ---
# Equilibrium prevalence, per-index notified, per-notified treated,
# PN-attributable share of clearance, and reinfection rate.
totals <- as.data.frame(t(casc))
totals$Scenario <- labels
totals$Per_index_notified  <- ifelse(totals$indices > 0,
                                     totals$notified / totals$indices, NA)
totals$Per_notified_treated <- ifelse(totals$notified > 0,
                                      totals$treated / totals$notified, NA)
totals$PN_share_of_clearance <- with(totals,
                                     ifelse((cleared_pn + index_clear +
                                             nat_clear) > 0,
                                            cleared_pn /
                                            (cleared_pn + index_clear +
                                             nat_clear),
                                            NA))
totals$Reinf_per_index <- sapply(ri_by_scen, function(x) mean(x, na.rm = TRUE))

out <- data.frame(
  Scenario = totals$Scenario,
  EqPrev   = round(totals$eq_prev, 3),
  Indices  = round(totals$indices),
  Notified_per_idx = round(totals$Per_index_notified, 2),
  Treated_per_ntfd = round(totals$Per_notified_treated, 2),
  PN_share_clear   = round(totals$PN_share_of_clearance, 2),
  Reinf_per_idx    = round(totals$Reinf_per_index, 2),
  row.names = NULL
)
print(out)
