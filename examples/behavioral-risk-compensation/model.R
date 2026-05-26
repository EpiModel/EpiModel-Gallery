
##
## SIR with Behavioral Risk Compensation During Illness
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##
## Author: Samuel M. Jenness (Emory University)
## Date: May 2026
##

suppressMessages(library(EpiModel))

# Standard Gallery unit test lines
rm(list = ls())
eval(parse(text = print(commandArgs(TRUE)[1])))

# Two run modes:
#   interactive(): full network, 5 sims, longer horizon, full bisection
#   non-interactive (CI): small network, 1 sim, short horizon, light bisection
if (interactive()) {
  nsims    <- 5
  ncores   <- 5
  nsteps   <- 300
  n        <- 500
  cal_iter <- 6
} else {
  nsims    <- 1
  ncores   <- 1
  nsteps   <- 120
  n        <- 250
  cal_iter <- 3
}


# 1. Network Model Estimation ----------------------------------------------

# A simple contact network. Mean degree 2, partnership duration 30
# timesteps (treat each step as a day). Network structure is intentionally
# uniform so that the pedagogical focus stays on the behavior-during-
# infection mechanism rather than on network heterogeneity.
nw <- network_initialize(n)
formation <- ~edges
target.stats <- c(round(2 * n / 2))
coef.diss <- dissolution_coefs(~offset(edges), duration = 30)
est <- netest(nw, formation, target.stats, coef.diss)

dx <- netdx(est, nsims = nsims, ncores = ncores, nsteps = nsteps,
            nwstats.formula = ~edges + degree(0:3))
print(dx)
plot(dx)


# 2. Parameters and Modules -----------------------------------------------

source("examples/behavioral-risk-compensation/module-fx.R")

# Disease parameterization (each timestep = 1 day):
#   el.rate = 1/3   ->  mean early phase  ~3 days
#   lr.rate = 1/5   ->  mean late  phase  ~5 days
#   Total infectious period ~8 days.
#
# Behavior multipliers applied to act.rate when the infected partner is in
# each sub-stage:
#   mult.early = 0.3   (heavy contact reduction at peak symptoms)
#   mult.late  = 0.6   (partial recovery toward baseline)
# Setting both = 1.0 recovers a behavior-naive SIR.

make_param <- function(inf.prob,
                       mult.early = 1.0,
                       mult.late  = 1.0) {
  param.net(
    inf.prob   = inf.prob,
    act.rate   = 1,
    el.rate    = 1 / 3,
    lr.rate    = 1 / 5,
    mult.early = mult.early,
    mult.late  = mult.late
  )
}

init <- init.net(i.num = 10)

control <- control.net(
  type = NULL,
  nsims = nsims,
  ncores = ncores,
  nsteps = nsteps,
  infection.FUN = infect,
  progress.FUN  = progress,
  verbose = FALSE
)


# 3. Calibration -----------------------------------------------------------

# Both models calibrated to the same target cumulative attack rate so the
# downstream NPI comparison is like-for-like.

# Cumulative attack rate = total new infections / population, averaged
# across simulations.
cum_attack <- function(sim) {
  df  <- as.data.frame(sim)
  inc <- tapply(df$si.flow, df$sim, sum, na.rm = TRUE)
  pop <- tapply(df$num,     df$sim, function(x) x[length(x)])
  mean(inc / pop, na.rm = TRUE)
}

# Bisection on inf.prob (per-act transmission probability). Closed
# population, monotone increasing relationship between inf.prob and final
# attack rate, so bisection converges reliably.
calibrate_beta <- function(target, mult.early, mult.late,
                           lo = 0.01, hi = 0.60,
                           max_iter = cal_iter) {
  mid <- (lo + hi) / 2
  for (i in seq_len(max_iter)) {
    mid <- (lo + hi) / 2
    p   <- make_param(inf.prob = mid,
                      mult.early = mult.early,
                      mult.late  = mult.late)
    sim <- netsim(est, p, init, control)
    ar  <- cum_attack(sim)
    cat(sprintf("  iter %d  inf.prob=%.4f  attack=%.3f\n",
                i, mid, ar))
    if (is.na(ar) || ar < target) lo <- mid else hi <- mid
  }
  (lo + hi) / 2
}

target_attack <- 0.40

cat("\nCalibrating DYNAMIC model (mult.early=0.3, mult.late=0.6) ",
    "to attack=", target_attack, "\n", sep = "")
beta_dyn <- calibrate_beta(target_attack,
                           mult.early = 0.3, mult.late = 0.6)

cat("\nCalibrating NAIVE model (mult.early=1.0, mult.late=1.0) ",
    "to attack=", target_attack, "\n", sep = "")
beta_naive <- calibrate_beta(target_attack,
                             mult.early = 1.0, mult.late = 1.0)

cat("\nCalibrated per-act transmission probability:\n",
    "  dynamic = ", round(beta_dyn,   4), "\n",
    "  naive   = ", round(beta_naive, 4), "\n",
    "  ratio   = ", round(beta_dyn / beta_naive, 3),
    "  (dynamic / naive)\n", sep = "")


# 4. Calibrated Baseline Scenarios ----------------------------------------

sim_dyn   <- netsim(est,
                    make_param(beta_dyn,
                               mult.early = 0.3, mult.late = 0.6),
                    init, control)
sim_naive <- netsim(est,
                    make_param(beta_naive,
                               mult.early = 1.0, mult.late = 1.0),
                    init, control)

print(sim_dyn)
print(sim_naive)


# 5. NPI Scenario: Isolation of Symptomatic Cases -------------------------

# Intervention drives the early-stage multiplier down to a fixed target
# level iso.mult (household-only contacts during the most symptomatic
# phase). The naive model's early-stage baseline is 1.0; the dynamic
# model's is already 0.3. Same intervention, different starting point,
# different absolute reduction.

iso.mult <- 0.1

sim_dyn_npi   <- netsim(est,
                        make_param(beta_dyn,
                                   mult.early = iso.mult,
                                   mult.late  = 0.6),
                        init, control)
sim_naive_npi <- netsim(est,
                        make_param(beta_naive,
                                   mult.early = iso.mult,
                                   mult.late  = 1.0),
                        init, control)


# 6. Analysis --------------------------------------------------------------

ar_dyn       <- cum_attack(sim_dyn)
ar_dyn_npi   <- cum_attack(sim_dyn_npi)
ar_naive     <- cum_attack(sim_naive)
ar_naive_npi <- cum_attack(sim_naive_npi)

# Percent reduction in cumulative incidence attributable to the NPI.
red_dyn   <- 100 * (ar_dyn   - ar_dyn_npi)   / ar_dyn
red_naive <- 100 * (ar_naive - ar_naive_npi) / ar_naive

summary_tbl <- data.frame(
  Model       = c("Naive", "Dynamic"),
  inf.prob    = round(c(beta_naive,    beta_dyn),    4),
  Attack_base = round(c(ar_naive,      ar_dyn),      3),
  Attack_NPI  = round(c(ar_naive_npi,  ar_dyn_npi),  3),
  Pct_reduced = round(c(red_naive,     red_dyn),     1),
  row.names = NULL
)
print(summary_tbl)


## --- Plot 1: Prevalence over time (calibrated baselines) ---
# Same calibration target, very different per-act transmissibility.
par(mfrow = c(1, 1), mar = c(3, 3, 2, 1), mgp = c(2, 1, 0))
plot(sim_naive, y = "i.num",
     main = "Calibrated Baseline Prevalence",
     ylab = "Number infectious",
     xlab = "Time step (days)",
     mean.col = "firebrick", mean.lwd = 2, mean.smooth = TRUE,
     qnts = FALSE, legend = FALSE)
plot(sim_dyn, y = "i.num",
     mean.col = "steelblue", mean.lwd = 2, mean.smooth = TRUE,
     qnts = FALSE, legend = FALSE, add = TRUE)
legend("topright",
       legend = c(sprintf("Naive (beta=%.3f)",   beta_naive),
                  sprintf("Dynamic (beta=%.3f)", beta_dyn)),
       col = c("firebrick", "steelblue"), lwd = 2, bty = "n")


## --- Plot 2: Calibrated transmissibility and NPI impact ---
par(mfrow = c(1, 2), mar = c(5, 4, 3, 1), mgp = c(2.4, 1, 0))

# Left: calibrated per-act transmission probability
bp1 <- barplot(c(Naive = beta_naive, Dynamic = beta_dyn),
               col = c("firebrick", "steelblue"),
               ylab = "Calibrated per-act transmission probability",
               main = "Calibrated Transmissibility",
               ylim = c(0, max(beta_naive, beta_dyn) * 1.25))
text(bp1, c(beta_naive, beta_dyn) * 1.07,
     sprintf("%.3f", c(beta_naive, beta_dyn)), font = 2)

# Right: baseline vs NPI attack rate, grouped by model. Plotting the
# actual attack rates (always non-negative) avoids the conceptual oddity
# of a "negative percent reduction" when stochastic noise dominates a
# small effect, and visually preserves the calibration target.
attack_mat <- matrix(c(ar_naive, ar_naive_npi,
                       ar_dyn,   ar_dyn_npi),
                     nrow = 2,
                     dimnames = list(c("Baseline", "NPI"),
                                     c("Naive", "Dynamic")))
cols_pair <- c(adjustcolor("firebrick", alpha.f = 0.45), "firebrick",
               adjustcolor("steelblue", alpha.f = 0.45), "steelblue")
ylim2 <- c(0, max(attack_mat, target_attack, na.rm = TRUE) * 1.3)

bp2 <- barplot(attack_mat, beside = TRUE,
               col = cols_pair, las = 1,
               ylab = "Cumulative attack rate",
               main = "Attack Rate: Baseline vs NPI",
               ylim = ylim2)

text(bp2, attack_mat + diff(ylim2) * 0.025,
     sprintf("%.2f", attack_mat), font = 2, cex = 0.9)

abline(h = target_attack, lty = 2, col = "gray40")

pct_change <- -c(red_naive, red_dyn)
mtext(side = 1, line = 2.4, at = colMeans(bp2),
      text = sprintf("%+.1f%%", pct_change),
      font = 2, cex = 0.95,
      col = ifelse(pct_change <= 0, "darkgreen", "firebrick"))

legend("topright",
       legend = c("Baseline", "+ NPI"),
       fill = c(adjustcolor("gray30", alpha.f = 0.45), "gray30"),
       border = NA, bty = "n", cex = 0.85)
