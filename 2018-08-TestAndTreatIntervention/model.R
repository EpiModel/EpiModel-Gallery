
##
## Test and Treat Intervention for an SIS Epidemic
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##
## Author: Samuel M. Jenness (Emory University)
## Date: August 2018
##

# Load EpiModel
suppressMessages(library(EpiModel))

# Standard Gallery unit test lines
rm(list = ls())
eval(parse(text = print(commandArgs(TRUE)[1])))

if (interactive()) {
  nsims <- 5
  ncores <- 5
  nsteps <- 500
} else {
  nsims <- 1
  ncores <- 1
  nsteps <- 50
}


# 1. Network Model Estimation ------------------------------------------------

# Initialize a 500-node network representing a sexual contact network
n <- 500
nw <- network_initialize(n)

# Formation model for an STI contact network:
#   edges: controls mean degree (175 edges → mean degree = 2*175/500 = 0.7)
#   concurrent: number of nodes with degree > 1 (110 = 22% of nodes have
#     overlapping partnerships, a key driver of STI transmission)
#   degrange(from = 4): constrains maximum degree to 3 (no node has 4+
#     partners simultaneously)
formation <- ~edges + concurrent + degrange(from = 4)
target.stats <- c(175, 110, 0)

# Partnership duration: mean of 50 weeks (~1 year)
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 50)
coef.diss

# Fit the ERGM
est <- netest(nw, formation, target.stats, coef.diss)

# Diagnostics
dx <- netdx(est, nsims = nsims, ncores = ncores, nsteps = nsteps,
            nwstats.formula = ~edges + concurrent + degree(0:5))
print(dx)
plot(dx)


# 2. Epidemic Model Parameters ------------------------------------------------

# Load custom module functions (testing/diagnosis, recovery)
source("2018-08-TestAndTreatIntervention/module-fx.R")

# Disease parameters — framed as a bacterial STI (e.g., gonorrhea):
#
#   Transmission (per-act and per-partnership probability):
#     inf.prob = 0.4: per-act transmission probability
#       (high, consistent with gonorrhea for receptive anal sex)
#     act.rate = 2: two sex acts per partnership per week
#     Per-partnership per-week probability: 1 - (1 - 0.4)^2 = 0.64
#
#   Recovery:
#     rec.rate = 1/20 = 0.05: natural clearance rate without treatment.
#       Mean infectious duration ~20 weeks. Many bacterial STIs can
#       self-resolve over months if untreated.
#     rec.rate.tx = 0.5: treatment recovery rate once diagnosed.
#       Mean duration on treatment ~2 weeks. Antibiotic treatment is
#       highly effective: 10x faster clearance than natural recovery.
#
#   Screening:
#     test.rate: proportion of undiagnosed who test each week.
#       This is the key intervention parameter, varied across scenarios.
#     test.dur = 2: diagnosis/treatment course lasts 2 weeks. After
#       this window, if infection has not cleared, the individual
#       returns to the untreated recovery rate until re-tested.

# Initial conditions: 10 infected nodes (2%) to seed the epidemic
init <- init.net(i.num = 10)

# Control settings: built-in infection module + custom tnt and recovery
control <- control.net(
  type = NULL,
  nsims = nsims,
  ncores = ncores,
  nsteps = nsteps,
  infection.FUN = infection.net,
  recovery.FUN = recov,
  tnt.FUN = tnt,
  verbose = FALSE
)


# 3. Scenario 1: No Screening (SIS Baseline) ---------------------------------

# Without any screening, only natural clearance (rec.rate = 0.05) drives
# recovery. The SIS model reaches an endemic equilibrium where new infections
# balance natural recoveries. This is the counterfactual — what happens
# without intervention.
param_none <- param.net(
  inf.prob = 0.4, act.rate = 2,
  rec.rate = 1 / 20, rec.rate.tx = 0.5,
  test.rate = 0, test.dur = 2
)

sim_none <- netsim(est, param_none, init, control)
print(sim_none)


# 4. Scenario 2: Standard Screening ------------------------------------------

# 10% of undiagnosed individuals test each week (~10-week average interval
# between tests). This is a moderate screening program — feasible for
# clinical settings with annual or semi-annual testing recommendations.
param_std <- param.net(
  inf.prob = 0.4, act.rate = 2,
  rec.rate = 1 / 20, rec.rate.tx = 0.5,
  test.rate = 0.1, test.dur = 2
)

sim_std <- netsim(est, param_std, init, control)
print(sim_std)


# 5. Scenario 3: Intensive Screening -----------------------------------------

# 30% of undiagnosed individuals test each week (~3-week average interval).
# This represents an aggressive public health intervention — frequent
# screening in high-prevalence populations. Shows how far screening
# can push down prevalence before diminishing returns.
param_int <- param.net(
  inf.prob = 0.4, act.rate = 2,
  rec.rate = 1 / 20, rec.rate.tx = 0.5,
  test.rate = 0.3, test.dur = 2
)

sim_int <- netsim(est, param_int, init, control)
print(sim_int)


# 6. Analysis -----------------------------------------------------------------

# Compute prevalence
sim_none <- mutate_epi(sim_none, prev = i.num / num)
sim_std <- mutate_epi(sim_std, prev = i.num / num)
sim_int <- mutate_epi(sim_int, prev = i.num / num)


## --- Plot 1: Prevalence Comparison (3 Scenarios) ---
# The central question: how does screening intensity affect population-level
# STI prevalence? No screening → endemic equilibrium at ~45%.
# Standard screening → reduced to ~30%. Intensive → near elimination.
par(mfrow = c(1, 1), mar = c(3, 3, 2, 1), mgp = c(2, 1, 0))
plot(sim_none, y = "prev",
     main = "Prevalence by Screening Intensity",
     ylab = "Prevalence (I / N)", xlab = "Week",
     mean.col = "firebrick", mean.lwd = 2, mean.smooth = TRUE,
     qnts.col = "firebrick", qnts.alpha = 0.2, qnts.smooth = TRUE,
     legend = FALSE)
plot(sim_std, y = "prev",
     mean.col = "steelblue", mean.lwd = 2, mean.smooth = TRUE,
     qnts.col = "steelblue", qnts.alpha = 0.2, qnts.smooth = TRUE,
     legend = FALSE, add = TRUE)
plot(sim_int, y = "prev",
     mean.col = "seagreen", mean.lwd = 2, mean.smooth = TRUE,
     qnts.col = "seagreen", qnts.alpha = 0.2, qnts.smooth = TRUE,
     legend = FALSE, add = TRUE)
legend("topright",
       legend = c("No Screening", "Standard (10%/wk)", "Intensive (30%/wk)"),
       col = c("firebrick", "steelblue", "seagreen"), lwd = 2, bty = "n")


## --- Plot 2: SIS Compartment Dynamics (No Screening vs. Intensive) ---
# Side-by-side comparison of S and I compartments. Without screening,
# the SIS model reaches a stable endemic equilibrium with ~225 infected.
# Intensive screening suppresses I to near zero.
par(mfrow = c(1, 2), mar = c(3, 3, 2, 1), mgp = c(2, 1, 0))
plot(sim_none, y = c("s.num", "i.num"),
     main = "No Screening",
     mean.col = c("steelblue", "firebrick"), mean.lwd = 2, mean.smooth = TRUE,
     qnts = FALSE, legend = TRUE)
plot(sim_int, y = c("s.num", "i.num"),
     main = "Intensive Screening",
     mean.col = c("steelblue", "firebrick"), mean.lwd = 2, mean.smooth = TRUE,
     qnts = FALSE, legend = TRUE)


## --- Plot 3: Diagnosis Coverage ---
# How many individuals are currently diagnosed at each timestep?
# This shows the test-and-treat cascade "in action" — the fraction of
# infected individuals who are in the treatment pipeline.
par(mfrow = c(1, 1), mar = c(3, 3, 2, 1), mgp = c(2, 1, 0))
plot(sim_std, y = "nDiag",
     main = "Currently Diagnosed (Treatment Pipeline)",
     ylab = "Number Diagnosed", xlab = "Week",
     mean.col = "steelblue", mean.lwd = 2, mean.smooth = TRUE,
     qnts.col = "steelblue", qnts.alpha = 0.2, qnts.smooth = TRUE,
     legend = FALSE)
plot(sim_int, y = "nDiag",
     mean.col = "#27ae60", mean.lwd = 2, mean.smooth = TRUE,
     qnts.col = "#27ae60", qnts.alpha = 0.2, qnts.smooth = TRUE,
     legend = FALSE, add = TRUE)
legend("topright",
       legend = c("Standard (10%/wk)", "Intensive (30%/wk)"),
       col = c("steelblue", "#27ae60"), lwd = 2, bty = "n")


## --- Plot 4: New Infections (Incidence) ---
# Incidence (new infections per week) shows the epidemic's momentum.
# Without screening, incidence stabilizes at a high endemic level.
# With screening, incidence drops as fewer infected individuals are
# available to transmit.
par(mfrow = c(1, 1), mar = c(3, 3, 2, 1), mgp = c(2, 1, 0))
plot(sim_none, y = "si.flow",
     main = "New Infections per Week",
     ylab = "New Infections", xlab = "Week",
     mean.col = "firebrick", mean.lwd = 2, mean.smooth = TRUE,
     qnts.col = "firebrick", qnts.alpha = 0.2, qnts.smooth = TRUE,
     legend = FALSE)
plot(sim_std, y = "si.flow",
     mean.col = "steelblue", mean.lwd = 2, mean.smooth = TRUE,
     qnts.col = "steelblue", qnts.alpha = 0.2, qnts.smooth = TRUE,
     legend = FALSE, add = TRUE)
plot(sim_int, y = "si.flow",
     mean.col = "#27ae60", mean.lwd = 2, mean.smooth = TRUE,
     qnts.col = "#27ae60", qnts.alpha = 0.2, qnts.smooth = TRUE,
     legend = FALSE, add = TRUE)
legend("topright",
       legend = c("No Screening", "Standard (10%/wk)", "Intensive (30%/wk)"),
       col = c("firebrick", "steelblue", "#27ae60"), lwd = 2, bty = "n")


## --- Summary Table ---
# Note: in an SIS model, cumulative new infections can paradoxically INCREASE
# with screening — treated individuals recover faster, become susceptible,
# and get reinfected. The right burden metric is person-weeks infected
# (total time spent in the I compartment across all individuals).
df_none <- as.data.frame(sim_none)
df_std <- as.data.frame(sim_std)
df_int <- as.data.frame(sim_int)

pw_none <- mean(tapply(df_none$i.num, df_none$sim, sum, na.rm = TRUE))
pw_std <- mean(tapply(df_std$i.num, df_std$sim, sum, na.rm = TRUE))
pw_int <- mean(tapply(df_int$i.num, df_int$sim, sum, na.rm = TRUE))

data.frame(
  Metric = c("Equilibrium prevalence",
             "Person-weeks infected",
             "Burden reduction (vs. none)",
             "Mean diagnosed (nDiag)"),
  No_Screening = c(
    round(mean(df_none$prev[df_none$time > nsteps * 0.8], na.rm = TRUE), 3),
    round(pw_none),
    "—",
    0
  ),
  Standard = c(
    round(mean(df_std$prev[df_std$time > nsteps * 0.8], na.rm = TRUE), 3),
    round(pw_std),
    paste0(round((1 - pw_std / pw_none) * 100), "%"),
    round(mean(df_std$nDiag[df_std$time > nsteps * 0.8], na.rm = TRUE), 1)
  ),
  Intensive = c(
    round(mean(df_int$prev[df_int$time > nsteps * 0.8], na.rm = TRUE), 3),
    round(pw_int),
    paste0(round((1 - pw_int / pw_none) * 100), "%"),
    round(mean(df_int$nDiag[df_int$time > nsteps * 0.8], na.rm = TRUE), 1)
  )
)


# 7. Sensitivity Analysis: Prevalence vs. Testing Rate ------------------------

# How does equilibrium prevalence respond to increasing screening intensity?
# This dose-response curve shows the marginal return of additional testing
# effort and helps identify the "sweet spot" for intervention design.
if (interactive()) {

  test_rates <- seq(0, 0.5, by = 0.05)
  eq_prev <- numeric(length(test_rates))

  for (i in seq_along(test_rates)) {
    cat("Running test.rate =", test_rates[i],
        "(", i, "/", length(test_rates), ")\n")
    param_i <- param.net(
      inf.prob = 0.4, act.rate = 2,
      rec.rate = 1 / 20, rec.rate.tx = 0.5,
      test.rate = test_rates[i], test.dur = 2
    )
    sim_i <- netsim(est, param_i, init, control)
    sim_i <- mutate_epi(sim_i, prev = i.num / num)
    df_i <- as.data.frame(sim_i)
    # Use last 20% of timesteps as equilibrium estimate
    eq_prev[i] <- mean(df_i$prev[df_i$time > nsteps * 0.8], na.rm = TRUE)
  }

  ## --- Plot 5: Dose-Response Curve ---
  # Indices for the three named scenarios (0, 0.1, 0.3)
  idx_none <- which.min(abs(test_rates - 0))
  idx_std  <- which.min(abs(test_rates - 0.1))
  idx_int  <- which.min(abs(test_rates - 0.3))

  par(mfrow = c(1, 1), mar = c(3, 3, 2, 1), mgp = c(2, 1, 0))
  plot(test_rates, eq_prev, type = "b", pch = 19, col = "steelblue",
       lwd = 2, xlab = "Weekly Testing Rate",
       ylab = "Equilibrium Prevalence",
       main = "Prevalence vs. Screening Intensity",
       ylim = c(0, max(eq_prev) * 1.15))
  abline(h = 0, lty = 2, col = "gray50")
  points(test_rates[c(idx_none, idx_std, idx_int)],
         eq_prev[c(idx_none, idx_std, idx_int)],
         pch = 19, col = "firebrick", cex = 1.5)
  text(0.05, eq_prev[idx_none] + 0.02, "No screening", col = "firebrick", cex = 0.8)
  text(0.15, eq_prev[idx_std] + 0.02, "Standard", col = "firebrick", cex = 0.8)
  text(0.35, eq_prev[idx_int] + 0.02, "Intensive", col = "firebrick", cex = 0.8)

  # Dose-response summary
  data.frame(
    test.rate = test_rates,
    eq.prevalence = round(eq_prev, 3)
  )
}
