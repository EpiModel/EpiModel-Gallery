
##
## SEIR/SEIRS Model: Adding an Exposed State to an SIR
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##
## Authors: Samuel M. Jenness, Venkata R. Duvvuri
## Date: August 2018
##

## Load EpiModel
suppressMessages(library(EpiModel))

# Standard Gallery unit test lines
rm(list = ls())
eval(parse(text = print(commandArgs(TRUE)[1])))

if (interactive()) {
  nsims <- 5
  ncores <- 5
  nsteps <- 800
} else {
  nsims <- 2
  ncores <- 2
  nsteps <- 200
}


# 1. Network Model Estimation ------------------------------------------------

# Initialize a 500-node network representing a contact network
n <- 500
nw <- network_initialize(n)

# Formation model: edges + isolates
# With 150 edges (mean degree = 2 * 150/500 = 0.6) and 240 isolates (48%
# of nodes with no active partnership), this produces a sparse network
# where roughly half the population has no current contact partner.
formation <- ~edges + isolates
target.stats <- c(150, 240)

# Partnership duration: mean of 25 timesteps. A longer duration reduces
# bias from the edges dissolution approximation used to fit the TERGM.
# This is a closed population (no vital dynamics), so no mortality
# adjustment is needed.
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 25)
coef.diss

# Fit the ERGM
est <- netest(nw, formation, target.stats, coef.diss)

# Diagnostics: verify the network simulation reproduces target statistics
dx <- netdx(est, nsims = nsims, ncores = ncores, nsteps = nsteps,
            nwstats.formula = ~edges + isolates + degree(0:5))
print(dx)
plot(dx)


# 2. Epidemic Model Parameters -----------------------------------------------

# Load custom module functions (infection, progression)
source("examples/seir-exposed-state/module-fx.R")

# Disease parameters (stylized for pedagogical clarity):
#
#   Transmission:
#     inf.prob = 0.5: per-act transmission probability (high, ensures
#       visible epidemic dynamics in a sparse network)
#     act.rate = 2: two acts per partnership per timestep
#     Per-timestep transmission probability per partnership:
#       finalProb = 1 - (1 - 0.5)^2 = 0.75
#
#   Progression rates (expressed as 1 / mean duration in compartment):
#     ei.rate = 1/50: mean latent (E) duration of 50 timesteps
#     ir.rate = 1/75: mean infectious (I) duration of 75 timesteps
#     rs.rate: 0 for SEIR (permanent immunity), 1/100 for SEIRS
#       (mean immune duration of 100 timesteps before returning to S)
#
#   These rates are chosen for pedagogical clarity: the progression
#   produces dynamics that are easy to visualize over the 800-step window.
#   For a specific pathogen, substitute published estimates (e.g.,
#   influenza: ei.rate ~ 1/2, ir.rate ~ 1/5).

# Initial conditions: 10 infected nodes (2% prevalence) to seed the epidemic
init <- init.net(i.num = 10)


# 3. Scenario 1: SEIR (Permanent Immunity) -----------------------------------

# In the SEIR model, recovered individuals gain permanent immunity.
# The epidemic burns through the susceptible population and dies out
# once the effective reproduction number drops below 1.
param_seir <- param.net(
  inf.prob = 0.5,
  act.rate = 2,
  ei.rate = 1 / 50,
  ir.rate = 1 / 75,
  rs.rate = 0
)

control <- control.net(
  type = NULL,
  nsteps = nsteps,
  nsims = nsims,
  ncores = ncores,
  infection.FUN = infect,
  progress.FUN = progress
)

sim_seir <- netsim(est, param_seir, init, control)
print(sim_seir)


# 4. Scenario 2: SEIRS (Waning Immunity) -------------------------------------

# Adding waning immunity (R -> S at rate 1/100) fundamentally changes the
# long-run dynamics. Instead of epidemic burnout, the disease persists
# indefinitely as recovered individuals return to the susceptible pool,
# creating sustained endemic transmission with periodic fluctuations.
param_seirs <- param.net(
  inf.prob = 0.5,
  act.rate = 2,
  ei.rate = 1 / 50,
  ir.rate = 1 / 75,
  rs.rate = 1 / 100
)

sim_seirs <- netsim(est, param_seirs, init, control)
print(sim_seirs)


# 5. Analysis -----------------------------------------------------------------

# Compute prevalence (proportion infectious)
sim_seir <- mutate_epi(sim_seir, prev = i.num / num)
sim_seirs <- mutate_epi(sim_seirs, prev = i.num / num)


## --- Plot 1: Prevalence Comparison ---
# The central question: how does waning immunity change disease burden?
# SEIR: prevalence rises then falls to zero as immunity accumulates.
# SEIRS: prevalence stabilizes at an endemic level as immunity wanes.
par(mfrow = c(1, 1), mar = c(3, 3, 2, 1), mgp = c(2, 1, 0))
plot(sim_seir, y = "prev",
     main = "Prevalence: SEIR vs. SEIRS",
     ylab = "Prevalence (I / N)", xlab = "Time Step",
     mean.col = "firebrick", mean.lwd = 2, mean.smooth = TRUE,
     qnts.col = "firebrick", qnts.alpha = 0.2, qnts.smooth = TRUE,
     legend = FALSE)
plot(sim_seirs, y = "prev",
     mean.col = "steelblue", mean.lwd = 2, mean.smooth = TRUE,
     qnts.col = "steelblue", qnts.alpha = 0.2, qnts.smooth = TRUE,
     legend = FALSE, add = TRUE)
legend("topright", legend = c("SEIR", "SEIRS"),
       col = c("firebrick", "steelblue"), lwd = 2, bty = "n")


## --- Plot 2: SEIR Compartment Dynamics ---
# All four compartments over time. Classic SEIR trajectory: S declines,
# E and I rise and fall, R accumulates toward the full population.
compartments <- c("s.num", "e.num", "i.num", "r.num")
titles <- c("Susceptible (S)", "Exposed (E)", "Infectious (I)", "Recovered (R)")
colors <- c("steelblue", "#f39c12", "firebrick", "#27ae60")

par(mfrow = c(2, 2), mar = c(3, 3, 2, 1), mgp = c(2, 1, 0))
for (j in seq_along(compartments)) {
  plot(sim_seir, y = compartments[j],
       main = titles[j], ylab = "Count", xlab = "Time Step",
       mean.col = colors[j], mean.lwd = 2, mean.smooth = TRUE,
       qnts.col = colors[j], qnts.alpha = 0.2, qnts.smooth = TRUE,
       legend = FALSE)
}


## --- Plot 3: SEIRS Compartment Dynamics ---
# With waning immunity, S does not deplete permanently and R does not
# accumulate forever. The system reaches a dynamic equilibrium with all
# four compartments maintaining nonzero populations.
par(mfrow = c(2, 2), mar = c(3, 3, 2, 1), mgp = c(2, 1, 0))
for (j in seq_along(compartments)) {
  plot(sim_seirs, y = compartments[j],
       main = titles[j], ylab = "Count", xlab = "Time Step",
       mean.col = colors[j], mean.lwd = 2, mean.smooth = TRUE,
       qnts.col = colors[j], qnts.alpha = 0.2, qnts.smooth = TRUE,
       legend = FALSE)
}


## --- Plot 4: New Infections (Incidence) ---
# SEIR: incidence peaks then declines to zero.
# SEIRS: incidence is sustained by the continuous resupply of susceptibles.
par(mfrow = c(1, 1), mar = c(3, 3, 2, 1), mgp = c(2, 1, 0))
plot(sim_seir, y = "se.flow",
     main = "New Infections: SEIR vs. SEIRS",
     ylab = "New Infections per Timestep", xlab = "Time Step",
     mean.col = "firebrick", mean.lwd = 2, mean.smooth = TRUE,
     qnts.col = "firebrick", qnts.alpha = 0.2, qnts.smooth = TRUE,
     legend = FALSE)
plot(sim_seirs, y = "se.flow",
     mean.col = "steelblue", mean.lwd = 2, mean.smooth = TRUE,
     qnts.col = "steelblue", qnts.alpha = 0.2, qnts.smooth = TRUE,
     legend = FALSE, add = TRUE)
legend("topright", legend = c("SEIR", "SEIRS"),
       col = c("firebrick", "steelblue"), lwd = 2, bty = "n")


## --- Summary Table ---
df_seir <- as.data.frame(sim_seir)
df_seirs <- as.data.frame(sim_seirs)

data.frame(
  Metric = c("Peak prevalence",
             "Mean prevalence",
             "Cumulative infections",
             "Final recovered count"),
  SEIR = c(
    round(max(tapply(df_seir$prev, df_seir$sim, max, na.rm = TRUE)), 3),
    round(mean(df_seir$prev, na.rm = TRUE), 3),
    round(mean(tapply(df_seir$se.flow, df_seir$sim, sum, na.rm = TRUE))),
    round(mean(df_seir$r.num[df_seir$time == max(df_seir$time)], na.rm = TRUE))
  ),
  SEIRS = c(
    round(max(tapply(df_seirs$prev, df_seirs$sim, max, na.rm = TRUE)), 3),
    round(mean(df_seirs$prev, na.rm = TRUE), 3),
    round(mean(tapply(df_seirs$se.flow, df_seirs$sim, sum, na.rm = TRUE))),
    round(mean(df_seirs$r.num[df_seirs$time == max(df_seirs$time)], na.rm = TRUE))
  )
)
