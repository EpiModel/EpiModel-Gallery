
##
## SEIRS Model with Leaky Vaccination and Vital Dynamics
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##
## Authors: Connor M. Van Meter (Emory University)
## Date: December 2018
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

# Initialize a 500-node network
n <- 500
nw <- network_initialize(n)

# Formation model: edges only (mean degree = 0.8)
formation <- ~edges
mean_degree <- 0.8
n_edges <- mean_degree * (n / 2)
target.stats <- c(n_edges)

# Departure rate: 0.008/week. Used in dissolution_coefs to adjust for
# population turnover -- without this correction, partnerships with
# departed nodes would artificially inflate the observed dissolution rate.
departure_rate <- 0.008

# Partnership duration: 50 weeks (~1 year)
coef.diss <- dissolution_coefs(dissolution = ~offset(edges),
                               duration = 50, d.rate = departure_rate)
coef.diss

# Fit the ERGM
est <- netest(nw, formation, target.stats, coef.diss)

# Diagnostics
dx <- netdx(est, nsims = nsims, ncores = ncores, nsteps = nsteps,
            nwstats.formula = ~edges + isolates + degree(0:5))
print(dx)
plot(dx)


# 2. Epidemic Model Parameters ------------------------------------------------

# Load custom module functions (infection with leaky protection, progression,
# departures, arrivals with vaccination)
source("examples/seirs-leaky-vaccination/module-fx.R")

# Disease and demographic parameters:
#
#   Transmission:
#     inf.prob = 0.5: high per-act probability
#     act.rate = 1: one contact per partnership per week
#     vaccine.efficacy = 0.8: vaccinated-protected individuals face
#       (1 - 0.8) * 0.5 = 0.1 per-act probability (80% reduction)
#
#   Disease progression (SEIRS):
#     ei.rate = 0.05: mean latent period ~20 weeks
#     ir.rate = 0.05: mean infectious period ~20 weeks
#     rs.rate = 0.05: mean immune period ~20 weeks (R -> S waning)
#     The SEIRS cycle allows reinfection -- recovered individuals
#     eventually lose immunity and return to the susceptible pool.
#     Vaccine-protected individuals who are reinfected retain their
#     protection through the cycle.
#
#   Vital dynamics:
#     departure.rate = 0.008: baseline weekly mortality
#     departure.disease.mult = 2: infected face 2x mortality
#     arrival.rate = 0.008: weekly birth rate (matched to departure rate
#       for approximately stable population size)
#
#   Vaccination (leaky):
#     Three routes: initialization, progression (campaigns), arrivals
#     Leaky vaccine REDUCES transmission probability rather than
#     preventing infection entirely (contrast with AON model).

# Initial conditions: 20 infected nodes (4% prevalence)
init <- init.net(i.num = 20)

# Control settings: custom modules with explicit execution order.
# Departures and arrivals run before infection so the network reflects
# the current population when transmission is simulated.
control <- control.net(
  type = NULL,
  nsteps = nsteps,
  nsims = nsims,
  ncores = ncores,
  infection.FUN = infect,
  progress.FUN = progress,
  departures.FUN = dfunc,
  arrivals.FUN = afunc,
  resimulate.network = TRUE,
  verbose = FALSE,
  module.order = c("resim_nets.FUN",
                   "departures.FUN", "arrivals.FUN",
                   "infection.FUN", "progress.FUN",
                   "prevalence.FUN")
)


# 3. Scenario 1: No Vaccination (SEIRS Baseline) -----------------------------

# No vaccination at any route. The SEIRS model produces endemic dynamics:
# waning immunity (R -> S at rs.rate) continuously replenishes the
# susceptible pool, sustaining ongoing transmission indefinitely.
param_novax <- param.net(
  inf.prob = 0.5, act.rate = 1,
  ei.rate = 0.05, ir.rate = 0.05, rs.rate = 0.05,
  departure.rate = departure_rate,
  departure.disease.mult = 2,
  arrival.rate = 0.008,
  vaccination.rate.initialization = 0,
  protection.rate.initialization = 0,
  vaccination.rate.progression = 0,
  protection.rate.progression = 0,
  vaccination.rate.arrivals = 0,
  protection.rate.arrivals = 0,
  vaccine.efficacy = 0
)

sim_novax <- netsim(est, param_novax, init, control)
print(sim_novax)


# 4. Scenario 2: Leaky Vaccination Program -----------------------------------

# Leaky vaccine with 80% efficacy: protected individuals face only 20%
# of the baseline transmission probability. Unlike AON vaccination,
# protected individuals CAN still be infected -- just at a reduced rate.
# If infected, they progress through E -> I -> R -> S normally, retaining
# their protection attribute for subsequent exposures.
param_vax <- param.net(
  inf.prob = 0.5, act.rate = 1,
  ei.rate = 0.05, ir.rate = 0.05, rs.rate = 0.05,
  departure.rate = departure_rate,
  departure.disease.mult = 2,
  arrival.rate = 0.008,
  vaccination.rate.initialization = 0.05,
  protection.rate.initialization = 0.8,
  vaccination.rate.progression = 0.05,
  protection.rate.progression = 0.8,
  vaccination.rate.arrivals = 0.6,
  protection.rate.arrivals = 0.8,
  vaccine.efficacy = 0.8
)

sim_vax <- netsim(est, param_vax, init, control)
print(sim_vax)


# 5. Analysis -----------------------------------------------------------------

# Compute derived epidemiological measures
sim_novax <- mutate_epi(sim_novax, prev = i.num / num)
sim_vax <- mutate_epi(sim_vax, prev = i.num / num)


## --- Plot 1: Prevalence Comparison ---
# The headline result: leaky vaccination reduces but does not eliminate
# disease. Unlike AON vaccination (which can achieve herd immunity),
# leaky vaccines allow breakthrough infections, so some transmission
# persists even at high coverage.
par(mfrow = c(1, 1), mar = c(3, 3, 2, 1), mgp = c(2, 1, 0))
plot(sim_novax, y = "prev",
     main = "Prevalence: No Vaccination vs. Leaky Vaccination",
     ylab = "Prevalence (I / N)", xlab = "Week",
     mean.col = "firebrick", mean.lwd = 2, mean.smooth = TRUE,
     qnts.col = "firebrick", qnts.alpha = 0.2, qnts.smooth = TRUE,
     legend = FALSE)
plot(sim_vax, y = "prev",
     mean.col = "#27ae60", mean.lwd = 2, mean.smooth = TRUE,
     qnts.col = "#27ae60", qnts.alpha = 0.2, qnts.smooth = TRUE,
     legend = FALSE, add = TRUE)
legend("topright",
       legend = c("No Vaccination", "Leaky Vaccination (80% eff.)"),
       col = c("firebrick", "#27ae60"), lwd = 2, bty = "n")


## --- Plot 2: SEIRS Compartments (No Vaccination) ---
# Without vaccination, SEIRS dynamics show endemic equilibrium: the
# R -> S waning immunity transition continuously replenishes susceptibles,
# sustaining ongoing transmission.
par(mfrow = c(2, 2), mar = c(3, 3, 2, 1), mgp = c(2, 1, 0))
plot(sim_novax, y = "s.num",
     main = "Susceptible (S)", ylab = "Count", xlab = "Week",
     mean.col = "steelblue", mean.lwd = 2, mean.smooth = TRUE,
     qnts.col = "steelblue", qnts.alpha = 0.2, qnts.smooth = TRUE,
     legend = FALSE)
plot(sim_novax, y = "e.num",
     main = "Exposed (E)", ylab = "Count", xlab = "Week",
     mean.col = "#f39c12", mean.lwd = 2, mean.smooth = TRUE,
     qnts.col = "#f39c12", qnts.alpha = 0.2, qnts.smooth = TRUE,
     legend = FALSE)
plot(sim_novax, y = "i.num",
     main = "Infectious (I)", ylab = "Count", xlab = "Week",
     mean.col = "firebrick", mean.lwd = 2, mean.smooth = TRUE,
     qnts.col = "firebrick", qnts.alpha = 0.2, qnts.smooth = TRUE,
     legend = FALSE)
plot(sim_novax, y = "r.num",
     main = "Recovered (R)", ylab = "Count", xlab = "Week",
     mean.col = "#27ae60", mean.lwd = 2, mean.smooth = TRUE,
     qnts.col = "#27ae60", qnts.alpha = 0.2, qnts.smooth = TRUE,
     legend = FALSE)


## --- Plot 3: SEIRS Compartments (With Leaky Vaccination) ---
# With vaccination, all compartments are reduced (fewer susceptibles
# become infected). Note: there is no "V" compartment as in AON --
# vaccine-protected individuals remain in S with reduced transmission.
par(mfrow = c(2, 2), mar = c(3, 3, 2, 1), mgp = c(2, 1, 0))
plot(sim_vax, y = "s.num",
     main = "Susceptible (S)", ylab = "Count", xlab = "Week",
     mean.col = "steelblue", mean.lwd = 2, mean.smooth = TRUE,
     qnts.col = "steelblue", qnts.alpha = 0.2, qnts.smooth = TRUE,
     legend = FALSE)
plot(sim_vax, y = "e.num",
     main = "Exposed (E)", ylab = "Count", xlab = "Week",
     mean.col = "#f39c12", mean.lwd = 2, mean.smooth = TRUE,
     qnts.col = "#f39c12", qnts.alpha = 0.2, qnts.smooth = TRUE,
     legend = FALSE)
plot(sim_vax, y = "i.num",
     main = "Infectious (I)", ylab = "Count", xlab = "Week",
     mean.col = "firebrick", mean.lwd = 2, mean.smooth = TRUE,
     qnts.col = "firebrick", qnts.alpha = 0.2, qnts.smooth = TRUE,
     legend = FALSE)
plot(sim_vax, y = "r.num",
     main = "Recovered (R)", ylab = "Count", xlab = "Week",
     mean.col = "#27ae60", mean.lwd = 2, mean.smooth = TRUE,
     qnts.col = "#27ae60", qnts.alpha = 0.2, qnts.smooth = TRUE,
     legend = FALSE)


## --- Plot 4: Vaccine Coverage Over Time ---
# Number of active nodes with vaccine protection. In the leaky model,
# these individuals are still counted in S (not moved to V).
par(mfrow = c(1, 1), mar = c(3, 3, 2, 1), mgp = c(2, 1, 0))
plot(sim_vax, y = "prt.num",
     main = "Vaccine-Protected Population Over Time",
     ylab = "Protected Individuals", xlab = "Week",
     mean.col = "#8e44ad", mean.lwd = 2, mean.smooth = TRUE,
     qnts.col = "#8e44ad", qnts.alpha = 0.2, qnts.smooth = TRUE,
     legend = FALSE)


## --- Plot 5: Incidence Comparison ---
# New infections (S -> E) per week. Leaky vaccination reduces but does
# not eliminate new infections -- breakthrough infections occur at a
# reduced rate.
par(mfrow = c(1, 1), mar = c(3, 3, 2, 1), mgp = c(2, 1, 0))
plot(sim_novax, y = "se.flow",
     main = "New Infections per Week",
     ylab = "New Infections (S -> E)", xlab = "Week",
     mean.col = "firebrick", mean.lwd = 2, mean.smooth = TRUE,
     qnts.col = "firebrick", qnts.alpha = 0.2, qnts.smooth = TRUE,
     legend = FALSE)
plot(sim_vax, y = "se.flow",
     mean.col = "#27ae60", mean.lwd = 2, mean.smooth = TRUE,
     qnts.col = "#27ae60", qnts.alpha = 0.2, qnts.smooth = TRUE,
     legend = FALSE, add = TRUE)
legend("topright",
       legend = c("No Vaccination", "Leaky Vaccination (80% eff.)"),
       col = c("firebrick", "#27ae60"), lwd = 2, bty = "n")


## --- Summary Table ---
df_novax <- as.data.frame(sim_novax)
df_vax <- as.data.frame(sim_vax)

late_novax <- df_novax$time > nsteps * 0.8
late_vax <- df_vax$time > nsteps * 0.8

data.frame(
  Metric = c("Mean prevalence",
             "Final prevalence",
             "Cumulative infections",
             "Mean population size",
             "Mean vaccine-protected"),
  No_Vaccination = c(
    round(mean(df_novax$prev[late_novax], na.rm = TRUE), 3),
    round(mean(df_novax$prev[df_novax$time == max(df_novax$time)],
               na.rm = TRUE), 3),
    round(mean(tapply(df_novax$se.flow, df_novax$sim, sum, na.rm = TRUE))),
    round(mean(df_novax$num[late_novax], na.rm = TRUE)),
    0
  ),
  Leaky_Vaccination = c(
    round(mean(df_vax$prev[late_vax], na.rm = TRUE), 3),
    round(mean(df_vax$prev[df_vax$time == max(df_vax$time)],
               na.rm = TRUE), 3),
    round(mean(tapply(df_vax$se.flow, df_vax$sim, sum, na.rm = TRUE))),
    round(mean(df_vax$num[late_vax], na.rm = TRUE)),
    round(mean(df_vax$prt.num[late_vax], na.rm = TRUE))
  )
)
