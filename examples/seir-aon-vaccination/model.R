
##
## SEIR Model with All-or-Nothing Vaccination and Vital Dynamics
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##
## Authors: Connor M. Van Meter (Emory University)
## Date: October 2018
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

# Departure rate: 0.008/week (~2.4 year mean lifespan).
# This is high relative to real populations but creates visible turnover
# within the simulation window, demonstrating how vital dynamics interact
# with vaccination and disease transmission.
departure_rate <- 0.008

# Partnership duration: 50 weeks (~1 year).
# d.rate adjusts the dissolution coefficient for population turnover.
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

# Load custom module functions (vaccine init, infection, progression,
# departures, arrivals with vaccination)
source("examples/seir-aon-vaccination/module-fx.R")

# Disease and demographic parameters:
#
#   Transmission:
#     inf.prob = 0.5: high per-act probability (e.g., measles-like)
#     act.rate = 1: one contact per partnership per week
#
#   Disease progression:
#     ei.rate = 0.05: mean latent period ~20 weeks
#     ir.rate = 0.05: mean infectious period ~20 weeks
#     (No R -> S: SEIR model with permanent natural immunity)
#
#   Vital dynamics:
#     departure.rate = 0.008: baseline weekly mortality (~2.4 years)
#     departure.disease.mult = 2: infected individuals face 2x mortality
#     arrival.rate = 0.008: weekly birth rate (matched to departure rate
#       for approximately stable population size)
#
#   Vaccination (all-or-nothing):
#     Three routes: initialization, progression (campaigns), arrivals (newborns)
#     Each route has a vaccination rate and a protection rate.
#     vaccination.rate = probability of receiving the vaccine
#     protection.rate = probability the vaccine confers immunity
#     Effective protection rate = vaccination.rate * protection.rate

# Initial conditions: 20 infected nodes (4% prevalence)
init <- init.net(i.num = 20)

# Control settings: custom modules + resimulate.network for vital dynamics
control <- control.net(
  type = NULL,
  nsims = nsims,
  ncores = ncores,
  nsteps = nsteps,
  initAttr.FUN = init_vaccine_attrs,
  infection.FUN = infect,
  progress.FUN = progress,
  departures.FUN = dfunc,
  arrivals.FUN = afunc,
  resimulate.network = TRUE,
  verbose = FALSE,
  module.order = c("resim_nets.FUN", "initAttr.FUN",
                   "departures.FUN", "arrivals.FUN",
                   "infection.FUN", "progress.FUN",
                   "prevalence.FUN")
)


# 3. Scenario 1: No Vaccination (SEIR Baseline) ------------------------------

# No vaccination at any route. The SEIR epidemic spreads through the
# population limited only by natural immunity (R compartment) and
# demographic turnover. New susceptible arrivals fuel ongoing transmission.
param_novax <- param.net(
  inf.prob = 0.5, act.rate = 1,
  ei.rate = 0.05, ir.rate = 0.05,
  departure.rate = departure_rate,
  departure.disease.mult = 2,
  arrival.rate = 0.008,
  vaccination.rate.initialization = 0,
  protection.rate.initialization = 0,
  vaccination.rate.progression = 0,
  protection.rate.progression = 0,
  vaccination.rate.arrivals = 0,
  protection.rate.arrivals = 0
)

sim_novax <- netsim(est, param_novax, init, control)
print(sim_novax)


# 4. Scenario 2: Strong AON Vaccination Program ------------------------------

# Aggressive vaccination: 60% of newborns vaccinated (80% protection rate),
# 5% of unvaccinated population vaccinated per week via campaigns (80%
# protection rate). The effective weekly protection rate for newborns is
# 0.6 * 0.8 = 48%. Over time, the V (vaccine-immune) compartment grows,
# reducing the susceptible pool and creating herd immunity effects.
param_vax <- param.net(
  inf.prob = 0.5, act.rate = 1,
  ei.rate = 0.05, ir.rate = 0.05,
  departure.rate = departure_rate,
  departure.disease.mult = 2,
  arrival.rate = 0.008,
  vaccination.rate.initialization = 0.05,
  protection.rate.initialization = 0.8,
  vaccination.rate.progression = 0.05,
  protection.rate.progression = 0.8,
  vaccination.rate.arrivals = 0.6,
  protection.rate.arrivals = 0.8
)

sim_vax <- netsim(est, param_vax, init, control)
print(sim_vax)


# 5. Analysis -----------------------------------------------------------------

# Compute derived epidemiological measures
sim_novax <- mutate_epi(sim_novax, prev = i.num / num)
sim_vax <- mutate_epi(sim_vax, prev = i.num / num)


## --- Plot 1: Prevalence Comparison ---
# The headline result: vaccination dramatically reduces disease prevalence.
# Without vaccination, the epidemic saturates the population. With strong
# AON vaccination, herd immunity limits transmission even among the
# unvaccinated.
par(mfrow = c(1, 1), mar = c(3, 3, 2, 1), mgp = c(2, 1, 0))
plot(sim_novax, y = "prev",
     main = "Prevalence: No Vaccination vs. AON Vaccination",
     ylab = "Prevalence (I / N)", xlab = "Week",
     mean.col = "firebrick", mean.lwd = 2, mean.smooth = TRUE,
     qnts.col = "firebrick", qnts.alpha = 0.2, qnts.smooth = TRUE,
     legend = FALSE)
plot(sim_vax, y = "prev",
     mean.col = "seagreen", mean.lwd = 2, mean.smooth = TRUE,
     qnts.col = "seagreen", qnts.alpha = 0.2, qnts.smooth = TRUE,
     legend = FALSE, add = TRUE)
legend("topright",
       legend = c("No Vaccination", "Strong AON Vaccination"),
       col = c("firebrick", "seagreen"), lwd = 2, bty = "n")


## --- Plot 2: SEIR Compartments (No Vaccination) ---
# Without vaccination, the SEIR dynamics show the classic epidemic curve:
# susceptibles decline as exposed and infectious rise, then recovery
# transfers individuals to R. Vital dynamics continuously introduce
# new susceptibles via births, potentially fueling recurrent waves.
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
     mean.col = "seagreen", mean.lwd = 2, mean.smooth = TRUE,
     qnts.col = "seagreen", qnts.alpha = 0.2, qnts.smooth = TRUE,
     legend = FALSE)


## --- Plot 3: SEIR-V Compartments (With Vaccination) ---
# With vaccination, the V (vaccine-immune) compartment absorbs a growing
# fraction of the population. This reduces the effective susceptible pool,
# lowering the reproductive number and suppressing the epidemic.
par(mfrow = c(2, 3), mar = c(3, 3, 2, 1), mgp = c(2, 1, 0))
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
     mean.col = "seagreen", mean.lwd = 2, mean.smooth = TRUE,
     qnts.col = "seagreen", qnts.alpha = 0.2, qnts.smooth = TRUE,
     legend = FALSE)
plot(sim_vax, y = "v.num",
     main = "Vaccine-Immune (V)", ylab = "Count", xlab = "Week",
     mean.col = "#8e44ad", mean.lwd = 2, mean.smooth = TRUE,
     qnts.col = "#8e44ad", qnts.alpha = 0.2, qnts.smooth = TRUE,
     legend = FALSE)


## --- Plot 4: Vaccine Coverage Over Time ---
# The vaccine-immune compartment grows as vaccination campaigns and newborn
# vaccination accumulate protected individuals in the population.
par(mfrow = c(1, 1), mar = c(3, 3, 2, 1), mgp = c(2, 1, 0))
plot(sim_vax, y = "v.num",
     main = "Vaccine-Immune Population Over Time",
     ylab = "Count (V)", xlab = "Week",
     mean.col = "#8e44ad", mean.lwd = 2, mean.smooth = TRUE,
     qnts.col = "#8e44ad", qnts.alpha = 0.2, qnts.smooth = TRUE,
     legend = FALSE)


## --- Plot 5: Incidence Comparison ---
# New infections (S -> E) per week. Vaccination reduces the rate of new
# infections by shrinking the susceptible pool.
par(mfrow = c(1, 1), mar = c(3, 3, 2, 1), mgp = c(2, 1, 0))
plot(sim_novax, y = "se.flow",
     main = "New Infections per Week",
     ylab = "New Infections (S -> E)", xlab = "Week",
     mean.col = "firebrick", mean.lwd = 2, mean.smooth = TRUE,
     qnts.col = "firebrick", qnts.alpha = 0.2, qnts.smooth = TRUE,
     legend = FALSE)
plot(sim_vax, y = "se.flow",
     mean.col = "seagreen", mean.lwd = 2, mean.smooth = TRUE,
     qnts.col = "seagreen", qnts.alpha = 0.2, qnts.smooth = TRUE,
     legend = FALSE, add = TRUE)
legend("topright",
       legend = c("No Vaccination", "Strong AON Vaccination"),
       col = c("firebrick", "seagreen"), lwd = 2, bty = "n")


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
             "Mean vaccine-immune (V)"),
  No_Vaccination = c(
    round(mean(df_novax$prev[late_novax], na.rm = TRUE), 3),
    round(mean(df_novax$prev[df_novax$time == max(df_novax$time)],
               na.rm = TRUE), 3),
    round(mean(tapply(df_novax$se.flow, df_novax$sim, sum, na.rm = TRUE))),
    round(mean(df_novax$num[late_novax], na.rm = TRUE)),
    0
  ),
  Strong_AON = c(
    round(mean(df_vax$prev[late_vax], na.rm = TRUE), 3),
    round(mean(df_vax$prev[df_vax$time == max(df_vax$time)],
               na.rm = TRUE), 3),
    round(mean(tapply(df_vax$se.flow, df_vax$sim, sum, na.rm = TRUE))),
    round(mean(df_vax$num[late_vax], na.rm = TRUE)),
    round(mean(df_vax$v.num[late_vax], na.rm = TRUE))
  )
)
