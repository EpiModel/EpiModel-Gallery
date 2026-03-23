
##
## HIV Transmission Model
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##
## Authors: Samuel M. Jenness, Connor Van Meter, Yuan Zhao, Emeli Anderson
## Date: March 2019
##

## Load EpiModel
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

# Initialize a network of 1000 nodes representing a sexually active population.
n <- 1000
nw <- network_initialize(n)

# Formation model: edges-only (simplest ERGM, producing a Bernoulli random
# graph). Mean degree of 0.8 means each person has on average 0.8 ongoing
# partnerships.
formation <- ~edges
mean_degree <- 0.8
target.stats <- c(mean_degree * (n / 2))

# Background departure rate: 0.003/person/week (~333 weeks or ~6.4 years
# average time in the population). This represents all-cause exit from the
# sexually active population (aging out, migration, background mortality).
departure_rate <- 0.003

# Average partnership duration of 52 weeks (1 year). The d.rate argument
# adjusts the dissolution coefficient so that the observed mean duration
# matches the target after accounting for competing risks from departures.
coef.diss <- dissolution_coefs(dissolution = ~offset(edges),
                               duration = 52, d.rate = departure_rate)
coef.diss

# Fit the ERGM
est <- netest(nw, formation, target.stats, coef.diss)

# Diagnostics: verify the simulated network reproduces target statistics
dx <- netdx(est, nsims = nsims, ncores = ncores, nsteps = nsteps)
print(dx)
plot(dx)


# 2. Epidemic Model Setup ----------------------------------------------------

# Load custom module functions (infection, progression, departure, arrival)
source("2019-03-HIV/module-fx.R")

# Transmission parameters (based on Granich et al. 2009):
#   - Per-act probability during chronic infection: 0.01
#   - Acute stage: 10x more infectious (high viral load)
#   - AIDS stage: 5x more infectious (rising viral load)
#   - ART: 95% reduction in infectiousness (multiplier = 0.05)
#   - 4 acts per partnership per week
#
# Progression rates (1 / mean weeks in stage):
#   - Acute:     ~12 weeks
#   - Chronic 1: ~260 weeks (~5 years)
#   - Chronic 2: ~260 weeks (~5 years)
#   - AIDS to departure: ~104 weeks (~2 years) without ART
#
# ART parameters:
#   - 1% per week of untreated infected start ART
#   - 0.5% per week of those on ART discontinue
#   - ART halves progression rate (multiplier = 0.5)

param_art <- param.net(
  inf.prob.chronic = 0.01,
  relative.inf.prob.acute = 10,
  relative.inf.prob.AIDS = 5,
  relative.inf.prob.ART = 0.05,
  act.rate = 4,
  AcuteToChronic1.Rate = 1 / 12,
  Chronic1ToChronic2.Rate = 1 / 260,
  Chronic2ToAIDS.Rate = 1 / 260,
  AIDSToDepart.Rate = 1 / 104,
  ART.Treatment.Rate = 0.01,
  ART.Discontinuance.Rate = 0.005,
  ART.Progression.Reduction.Rate = 0.5,
  arrival.rate = 0.002,
  departure.rate = departure_rate
)

# Initial conditions: 5% prevalence, all initially in the acute stage
init <- init.net(i.num = round(0.05 * n))

# Control settings:
#   type = NULL: fully custom module set (no built-in SIS/SIR/SI logic)
#   resimulate.network = TRUE: redraw the network from the fitted ERGM each
#     step, appropriate when vital dynamics change the active node set
#   module.order: progression runs before infection so that disease stage
#     and ART status are current before transmission probabilities are computed
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
  verbose = TRUE,
  module.order = c("resim_nets.FUN",
                   "progress.FUN",
                   "infection.FUN",
                   "departures.FUN",
                   "arrivals.FUN",
                   "prevalence.FUN")
)


# 3. Scenario 1: With ART ---------------------------------------------------

sim_art <- netsim(est, param_art, init, control)
print(sim_art)


# 4. Scenario 2: No ART (Baseline) ------------------------------------------

# Same parameters but ART treatment rate = 0. This counterfactual lets us
# evaluate the impact of ART on prevalence and incidence -- the central
# question of the Granich et al. model.
param_noart <- param.net(
  inf.prob.chronic = 0.01,
  relative.inf.prob.acute = 10,
  relative.inf.prob.AIDS = 5,
  relative.inf.prob.ART = 0.05,
  act.rate = 4,
  AcuteToChronic1.Rate = 1 / 12,
  Chronic1ToChronic2.Rate = 1 / 260,
  Chronic2ToAIDS.Rate = 1 / 260,
  AIDSToDepart.Rate = 1 / 104,
  ART.Treatment.Rate = 0,
  ART.Discontinuance.Rate = 0,
  ART.Progression.Reduction.Rate = 0.5,
  arrival.rate = 0.002,
  departure.rate = departure_rate
)

sim_noart <- netsim(est, param_noart, init, control)
print(sim_noart)


# 5. Analysis ----------------------------------------------------------------

# Compute derived epidemiological measures
sim_art <- mutate_epi(sim_art,
  prev = i.num / num,
  ir.rate = acute.flow / s.num,
  ART.num = acute.ART.num + chronic1.ART.num +
            chronic2.ART.num + AIDS.ART.num
)
sim_art <- mutate_epi(sim_art,
  ART.prev = ART.num / i.num
)

sim_noart <- mutate_epi(sim_noart,
  prev = i.num / num,
  ir.rate = acute.flow / s.num
)

# Total counts per disease stage (combining ART + no-ART sub-states)
sim_art <- mutate_epi(sim_art,
  acute.num = acute.ART.num + acute.NoART.num,
  chronic1.num = chronic1.ART.num + chronic1.NoART.num,
  chronic2.num = chronic2.ART.num + chronic2.NoART.num,
  AIDS.num = AIDS.ART.num + AIDS.NoART.num
)
sim_noart <- mutate_epi(sim_noart,
  acute.num = acute.ART.num + acute.NoART.num,
  chronic1.num = chronic1.ART.num + chronic1.NoART.num,
  chronic2.num = chronic2.ART.num + chronic2.NoART.num,
  AIDS.num = AIDS.ART.num + AIDS.NoART.num
)


## --- Plot 1: HIV Prevalence Comparison ---
# The central question: does ART reduce population-level prevalence?
par(mfrow = c(1, 1), mar = c(3, 3, 2, 1), mgp = c(2, 1, 0))
plot(sim_noart, y = "prev",
     main = "HIV Prevalence: ART vs. No ART",
     ylab = "Prevalence", xlab = "Weeks",
     mean.col = "firebrick", mean.lwd = 2, mean.smooth = TRUE,
     qnts.col = "firebrick", qnts.alpha = 0.2, qnts.smooth = TRUE,
     legend = FALSE)
plot(sim_art, y = "prev", add = TRUE,
     mean.col = "steelblue", mean.lwd = 2, mean.smooth = TRUE,
     qnts.col = "steelblue", qnts.alpha = 0.2, qnts.smooth = TRUE,
     legend = FALSE)
legend("topleft", legend = c("No ART", "With ART"),
       col = c("firebrick", "steelblue"), lwd = 2, bty = "n")


## --- Plot 2: HIV Incidence Rate Comparison ---
# Incidence rate = new infections per susceptible per week.
# ART reduces infectiousness, which should lower transmission.
plot(sim_noart, y = "ir.rate",
     main = "HIV Incidence Rate: ART vs. No ART",
     ylab = "Incidence Rate", xlab = "Weeks",
     mean.col = "firebrick", mean.lwd = 2, mean.smooth = TRUE,
     qnts.col = "firebrick", qnts.alpha = 0.2, qnts.smooth = TRUE,
     legend = FALSE)
plot(sim_art, y = "ir.rate", add = TRUE,
     mean.col = "steelblue", mean.lwd = 2, mean.smooth = TRUE,
     qnts.col = "steelblue", qnts.alpha = 0.2, qnts.smooth = TRUE,
     legend = FALSE)
legend("topleft", legend = c("No ART", "With ART"),
       col = c("firebrick", "steelblue"), lwd = 2, bty = "n")


## --- Plot 3: Disease Stage Counts by Scenario (2x2 Panel) ---
# One panel per disease stage, each comparing ART (blue) vs. No ART (red).
# With ART slowing progression, we expect fewer people reaching AIDS and
# more accumulating in earlier stages.
par(mfrow = c(2, 2), mar = c(3, 3, 2, 1), mgp = c(2, 1, 0))

stages <- c("acute.num", "chronic1.num", "chronic2.num", "AIDS.num")
titles <- c("Acute", "Chronic 1", "Chronic 2", "AIDS")

for (i in seq_along(stages)) {
  plot(sim_noart, y = stages[i],
       main = titles[i], ylab = "Count", xlab = "Weeks",
       mean.col = "firebrick", mean.lwd = 2, mean.smooth = TRUE,
       qnts.col = "firebrick", qnts.alpha = 0.2, qnts.smooth = TRUE,
       legend = FALSE)
  plot(sim_art, y = stages[i], add = TRUE,
       mean.col = "steelblue", mean.lwd = 2, mean.smooth = TRUE,
       qnts.col = "steelblue", qnts.alpha = 0.2, qnts.smooth = TRUE,
       legend = FALSE)
  legend("topleft", legend = c("No ART", "With ART"),
         col = c("firebrick", "steelblue"), lwd = 2, cex = 0.8, bty = "n")
}


## --- Plot 4: ART Coverage Among PLHIV ---
# What fraction of people living with HIV are on ART? This tracks toward
# the UNAIDS treatment coverage targets.
plot(sim_art, y = "ART.prev",
     main = "ART Coverage Among PLHIV",
     ylab = "Proportion on ART", xlab = "Weeks",
     mean.col = "steelblue", mean.lwd = 2, mean.smooth = TRUE,
     qnts.col = "steelblue", qnts.alpha = 0.2, qnts.smooth = TRUE,
     ylim = c(0, 1))


## --- Summary Table ---

# Use trajectory-level summaries (cumulative, mean) rather than endpoint
# values, since both scenarios converge to similar endemic equilibria.
df_art <- as.data.frame(sim_art)
df_noart <- as.data.frame(sim_noart)
last_t <- max(df_art$time)

cum_inf_art <- sum(df_art$acute.flow, na.rm = TRUE)
cum_inf_noart <- sum(df_noart$acute.flow, na.rm = TRUE)

# Summary statistics
data.frame(
  Metric = c("Cumulative new infections",
             "Infections averted by ART",
             "Mean prevalence",
             "Peak prevalence",
             "ART coverage at end (among PLHIV)"),
  With_ART = c(cum_inf_art,
               cum_inf_noart - cum_inf_art,
               round(mean(df_art$prev, na.rm = TRUE), 3),
               round(max(df_art$prev, na.rm = TRUE), 3),
               round(mean(df_art$ART.prev[df_art$time == last_t]), 3)),
  No_ART = c(cum_inf_noart,
             NA,
             round(mean(df_noart$prev, na.rm = TRUE), 3),
             round(max(df_noart$prev, na.rm = TRUE), 3),
             NA)
)
