
##
## Simple SI Model with Vital Dynamics
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
  nsteps <- 256
} else {
  nsims <- 1
  ncores <- 1
  nsteps <- 52
}


# 1. Vital Dynamics Setup -----------------------------------------------------

# Age range: 0 to 85 (one value per year of age)
ages <- 0:85

# Age-specific mortality rates per 100,000 per year for US population.
# Source: https://www.statista.com/statistics/241572/death-rate-by-age-and-sex-in-the-us/
# Age groups: <1, 1-4, 5-9, 10-14, 15-19, 20-24, 25-29, 30-34, 35-39,
#             40-44, 45-49, 50-54, 55-59, 60-64, 65-69, 70-74, 75-79,
#             80-84, 85+
departure_rate <- c(588.45, 24.8, 11.7, 14.55, 47.85, 88.2, 105.65, 127.2,
                    154.3, 206.5, 309.3, 495.1, 736.85, 1051.15, 1483.45,
                    2294.15, 3642.95, 6139.4, 13938.3)

# Convert from "per 100,000 per year" to "per person per week" (our timestep)
dep_rate_pw <- departure_rate / 1e5 / 52

# Expand the 19 age-group rates into a vector of 86 yearly rates.
# Each rate is repeated for the number of years in its age group:
#   1 year for <1, 4 years for 1-4, then 5-year groups through 80-84,
#   and 1 year for 85+.
age_spans <- c(1, 4, rep(5, 16), 1)
dr_vec <- rep(dep_rate_pw, times = age_spans)

# Visualize: mortality is low from ages 1-40, then rises exponentially
par(mar = c(3, 3, 2, 1), mgp = c(2, 1, 0), mfrow = c(1, 1))
barplot(dr_vec, col = "steelblue", border = NA,
        xlab = "Age", ylab = "Weekly Departure Rate",
        main = "Age-Specific Mortality Rates (US)")


# 2. Network Model Estimation ------------------------------------------------

# Initialize a 500-node network
n <- 500
nw <- network_initialize(n)

# Assign ages uniformly across 0-85 (a simplification; real populations
# have non-uniform age distributions)
ageVec <- sample(ages, n, replace = TRUE)
nw <- set_vertex_attribute(nw, "age", ageVec)

# Formation model: edges + absdiff("age")
# The absdiff term controls age-assortative mixing: the target statistic
# is the sum of |age_i - age_j| across all edges. A lower value means
# partnerships form preferentially between similarly-aged individuals.
formation <- ~edges + absdiff("age")

# Target statistics:
#   mean_degree = 0.8: average number of concurrent partnerships per node
#   avg.abs.age.diff = 1.5: mean absolute age difference within partnerships
#     (mild age assortativity -- partners tend to be close in age)
mean_degree <- 0.8
n_edges <- mean_degree * (n / 2)
avg.abs.age.diff <- 1.5
target.stats <- c(n_edges, n_edges * avg.abs.age.diff)

# Partnership duration: mean of 60 weeks (~1.2 years).
# The d.rate argument adjusts the dissolution coefficient for population
# turnover -- without this, partnerships with departed nodes would
# artificially inflate the observed dissolution rate.
coef.diss <- dissolution_coefs(~offset(edges), duration = 60,
                               d.rate = mean(dr_vec))
coef.diss

# Fit the ERGM
est <- netest(nw, formation, target.stats, coef.diss)

# Diagnostics: verify the network simulation reproduces target statistics
dx <- netdx(est, nsims = nsims, ncores = ncores, nsteps = nsteps,
            nwstats.formula = ~edges + absdiff("age") + isolates + degree(0:5))
print(dx)
plot(dx)


# 3. Epidemic Model Parameters ------------------------------------------------

# Load custom module functions (aging, departures, arrivals)
source("2018-08-SIwithVitalDynamics/module-fx.R")

# Epidemic and demographic parameters:
#
#   Transmission:
#     inf.prob = 0.15: per-act transmission probability
#     (uses built-in infection.net with default act.rate = 1)
#
#   Vital dynamics:
#     departure.rates = dr_vec: age-specific weekly mortality rates (86 values)
#     departure.disease.mult: multiplier applied to departure rates for
#       infected individuals. At mult=1, disease has no effect on survival.
#       At mult=50, infected individuals face 50x higher mortality --
#       exaggerated to show the population-level impact of disease-induced
#       mortality clearly.
#     arrival.rate = mean(dr_vec): weekly birth rate per person, set equal
#       to the mean departure rate so the population is approximately stable
#       in the absence of disease-induced excess mortality.

# Initial conditions: 50 infected nodes (10% prevalence) to seed the epidemic
init <- init.net(i.num = 50)


# 4. Scenario 1: No Disease-Induced Mortality (Baseline) ---------------------

# All-cause mortality only. Disease has no effect on survival, so the
# population size remains stable and prevalence rises toward saturation.
param_base <- param.net(
  inf.prob = 0.15,
  departure.rates = dr_vec,
  departure.disease.mult = 1,
  arrival.rate = mean(dr_vec)
)

# Control settings: custom vital dynamics modules + built-in infection module.
# resimulate.network = TRUE is required because ages change each timestep
# and the ERGM formation model includes absdiff("age").
control <- control.net(
  type = NULL,
  nsims = nsims,
  ncores = ncores,
  nsteps = nsteps,
  aging.FUN = aging,
  departures.FUN = dfunc,
  arrivals.FUN = afunc,
  infection.FUN = infection.net,
  resim_nets.FUN = resim_nets,
  resimulate.network = TRUE,
  verbose = FALSE
)

sim_base <- netsim(est, param_base, init, control)
print(sim_base)


# 5. Scenario 2: Disease-Induced Mortality (50x) -----------------------------

# Infected individuals face 50x the baseline age-specific mortality rate.
# This is exaggerated for pedagogical clarity: it shows how a highly lethal
# disease drives population decline when infected individuals die faster
# than births can replace them.
param_lethal <- param.net(
  inf.prob = 0.15,
  departure.rates = dr_vec,
  departure.disease.mult = 50,
  arrival.rate = mean(dr_vec)
)

sim_lethal <- netsim(est, param_lethal, init, control)
print(sim_lethal)


# 6. Analysis -----------------------------------------------------------------

# Compute prevalence (proportion infected)
sim_base <- mutate_epi(sim_base, prev = i.num / num)
sim_lethal <- mutate_epi(sim_lethal, prev = i.num / num)


## --- Plot 1: Prevalence Comparison ---
# Without disease mortality, SI prevalence saturates near 1 as everyone
# eventually becomes infected. With lethal disease, mortality removes
# infected individuals before they can transmit, slowing the epidemic.
par(mfrow = c(1, 1), mar = c(3, 3, 2, 1), mgp = c(2, 1, 0))
plot(sim_base, y = "prev",
     main = "Prevalence: Baseline vs. Disease Mortality",
     ylab = "Prevalence (I / N)", xlab = "Week",
     mean.col = "steelblue", mean.lwd = 2, mean.smooth = TRUE,
     qnts.col = "steelblue", qnts.alpha = 0.2, qnts.smooth = TRUE,
     legend = FALSE)
plot(sim_lethal, y = "prev",
     mean.col = "firebrick", mean.lwd = 2, mean.smooth = TRUE,
     qnts.col = "firebrick", qnts.alpha = 0.2, qnts.smooth = TRUE,
     legend = FALSE, add = TRUE)
legend("topleft", legend = c("Baseline (mult=1)", "Lethal (mult=50)"),
       col = c("steelblue", "firebrick"), lwd = 2, bty = "n")


## --- Plot 2: Population Size ---
# Baseline: population stable (arrivals ~ departures).
# Lethal: population declines because disease-induced mortality exceeds
# the birth rate calibrated to background mortality only.
par(mfrow = c(1, 1), mar = c(3, 3, 2, 1), mgp = c(2, 1, 0))
plot(sim_base, y = "num",
     main = "Population Size",
     ylab = "Total Population", xlab = "Week",
     mean.col = "steelblue", mean.lwd = 2, mean.smooth = TRUE,
     qnts.col = "steelblue", qnts.alpha = 0.2, qnts.smooth = TRUE,
     ylim = c(200, 600), legend = FALSE)
plot(sim_lethal, y = "num",
     mean.col = "firebrick", mean.lwd = 2, mean.smooth = TRUE,
     qnts.col = "firebrick", qnts.alpha = 0.2, qnts.smooth = TRUE,
     legend = FALSE, add = TRUE)
legend("topright", legend = c("Baseline (mult=1)", "Lethal (mult=50)"),
       col = c("steelblue", "firebrick"), lwd = 2, bty = "n")


## --- Plot 3: Vital Dynamics Detail (Disease Mortality Scenario) ---
# Four-panel view of the lethal disease scenario showing how disease-
# induced mortality affects the demographic processes.
par(mfrow = c(2, 2), mar = c(3, 3, 2, 1), mgp = c(2, 1, 0))

plot(sim_lethal, y = "meanAge",
     main = "Mean Age", ylab = "Years", xlab = "Week",
     mean.col = "firebrick", mean.lwd = 2, mean.smooth = TRUE,
     qnts.col = "firebrick", qnts.alpha = 0.2, qnts.smooth = TRUE,
     legend = FALSE)

plot(sim_lethal, y = "d.flow",
     main = "Departures per Week", ylab = "Count", xlab = "Week",
     mean.col = "firebrick", mean.lwd = 2, mean.smooth = TRUE,
     qnts.col = "firebrick", qnts.alpha = 0.2, qnts.smooth = TRUE,
     legend = FALSE)

plot(sim_lethal, y = "a.flow",
     main = "Arrivals per Week", ylab = "Count", xlab = "Week",
     mean.col = "#27ae60", mean.lwd = 2, mean.smooth = TRUE,
     qnts.col = "#27ae60", qnts.alpha = 0.2, qnts.smooth = TRUE,
     legend = FALSE)

plot(sim_lethal, y = "s.num",
     main = "Susceptible Count", ylab = "Count", xlab = "Week",
     mean.col = "steelblue", mean.lwd = 2, mean.smooth = TRUE,
     qnts.col = "steelblue", qnts.alpha = 0.2, qnts.smooth = TRUE,
     legend = FALSE)


## --- Plot 4: New Infections ---
par(mfrow = c(1, 1), mar = c(3, 3, 2, 1), mgp = c(2, 1, 0))
plot(sim_base, y = "si.flow",
     main = "New Infections: Baseline vs. Disease Mortality",
     ylab = "New Infections per Week", xlab = "Week",
     mean.col = "steelblue", mean.lwd = 2, mean.smooth = TRUE,
     qnts.col = "steelblue", qnts.alpha = 0.2, qnts.smooth = TRUE,
     legend = FALSE)
plot(sim_lethal, y = "si.flow",
     mean.col = "firebrick", mean.lwd = 2, mean.smooth = TRUE,
     qnts.col = "firebrick", qnts.alpha = 0.2, qnts.smooth = TRUE,
     legend = FALSE, add = TRUE)
legend("topright", legend = c("Baseline (mult=1)", "Lethal (mult=50)"),
       col = c("steelblue", "firebrick"), lwd = 2, bty = "n")


## --- Summary Table ---
df_base <- as.data.frame(sim_base)
df_lethal <- as.data.frame(sim_lethal)

data.frame(
  Metric = c("Mean prevalence",
             "Final prevalence",
             "Mean population size",
             "Final population size",
             "Cumulative departures",
             "Cumulative infections"),
  Baseline = c(
    round(mean(df_base$prev, na.rm = TRUE), 3),
    round(mean(df_base$prev[df_base$time == max(df_base$time)], na.rm = TRUE), 3),
    round(mean(df_base$num, na.rm = TRUE)),
    round(mean(df_base$num[df_base$time == max(df_base$time)], na.rm = TRUE)),
    round(mean(tapply(df_base$d.flow, df_base$sim, sum, na.rm = TRUE))),
    round(mean(tapply(df_base$si.flow, df_base$sim, sum, na.rm = TRUE)))
  ),
  Lethal = c(
    round(mean(df_lethal$prev, na.rm = TRUE), 3),
    round(mean(df_lethal$prev[df_lethal$time == max(df_lethal$time)], na.rm = TRUE), 3),
    round(mean(df_lethal$num, na.rm = TRUE)),
    round(mean(df_lethal$num[df_lethal$time == max(df_lethal$time)], na.rm = TRUE)),
    round(mean(tapply(df_lethal$d.flow, df_lethal$sim, sum, na.rm = TRUE))),
    round(mean(tapply(df_lethal$si.flow, df_lethal$sim, sum, na.rm = TRUE)))
  )
)
