
##
## Simple Cost-effectiveness Model (Simple SI model with cost and utility tracking)
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##
## Author: Gregory Knowlton (University of Minnesota)
## Date: October 2021
##

# Load EpiModel
suppressMessages(library(EpiModel))

# Set seed to ensure resulting CEA output has meaningful interpretation
# (i.e. the intervention is non-dominated)
set.seed(120792)

# Standard Gallery unit test lines
rm(list = ls())
eval(parse(text = print(commandArgs(TRUE)[1])))

if (interactive()) {
  nsims <- 5
  ncores <- 5
  nsteps <- 256
} else {
  nsims <- 2
  ncores <- 2
  nsteps <- 104
}

# Vital Dynamics Setup ----------------------------------------------------

# Yearly mortality rates by age: (1-84: 0.02, 85-99: 1.00, 100+: 2.00)
ages <- 0:100
dr <- 0.02
death_rate_pp_pw <- c(dr / 52,
                      50 * dr / 52,
                      100 * dr / 52)

# Build out a mortality rate vector
age_spans <- c(85, 15, 1)
dr_vec <- rep(death_rate_pp_pw, times = age_spans)
data.frame(ages, dr_vec)

# Network Model Estimation ------------------------------------------------

# Initialize the network
n <- 500
init.ages <- 16:85
ageVec <- sample(init.ages, n, replace = TRUE)

# Throughout simulation, individuals over age 65 are sexually inactive
active.sVec <- ifelse(ageVec <= 65, 1, 0)

# Initial population used to fit network must contain
# a mixture of sexually active and inactive individuals
# Here we assume age 65+ is inactive.
n.active.s <- length(which(ageVec <= 65))
n.inactive.s <- length(which(ageVec > 65))

# Initialize network
nw <- network_initialize(n)

# Set up age attribute
nw <- set_vertex_attribute(nw, "age", ageVec)

# active.s denotes which individuals are sexually active
nw <- set_vertex_attribute(nw, "active.s", active.sVec)

# Define the formation model: edges
# level 1 of the active.s corresponds to a value of 0,
# indicating non-participation in the sexual network
formation <- ~ edges + absdiff("age") + nodefactor("active.s", levels = 1)

# Input the appropriate target statistics for each term
mean_degree <- 0.8
edges <- mean_degree * (n.active.s / 2)
avg.abs.age.diff <- 1.5
absdiff <- edges * avg.abs.age.diff

# No edges should contain a sexually inactive ego
inactive.s.edges <- 0

target.stats <- c(edges, absdiff, inactive.s.edges)

# Parameterize the dissolution model
coef.diss <- dissolution_coefs(~ offset(edges), 60, mean(dr_vec))
coef.diss

# Fit the model
est <- netest(nw, formation, target.stats, coef.diss)

# Model diagnostics
dx <- netdx(est,
            nsims = nsims, ncores = ncores, nsteps = nsteps,
            nwstats.formula = ~ edges + absdiff("age") + isolates + degree(0:5) +
              nodefactor("active.s", levels = 1))
print(dx)
plot(dx)


# Epidemic model simulation -----------------------------------------------

# Epidemic model parameters for prophylaxis intervention scenario
end.horizon <- 52 + 14
param_inter <- param.net(inf.prob = 0.15,
                         death.rates = dr_vec,
                         cea.start = 14,
                         end.horizon = end.horizon,
                         arrival.rate = dr / 52,

                         # Intervention effectiveness/start time
                         inter.eff = 0.50,
                         inter.start = 14,
                         # Intervention Cost
                         inter.cost = 500000,

                         # Weekly costs by health status
                         sus.cost = 150,
                         inf.cost = 300,
                         # QALYs by health status
                         sus.qaly = 1.00,
                         inf.qaly = 0.75,

                         # QALY reduction per year of age
                         age.decrement = -0.003,
                         # Discount rate for costs and QALYs
                         disc.rate = 0.03)

# Initial conditions
init <- init.net(i.num = 50)

# Read in the module functions
if (interactive()) {
  source("2021-10-CostEffectivenessAnalysis/module-fx.R", echo = TRUE)
} else {
  source("module-fx.R")
}

# At the time step that marks the start of the end horizon, network
# resimulation will cease. Disease transmission is disabled at this point.
control.updater.list <- list(
  list(
    at = end.horizon,
    control = list(
      resimulate.network = FALSE
    )
  )
)

# Control settings
control <- control.net(type = NULL,
                       nsims = nsims,
                       ncores = ncores,
                       nsteps = nsteps,
                       aging.FUN = aging,
                       departures.FUN = dfunc,
                       arrivals.FUN = afunc,
                       infection.FUN = ifunc,
                       cea.FUN = costeffect,
                       resim_nets.FUN = resim_nets,
                       resimulate.network = TRUE,
                       control.updater.list = control.updater.list,
                       verbose = TRUE)

# Run the network model simulation with netsim
sim_inter <- netsim(est, param_inter, init, control)
print(sim_inter)

# Plot outcomes for prophylaxis intervention scenario
par(mfrow = c(1, 3))
plot(sim_inter, y = "d.flow", mean.smooth = TRUE, qnts = 1, main = "Departures")
plot(sim_inter, y = "a.flow", mean.smooth = TRUE, qnts = 1, main = "Arrivals")
plot(sim_inter, y = "si.flow", mean.smooth = TRUE, qnts = 1, main = "Infections")

par(mfrow = c(2, 2))
plot(sim_inter, y = "cost", mean.smooth = TRUE, qnts = 1, main = "Costs (undiscounted)")
plot(sim_inter, y = "qaly", mean.smooth = TRUE, qnts = 1, main = "QALYs (undiscounted)")
plot(sim_inter, y = "cost.disc", mean.smooth = TRUE, qnts = 1, main = "Costs (discounted)")
plot(sim_inter, y = "qaly.disc", mean.smooth = TRUE, qnts = 1, main = "QALYs (discounted)")


# Epidemic model parameters for no-intervention baseline (bl) scenario
param_bl <- param.net(inf.prob = 0.15,
                      death.rates = dr_vec,
                      end.horizon = end.horizon,
                      cea.start = 14,
                      arrival.rate = dr / 52,

                      # Weekly costs by health status
                      sus.cost = 150,
                      inf.cost = 300,
                      # QALYs by health status
                      sus.qaly = 1.00,
                      inf.qaly = 0.75,

                      # QALY reduction per year of age
                      age.decrement = -0.003,
                      # Discount rate for costs and QALYs
                      disc.rate = 0.03)

# Run the network model simulation with netsim for no intervention scenario
sim_bl <- netsim(est, param_bl, init, control)


# Cost-effectiveness Analysis -------------------------------------------------

# Load decision-analytic modeling package for cost-effectiveness analysis
suppressMessages(library(dampack))

# Cumulative discounted outcomes for each strategy are calculated
# The mean of these cumulative outcomes is the expected cost/effect
cost_bl <- as.data.frame(sim_bl, out = "mean")$cost.disc
qaly_bl <- as.data.frame(sim_bl, out = "mean")$qaly.disc
cost_inter <- as.data.frame(sim_inter, out = "mean")$cost.disc
qaly_inter <- as.data.frame(sim_inter, out = "mean")$qaly.disc
cost <- c(sum(cost_bl, na.rm = TRUE),
          sum(cost_inter, na.rm = TRUE))
effect <- c(sum(qaly_bl, na.rm = TRUE),
          sum(qaly_inter, na.rm = TRUE))
strategies <- c("No Intervention", "Universal Prophylaxis")


## Calculate ICER comparing competing strategies ##
icer_internal <- calculate_icers(cost = cost,
                                 effect = effect,
                                 strategies = strategies)
# The costs and QALYs used here were calculated at each time step using an
# explicit cost-effectiveness module within the EpiModel simulation.
icer_internal

# Plot ICER
plot(icer_internal)


## Calculating cumulative costs and effects external to EpiModel simulation ##
calc_outcomes <- function(sim, intervention) {

  # Define parameters
  end.horizon <- 52 + 14
  sus.cost <- 150
  inf.cost <- 300
  sus.qaly <- 1.00
  inf.qaly <- 0.75
  age.decrement <- -0.003
  disc.rate <- 0.03
  cea.start <- 14
  nsteps <- 104
  inter.cost <- 500000

  # Extract the number of susceptible/infected individuals across each simulation
  # and over time from the start of the analytic time horizon to the final time step.
  # Multiply these values by the per-person costs/QALYs accrued over one time step.
  pop.sus.cost <- unstack(as.data.frame(sim), s.num ~ sim)[(cea.start:nsteps) - 1, ] * sus.cost
  pop.inf.cost <- unstack(as.data.frame(sim), i.num ~ sim)[(cea.start:nsteps) - 1, ] * inf.cost
  pop.sus.qaly <- unstack(as.data.frame(sim), s.num ~ sim)[(cea.start:nsteps) - 1, ] * sus.qaly
  pop.inf.qaly <- unstack(as.data.frame(sim), i.num ~ sim)[(cea.start:nsteps) - 1, ] * inf.qaly

  # Mean age and population size across simulations and over time is used
  # to apply age-dependent reduction in QALYs
  meanAge <- unstack(as.data.frame(sim), meanAge ~ sim)[(cea.start:nsteps), ]
  pop.num <- unstack(as.data.frame(sim), num ~ sim)[(cea.start:nsteps) - 1, ]

  # Intervention cost is evenly distributed across analytic time horizon
  # The timing of these costs is important due to discounting
  if (intervention == TRUE) {
    inter.cost.vec <- c(rep(inter.cost / (end.horizon - cea.start),
                            end.horizon - cea.start),
                        rep(0,
                            nsteps - end.horizon + 1))
  } else {
    inter.cost.vec <- rep(0, nsteps - cea.start + 1)
  }

  # Combine age decrement and accrued QALYs from susceptible and infected
  pop.qaly <- ((meanAge * pop.num * age.decrement) +
                 pop.sus.qaly + pop.inf.qaly) / 52
  # Combine intervention costs and accrued costs from susceptible and infected
  pop.cost <- (pop.sus.cost + pop.inf.cost + inter.cost.vec)

  # Discount costs and QALYs
  pop.cost.disc <- pop.cost * (1 - disc.rate) ^ (0:(nsteps - cea.start) / 52)
  pop.qaly.disc <- pop.qaly * (1 - disc.rate) ^ (0:(nsteps - cea.start) / 52)

  # For each simulation, sum discounted costs and QALYs over time
  cuml.qaly.disc <- mean(colSums(pop.qaly.disc, na.rm = TRUE))
  cuml.cost.disc <- mean(colSums(pop.cost.disc, na.rm = TRUE))

  # cumulative costs and QALYs are returned
  return(list(cuml.qaly.disc = cuml.qaly.disc,
              cuml.cost.disc = cuml.cost.disc))
}

# Calculate discounted costs and QALYs for both strategies using function above
cost <- c(calc_outcomes(sim = sim_bl, intervention = FALSE)$cuml.cost.disc,
          calc_outcomes(sim = sim_inter, intervention = TRUE)$cuml.cost.disc)
effect <- c(calc_outcomes(sim = sim_bl, intervention = FALSE)$cuml.qaly.disc,
            calc_outcomes(sim = sim_inter, intervention = TRUE)$cuml.qaly.disc)
strategies <- c("No Intervention", "Universal Prophylaxis")


# Calculate incremental cost-effectiveness ratio comparing competing strategies
icer_external <- calculate_icers(cost = cost,
                                 effect = effect,
                                 strategies = strategies)
# Results are identical to those in icer_internal
icer_external
