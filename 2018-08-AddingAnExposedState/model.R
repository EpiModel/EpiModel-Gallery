
##
## SEIR Model: Adding an Exposed State to an SIR
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


# Network model estimation ------------------------------------------------

# Initialize the network
nw <- network_initialize(500)

# Define the formation model: edges + isolates (number with degree of 0)
formation = ~edges + isolates

# Input the appropriate target statistics for each term
target.stats <- c(150, 240)

# Parameterize the dissolution model
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 10)
coef.diss

# Fit the model
est <- netest(nw, formation, target.stats, coef.diss)

# Model diagnostics
dx <- netdx(est, nsims = 2, ncores = 2, nsteps = 500,
            nwstats.formula = ~edges + isolates + degree(0:5))
print(dx)
plot(dx, plots.joined = FALSE, qnts.alpha = 0.8)


# Epidemic model simulation -----------------------------------------------

# Model parameters
param <- param.net(inf.prob = 0.5, act.rate = 2,
                   ei.rate = 0.01, ir.rate = 0.01)

# Initial conditions
init <- init.net(i.num = 10)

# Read in the module functions
if (interactive()) {
  source("2018-08-AddingAnExposedState/module-fx.R", echo = TRUE)
} else {
  source("module-fx.R")
}

# Control settings
control <- control.net(type = NULL,
                       nsteps = nsteps,
                       nsims = nsims,
                       ncores = ncores,
                       infection.FUN = infect,
                       progress.FUN = progress,
                       recovery.FUN = recovery.net,
                       departures.FUN = departures.net,
                       arrivals.FUN = arrivals.net,
                       prevalence.FUN = prevalence.net)

# Run the network model simulation with netsim
sim <- netsim(est, param, init, control)
print(sim)

# Plot outcomes
par(mar = c(3,3,1,1), mgp = c(2,1,0))
plot(sim,
     mean.col = 1:4, mean.lwd = 1, mean.smooth = FALSE,
     qnts = 1, qnts.col = 1:4, qnts.alpha = 0.25, qnts.smooth = FALSE,
     legend = TRUE)

plot(sim, y = c("se.flow", "ei.flow", "ir.flow"),
     mean.col = 1:4, mean.lwd = 1, mean.smooth = TRUE,
     qnts.col = 1:4, qnts.alpha = 0.25, qnts.smooth = TRUE,
     ylim = c(0, 3), legend = TRUE)

# Average across simulations at beginning, middle, end
df <- as.data.frame(sim)
df[c(2, 100, 500), ]



# Extension #1: Adding an R --> S Transition (SEIRS) ----------------------

# Model parameters
param <- param.net(inf.prob = 0.5, act.rate = 2,
                   ei.rate = 0.01, ir.rate = 0.01,
                   rs.rate = 0.005)

# Initial conditions
init <- init.net(i.num = 10)

# Read in the module functions
source("module-fx.R", echo = TRUE)

# Control settings
control <- control.net(nsteps = nsteps,
                       nsims = nsims,
                       ncores = ncores,
                       infection.FUN = infect,
                       progress.FUN = progress2,
                       recovery.FUN = recovery.net,
                       departures.FUN = departures.net,
                       arrivals.FUN = arrivals.net,
                       prevalence.FUN = prevalence.net)

# Run the network model simulation with netsim
sim <- netsim(est, param, init, control)
print(sim)

# Plot outcomes
par(mar = c(3,3,1,1), mgp = c(2,1,0))
plot(sim,
     mean.col = 1:4, mean.lwd = 1, mean.smooth = FALSE,
     qnts = 1, qnts.col = 1:4, qnts.alpha = 0.25, qnts.smooth = FALSE,
     legend = TRUE)

plot(sim, y = c("se.flow", "ei.flow", "ir.flow", "rs.flow"),
     mean.col = 1:4, mean.lwd = 1, mean.smooth = TRUE,
     qnts.col = 1:4, qnts.alpha = 0.25, qnts.smooth = TRUE,
     ylim = c(0, 3), legend = TRUE)

# Average across simulations at beginning, middle, end
df <- as.data.frame(sim)
df[c(2, 100, 500), ]
