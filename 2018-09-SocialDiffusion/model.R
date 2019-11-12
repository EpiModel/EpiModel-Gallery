
##
## Modeling Social Diffusion in a Social Network
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##
## Author: Yuan Zhao
## Date: November 2018
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

# Network model estimation ------------------------------------------------

# Initialize the network
nw <- network.initialize(500, directed = FALSE)

# Define the formation model: edges + isolates (number with degree of 0)
formation = ~edges + isolates

# Input the appropriate target statistics for each term
target.stats <- c(450, 30)

# Parameterize the dissolution model
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 10)
coef.diss

# Fit the model
est <- netest(nw, formation, target.stats, coef.diss)

# Model diagnostics
dx <- netdx(est, nsims = nsims, ncores = ncores, nsteps = nsteps,
            nwstats.formula = ~edges + isolates + degree(0:5))
print(dx)
plot(dx, plots.joined = FALSE, qnts.alpha = 0.8)


# Scenario 1 --------------------------------------------------------------

# Model parameters
param <- param.net(inf.prob = 0.5,
                   act.rate = 2,
                   min.degree = 3)

# Initial conditions
init <- init.net(i.num = 100)

# Read in the module functions
if (interactive()) {
  source("2018-09-SocialDiffusion/module-fx.R", echo = TRUE)
} else {
  source("module-fx.R", echo = TRUE)
}

# Control settings
control <- control.net(type = "SI",
                       nsteps = nsteps,
                       nsims = nsims,
                       ncores = ncores,
                       infection.FUN = diffuse_mod)


# Run the network model simulation with netsim
sim <- netsim(est, param, init, control)
print(sim)

# Plot outcomes
par(mar = c(3,3,1,1), mgp = c(2,1,0))
plot(sim,
     mean.col = 1:2, mean.lwd = 1, mean.smooth = FALSE,
     qnts = 1, qnts.col = 1:2, qnts.alpha = 0.25, qnts.smooth = FALSE,
     legend = TRUE)

plot(sim, y = "si.flow",
     mean.col = 1, mean.lwd = 1, mean.smooth = TRUE,
     qnts.col = 1, qnts.alpha = 0.25, qnts.smooth = TRUE,
     ylim = c(0, 6), legend = TRUE)

# Average across simulations at beginning, middle, and end
df <- as.data.frame(sim)
dplyr::filter(df, sim == 1 & time %in% c(2, 50, 100, 290))


# Scenario 2 --------------------------------------------------------------

# Model parameters
param <- param.net(act.rate = 2,
                   beta0 = -7,
                   beta1 = 0.5)

# Initial conditions
init <- init.net(i.num = 100)

# Control settings
control <- control.net(nsteps = nsteps,
                       nsims = nsims,
                       ncores = ncores,
                       infection.FUN = diffuse_mod2)

# Run the network model simulation with netsim
sim <- netsim(est, param, init, control)
print(sim)

# Plot outcomes
par(mar = c(3,3,1,1), mgp = c(2,1,0))
plot(sim,
     mean.col = 1:2, mean.lwd = 1, mean.smooth = FALSE,
     qnts = 1, qnts.col = 1:2, qnts.alpha = 0.25, qnts.smooth = FALSE,
     legend = TRUE)

plot(sim, y = "si.flow",
     mean.col = 1, mean.lwd = 1, mean.smooth = TRUE,
     qnts.col = 1, qnts.alpha = 0.25, qnts.smooth = TRUE,
     ylim = c(0, 2), legend = TRUE)

# Average across simulations at beginning, middle, and the end
df <- as.data.frame(sim)
dplyr::filter(df, sim == 1 & time %in% c(2, 50, 100, 290))
