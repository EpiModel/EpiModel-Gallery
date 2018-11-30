##
## SEIR Model extension: Syphilis Progression Model
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##
## Authors: Samuel M. Jenness, Yuan Zhao
## Date: November 2018
##

## Load EpiModel
suppressMessages(library(EpiModel))

# Standard Gallery unit test lines
rm(list = ls())
eval(parse(text = print(commandArgs(TRUE)[1])))

# Network model estimation ------------------------------------------------

# Initialize the network
nw <- network.initialize(500, directed = FALSE)

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
dx <- netdx(est, nsims = 8, ncores = 8, nsteps = 500,
            nwstats.formula = ~edges + isolates + degree(0:5))
print(dx)
plot(dx, plots.joined = FALSE, qnts.alpha = 0.8)


# Epidemic model simulation -----------------------------------------------

# Model parameters
param <- param.net(inf.prob1 = 0.18, inf.prob2 = 0.09, act.rate = 2,
                   ipr.rate = 1/4, prse.rate = 1/9, seel.rate = 1/17,
                   elll.rate = 1/22,llter.rate= 1/1508)

# Initial conditions
init <- init.net(i.num=50)

# Read in the module functions
source("module-fx.R", echo = TRUE)

# Control settings
control <- control.net(nsteps = 600,
                       nsims = 4,
                       ncores = 1,
                       infection.FUN = infect,
                       progress.FUN = progress,
                       recovery.FUN = NULL)

# Run the network model simulation with netsim
sim <- netsim(est, param, init, control)
print(sim)

# Plot outcomes
par(mar = c(3,3,1,1), mgp = c(2,1,0))
plot(sim,
     mean.col = 1:8, mean.lwd = 1, mean.smooth = FALSE,
     qnts = 1, qnts.col = 1:8, qnts.alpha = 0.25, qnts.smooth = FALSE,
     legend = TRUE)

plot(sim, y = c("si.flow", "ipr.flow", "prse.flow","seel.flow", "elll.flow", "llter.flow"),
     mean.col = 1:6, mean.lwd = 1, mean.smooth = TRUE,
     qnts.col = 1:6, qnts.alpha = 0.25, qnts.smooth = TRUE,
     ylim=c(0,15),legend = TRUE)

# Average across simulations at beginning, middle, end
df <- as.data.frame(sim)
df[c(2, 100, 500), ]



