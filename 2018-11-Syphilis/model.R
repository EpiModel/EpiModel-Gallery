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
n <- 1000
nw <- network.initialize(n, directed = FALSE)

# Define the formation model: edges + isolates (number with degree of 0)
formation = ~edges + isolates

# Input the appropriate target statistics for each term
target.stats <- c(300, 480)

# Parameterize the dissolution model
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 10)
coef.diss

# Fit the model
est <- netest(nw, formation, target.stats, coef.diss)

# Model diagnostics
dx <- netdx(est, nsims = nsims, ncores = ncores, nsteps = nsteps,
            nwstats.formula = ~edges + isolates + degree(0:3))
print(dx)
plot(dx)


# Epidemic model simulation -----------------------------------------------

# Model parameters
param <- param.net(inf.prob1 = 0.18,
                   inf.prob2 = 0.09,
                   act.rate = 2,
                   ipr.rate = 1/4,
                   prse.rate = 1/9,
                   seel.rate = 1/17,
                   elll.rate = 1/22,
                   llter.rate = 1/1508,
                   pri.sym = 0.205,
                   sec.sym = 0.106,
                   early.trt = 0.8,
                   late.trt = 1.0,
                   scr.rate = 1/52)

# Initial conditions
init <- init.net(i.num = 10)

# Read in the module functions
if (interactive()) {
  source("2018-11-Syphilis/module-fx.R", echo = TRUE)
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
                       tnt.FUN = tnt)

# Run the network model simulation with netsim
sim <- netsim(est, param, init, control)
print(sim)

# Examine data from simulations
df <- as.data.frame(sim)
dplyr::filter(df, sim == 1 & time %in% c(2, 25, 50))

## Plot outcomes

# All compartments
par(mar = c(3,3,1,1), mgp = c(2,1,0))
plot(sim,y = c("s.num", "i.num", "inc.num", "pr.num",
               "se.num", "el.num", "ll.num", "ter.num"),
     mean.col = 1:8, mean.lwd = 1, mean.smooth = FALSE,
     qnts = 1, qnts.col = 1:8, qnts.alpha = 0.25, qnts.smooth = FALSE,
     legend = TRUE)

# SI compartment counts
par(mar = c(3,3,1,1), mgp = c(2,1,0))
plot(sim,y = c("s.num", "i.num"), mean.col = 1:2, mean.lwd = 1,
     mean.smooth = FALSE, qnts = 1, qnts.col = c(1:2), qnts.alpha = 0.25,
     qnts.smooth = FALSE, legend = TRUE)

# Syphilis stage compartment counts
par(mar = c(3,3,1,1), mgp = c(2,1,0))
plot(sim,y = c("inc.num", "pr.num", "se.num",
               "el.num", "ll.num", "ter.num"),
     mean.col = 3:8, mean.lwd = 1, mean.smooth = FALSE,
     qnts = 1, qnts.col = 3:8, qnts.alpha = 0.25, qnts.smooth = FALSE)
legend("topleft", c("inc.num","pr.num","se.num",
                    "el.num", "ll.num", "ter.num"),
       col = 3:8, lty = 1, cex = 0.8)

# Incidence/compartment flow rates
plot(sim, y = c("si.flow", "ipr.flow", "prse.flow",
                "seel.flow", "elll.flow", "llter.flow"),
     mean.col = 1:6, mean.lwd = 1, mean.smooth = TRUE,
     qnts.col = 1:6, qnts.alpha = 0.25, qnts.smooth = TRUE,
     ylim = c(0, 10), legend = TRUE)

# Plot screening rate
plot(sim, y = "scr.flow", mean.col = 1, mean.lwd = 1, mean.smooth = TRUE,
     qnts.col = 1, qnts.alpha = 0.25, qnts.smooth = TRUE, legend = TRUE)

# Duration spent in each stage of infection
plot(sim, y = c("syph.dur", "syph2.dur", "syph3.dur", "syph4.dur"),
     mean.col = 1:4, mean.lwd = 1, mean.smooth = TRUE,
     qnts.col = 1:4, qnts.alpha = 0.25, qnts.smooth = TRUE, ylim = c(0,40),
     legend = TRUE)

# Prevalence of symptoms
plot(sim, y = "sym.num",
     mean.col = 1, mean.lwd = 1, mean.smooth = TRUE, qnts.col = 1,
     qnts.alpha = 0.25, qnts.smooth = TRUE, legend = TRUE)
