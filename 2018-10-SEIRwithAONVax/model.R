##
## SEIR Model with Vital Dynamics and an All or Nothing Vaccine Implementation
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##
## Authors: Connor M. Van Meter
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
n <- 100
nw <- network_initialize(n)

# Define the formation model: edges
formation = ~edges

# Input the appropriate target statistics for each term
mean_degree <- 0.8
edges <- mean_degree * (n/2)

# Input the appropriate target statistics for each term
target.stats <- c(edges)

#Set mortality rate
mr_rate = 0.008

# Parameterize the dissolution model
coef.diss <- dissolution_coefs(dissolution = ~offset(edges),
                               duration = 20, d.rate = mr_rate)
coef.diss

# Fit the model
est <- netest(nw, formation, target.stats, coef.diss)

# Model diagnostics
dx <- netdx(est, nsims = nsims, ncores = ncores, nsteps = nsteps)

print(dx)
plot(dx)


# Epidemic model simulation -----------------------------------------------

# Model parameters
param <- param.net(inf.prob = 0.5,
                   birth.rate = 0.01,
                   mortality.rate = mr_rate,
                   mortality.disease.mult = 2,
                   act.rate = 1,
                   ei.rate = 0.05,
                   ir.rate = 0.05,
                   vaccination.rate.initialization = 0.05,
                   protection.rate.initialization = 0.8,
                   vaccination.rate.progression = 0.05,
                   protection.rate.progression = 0.8,
                   vaccination.rate.births = 0.6,
                   protection.rate.births = 0.8)

# Initial conditions
init <- init.net(i.num = 20)

# Read in the module functions
if (interactive()) {
  source("2018-10-SEIRwithAONVax/module-fx.R", echo = TRUE)
} else {
  source("module-fx.R", echo = TRUE)
}

# Control settings
control <- control.net(type = NULL,
                       nsteps = nsteps,
                       nsims = nsims,
                       ncores = ncores,
                       infection.FUN = infect,
                       progress.FUN = progress,
                       departures.FUN = dfunc,
                       arrivals.FUN = bfunc,
                       resimulate.network = TRUE,
                       verbose = TRUE)

# Run the network model simulation with netsim
sim <- netsim(est, param, init, control)
print(sim)


##################################################################

# Examine the data from the simulation
df <- as.data.frame(sim)
df

#Epidemic plot of SEIR-V compartment counts, entrances, and exits over simulation
par(mar = c(3,3,1,1), mgp = c(2,1,0), mfrow = c(1,1))
plot(sim, y = c("s.num","e.num","i.num","r.num", "v.num", "b.num", "d.num", "num"),
     mean.col = 1:8, mean.lwd = 1, mean.smooth = FALSE,
     qnts = 1, qnts.col = 1:8, qnts.alpha = 0.25, qnts.smooth = FALSE,
     legend = TRUE)

#Cumulative Incidence and Prevalence Plots of the Vaccine Model
sim <- mutate_epi(sim, ci = se.flow / s.num, prev = e.num / num)
plot(sim, y = c("ci", "prev"), mean.lwd = 1, mean.smooth = TRUE, legend = TRUE)


##################################################################
#VACCINE MODEL COMPARISON
##################################################################

# Update one or more vaccine parameters, simulate the new model,
# and compare with the original model

# Model parameters
param <- param.net(inf.prob = 0.5,
                   birth.rate = 0.01,
                   mortality.rate = mr_rate,
                   mortality.disease.mult = 2,
                   act.rate = 1,
                   ei.rate = 0.05,
                   ir.rate = 0.05,
                   vaccination.rate.initialization = 0.05,
                   protection.rate.initialization = 0.3,
                   vaccination.rate.progression = 0.05,
                   protection.rate.progression = 0.3,
                   vaccination.rate.births = 0.2,
                   protection.rate.births = 0.3)

# Initial conditions
init <- init.net(i.num = 20)

# Control settings
control <- control.net(type = NULL,
                       nsteps = nsteps,
                       nsims = nsims,
                       ncores = ncores,
                       infection.FUN = infect,
                       progress.FUN = progress,
                       recovery.FUN = NULL,
                       arrivals.FUN = bfunc,
                       departures.FUN = dfunc,
                       prevalence.FUN = prevalence.net,
                       resimulate.network = TRUE,
                       verbose = TRUE)

# Run the network model simulation with netsim
sim2 <- netsim(est, param, init, control)
print(sim2)

##Compare incidence and prevalence of simulation 1 to simulation 2

#Calculate cumulative incidence and prevalence of simulation 2
sim2 <- mutate_epi(sim2, ci2 = se.flow / s.num, prev2 = e.num / num)

par(mfrow = c(1,1))
plot(sim, y = c("ci", "prev"), mean.lwd = 1, mean.smooth = TRUE, legend = TRUE)
plot(sim2, y = c("ci2", "prev2"), mean.lwd = 1, mean.smooth = TRUE, add = TRUE,
     mean.col = c("steelblue", "firebrick"), qnts.col = c("steelblue", "firebrick"))
