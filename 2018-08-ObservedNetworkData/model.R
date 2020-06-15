
##
## Modeling Epidemics over Observed Networks
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


# Import Observed Network Data --------------------------------------------

# Use dynamic data from networkDynamicData package
library(networkDynamicData)

# Actually, this is a simulated dataset, but let's pretend it is observed
data(concurrencyComparisonNets)
nw <- base

# Examine the structure of a networkDynamic class object to the correct data form
print(nw)
head(as.data.frame(nw), 50)

# Removing the "observed" disease status attribute, as we'll be simulating that
nw <- delete.vertex.attribute(nw, "status.active")

if (interactive()) {
  nsims <- 5
  ncores <- 5
  nsteps <- 100
} else {
  nsims <- 1
  ncores <- 1
  nsteps <- 50
}

# Example 1: Epidemic Model Simulation ------------------------------------

# Epidemic model parameters
param <- param.net(inf.prob = 0.5)

# Initial conditions
init <- init.net(i.num = 10)

# Read in the module functions
if (interactive()) {
  source("2018-08-ObservedNetworkData/module-fx.R", echo = TRUE)
} else {
  source("module-fx.R")
}

# Control settings (must be link nsteps to number of observed time steps in network)
control <- control.net(type = NULL,
                       nsteps = nsteps,
                       nsims = nsims,
                       ncores = ncores,
                       initialize.FUN = new_init_mod,
                       infection.FUN = new_infect_mod,
                       prevalence.FUN = prevalence.net,
                       resimulate.network = FALSE,
                       skip.check = TRUE,
                       tergmLite = FALSE,
                       verbose = FALSE)

# Run the network model simulation with netsim
sim <- netsim(nw, param, init, control)
print(sim)

# Plot outcomes
par(mfrow = c(1,2), mar = c(3,3,1,1), mgp = c(2,1,0))
plot(sim, main = "Prevalence")
plot(sim, y = "si.flow", main = "Incidence")

# Examine the data
df <- as.data.frame(sim)
head(df, 25)


# Caveats -----------------------------------------------------------------

# Dyads in observed nD object remain active after the last observation
head(get.dyads.active(nw, at = 1), 10)
head(get.dyads.active(nw, at = 100), 10)
head(get.dyads.active(nw, at = 200), 10)
head(get.dyads.active(nw, at = Inf), 10)

# Nothing prevents you from simulating past the observations
control <- control.net(type = NULL,
                       nsteps = 200,
                       nsims = nsims,
                       ncores = ncores,
                       initialize.FUN = new_init_mod,
                       infection.FUN = new_infect_mod,
                       prevalence.FUN = prevalence.net,
                       resimulate.network = FALSE,
                       skip.check = TRUE,
                       tergmLite = FALSE,
                       verbose = FALSE)

# However, you won't get meaningful results past the observations
sim <- netsim(nw, param, init, control)
plot(sim, main = "Prevalence")
plot(sim, y = "si.flow", main = "Incidence")



# Example 2: Adding Networking Tracking and Time-Varying Risk -------------

# Epidemic model parameters, adding a two-disease-stage model, where the infection
# probability varies by stage and the first stage lasts 5 time steps
param <- param.net(inf.prob.stage1 = 0.05,
                   inf.prob.stage2 = 0.15,
                   dur.stage1 = 5)

# Initial conditions
init <- init.net(i.num = 10)

# Control settings (must be link nsteps to number of observed time steps in network)
control <- control.net(type = NULL,
                       nsteps = nsteps,
                       nsims = nsims,
                       ncores = ncores,
                       initialize.FUN = new_init_mod2,
                       infection.FUN = new_infect_mod2,
                       prevalence.FUN = prevalence.net,
                       resimulate.network = FALSE,
                       skip.check = TRUE,
                       tergmLite = FALSE,
                       verbose = FALSE)

# Run the network model simulation with netsim
sim <- netsim(nw, param, init, control)
print(sim)

# Plot outcomes
par(mfrow = c(1,2), mar = c(3,3,1,1), mgp = c(2,1,0))
plot(sim, main = "Prevalence")
plot(sim, y = "si.flow", main = "Incidence")

# Plot network at various time steps
par(mfrow = c(1,2), mar = c(1,1,1,1))
plot(sim, type = "network", col.status = TRUE, at = 2, sims = 1)
plot(sim, type = "network", col.status = TRUE, at = nsteps, sims = 1)

# Extract individual-level attributes over time
nwd <- get_network(sim, 1)
head(get.vertex.attribute.active(nwd, "testatus", at = 1), 25)
head(get.vertex.attribute.active(nwd, "testatus", at = nsteps), 25)

# Examine the data
df <- as.data.frame(sim)
head(df, 25)

# Mean statistics
df_mean <- as.data.frame(sim, out = "mean")
head(df_mean, 25)

# Cumulative incidence
sum(df_mean$si.flow, na.rm = TRUE)
