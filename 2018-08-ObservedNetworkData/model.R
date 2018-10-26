
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


# Epidemic model simulation -----------------------------------------------

# Epidemic model parameters
param <- param.net(inf.prob = 0.5)

# Initial conditions
init <- init.net(i.num = 10)

# Read in the module functions
source("module-fx.R", echo = TRUE)

# Control settings (must be link nsteps to number of observed time steps in network)
control <- control.net(type = "SI",
                       nsteps = 100,
                       nsims = 4,
                       ncores = 4,
                       initialize.FUN = new_init_mod,
                       infection.FUN = new_infect_mod,
                       module.order = c("infection.FUN", "get_prev.FUN"),
                       skip.check = TRUE,
                       save.nwstats = FALSE,
                       save.network = FALSE,
                       verbose.int = 0)

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
control <- control.net(type = "SI",
                       nsteps = 200,
                       nsims = 4,
                       ncores = 4,
                       initialize.FUN = new_init_mod,
                       infection.FUN = new_infect_mod,
                       module.order = c("infection.FUN", "get_prev.FUN"),
                       skip.check = TRUE,
                       save.nwstats = FALSE,
                       save.network = FALSE,
                       verbose.int = 0)

# However, you won't get meaningful results past the observations
sim <- netsim(nw, param, init, control)
plot(sim, main = "Prevalence")
plot(sim, y = "si.flow", main = "Incidence")
