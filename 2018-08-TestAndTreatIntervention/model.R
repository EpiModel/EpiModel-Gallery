
##
## Test and Treat Intervention for an SIS Epidemic
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##
## Author: Samuel M. Jenness (Emory University)
## Date: August 2018
##

# Load EpiModel
library(EpiModel)
rm(list = ls())
options(error = stop)

# Network model estimation ------------------------------------------------

# Initialize the network
nw <- network.initialize(n = 500, directed = FALSE)

# Define the formation model: edges,
#                             number concurrent (degree > 1),
#                             number with degree 4+
formation <- ~edges + concurrent + degrange(from = 4)

# Input the appropriate target statistics for each term
target.stats <- c(175, 110, 0)

# Parameterize the dissolution model
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 50)
coef.diss

# Fit the model
est <- netest(nw, formation, target.stats, coef.diss)

# Model diagnostics
dx <- netdx(est, nsims = 8, ncores = 8, nsteps = 500)
print(dx)
plot(dx, plots.joined = FALSE)


# Epidemic model simulation -----------------------------------------------

# Parameterizing an SIS epidemic
param <- param.net(inf.prob = 0.4, act.rate = 2,
                   rec.rate = 0.05, rec.rate.tx = 0.5,
                   test.rate = 0.1, test.dur = 2)

# Initial conditions
init <- init.net(i.num = 10)

# Read in the module functions
source("module-fx.R", echo = TRUE)

# Control settings
control <- control.net(nsims = 8,
                       ncores = 8,
                       nsteps = 500,
                       recovery.FUN = recov,
                       tnt.FUN = tnt)

# Run the network model simulation with netsim
sim <- netsim(est, param, init, control)
print(sim)

# Plot outcomes
par(mar = c(3,3,2,1), mgp = c(2,1,0))
plot(sim)
plot(sim, y = c("si.flow", "is.flow"), legend = TRUE)
plot(sim, y = c("nTest", "nReset"), legend = TRUE)

# Average across simulations at beginning, middle, end
df <- as.data.frame(sim)
df[c(2, 100, 500), ]

plot(df2)
