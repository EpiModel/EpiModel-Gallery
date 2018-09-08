
##
## Modeling Two Competing Strains in an SIS Epidemic
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##
## Author: Steven M. Goodreau (University of Washington)
## Date: September 2018
##

# Load EpiModel
library(EpiModel)
rm(list = ls())


# Network model estimation ------------------------------------------------

# Initialize the network
nw <- network.initialize(n = 500, directed = FALSE)

# Define the formation model: edges,
#                             number concurrent (degree > 1),
#                             number with degree 4+
#formation <- ~edges + concurrent + degrange(from = 4)
formation <- ~edges

# Input the appropriate target statistics for each term
#target.stats <- c(175, 110, 0)
target.stats <- c(175)

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
param <- param.net(inf.prob = 0.5, inf.prob.st2 = 0.6,
                   pct.st2 = 0.5, 
                   act.rate = 2, rec.rate = 0.1
                  )

# Initial conditions
init <- init.net(i.num = 50)

# Read in the module functions
source("module-fx.R", echo = TRUE)

# Control settings
control <- control.net(type = "SIS", 
                       nsims = 1,
                       ncores = 1,
                       nsteps = 500,
                       #recovery.FUN = recov,
                       infection.FUN = infection.2strains)

# Run the network model simulation with netsim
sim <- netsim(est, param, init, control)
print(sim)

sim <- mutate_epi(sim, st1.prev = i.num.st1 / (i.num.st1 + i.num.st2))
sim <- mutate_epi(sim, st2.prev = i.num.st2 / (i.num.st1 + i.num.st2))
plot(sim, y=c("st1.prev", "st2.prev"), sim.lines = TRUE, 
     sim.col=c('purple', 'green'), mean.line = FALSE)

# Plot outcomes
par(mar = c(3,3,2,1), mgp = c(2,1,0))
plot(sim)
plot(sim, y = c("si.flow", "is.flow"), legend = TRUE)
plot(sim, y = c("nTest", "nReset"), legend = TRUE)

# Average across simulations at beginning, middle, end
df <- as.data.frame(sim)
df[c(2, 100, 500), ]
