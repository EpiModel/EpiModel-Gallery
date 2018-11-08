##
## Information Diffusion Model in a Social Network: Adding a minimum number of degree to S in an SI Model
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##
## Author: Yuan Zhao
## Date: September 2018
##

# Load EpiModel
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
target.stats <- c(300, 30)

# Parameterize the dissolution model
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 10)
coef.diss

# Fit the model
est <- netest(nw, formation, target.stats, coef.diss)

# Model diagnostics
dx <- netdx(est, nsims = 8, ncores = 8, nsteps = 300,
            nwstats.formula = ~edges + isolates + degree(0:5))
print(dx)
plot(dx, plots.joined = FALSE, qnts.alpha = 0.8)


# Epidemic model simulation -----------------------------------------------

# Model parameters
param <- param.net(inf.prob = 0.5, act.rate = 2,
                    min.degree=3)

# Initial conditions
init <- init.net(i.num = 100)

# Read in the module functions
source("module-fx.R", echo = TRUE)

# Control settings
control <- control.net(type="SI",
                       nsteps = 300,
                       nsims = 8,
                       ncores = 4,
                       infection.FUN = infect_mod
)

# Run the network model simulation with netsim
set.seed(123456)
sim <- netsim(est, param, init, control)
print(sim)

# Plot outcomes
par(mar = c(3,3,1,1), mgp = c(2,1,0))
plot(sim,
     mean.col = 1:4, mean.lwd = 1, mean.smooth = FALSE,
     qnts = 1, qnts.col = 1:4, qnts.alpha = 0.25, qnts.smooth = FALSE,
     legend = TRUE)

plot(sim, y = c("si.flow"),
     mean.col = 1:4, mean.lwd = 1, mean.smooth = TRUE,
     qnts.col = 1:4, qnts.alpha = 0.25, qnts.smooth = TRUE,
     ylim = c(0, 3), legend = TRUE)

# Average across simulations at beginning, middle, 
df <- as.data.frame(sim)
df[c(2, 50, 100), ]



# Plot network
par(mar = c(3,3,1,1), mgp = c(2,1,0))
#plot(sim,type="network",at=1,sims="mean",
     #col.status=TRUE, main="Prevalence at t1")
plot(sim,type="network",at=2,sims="mean",
     col.status=TRUE, main="Prevalence at t2")
#plot(sim,type="network",at=3,sims="mean",
     #col.status=TRUE, main="Prevalence at t3")
#plot(sim,type="network",at=10,sims="mean",
     #col.status=TRUE, main="Prevalence at t10")
plot(sim,type="network",at=50,sims="mean",
     col.status=TRUE, main="Prevalence at t50")


##Scenario 2 the infection probility as the function of discordant degree
# Model parameters
param <- param.net(act.rate = 2, ei.rate = 0.01, ir.rate = 0.01, 
                   log_a=100, log_b=0.7)

# Initial conditions
init <- init.net(i.num = 150)

control <- control.net(nsteps = 300,
                       nsims = 8,
                       ncores = 1,
                       infection.FUN = infect_newmod
)

# Run the network model simulation with netsim
sim <- netsim(est, param, init, control)

print(sim)

# Plot outcomes
par(mar = c(3,3,1,1), mgp = c(2,1,0))
plot(sim,
     mean.col = 1:4, mean.lwd = 1, mean.smooth = FALSE,
     qnts = 1, qnts.col = 1:4, qnts.alpha = 0.25, qnts.smooth = FALSE,
     legend = TRUE)

plot(sim, y = c("si.flow"),
     mean.col = 1:4, mean.lwd = 1, mean.smooth = TRUE,
     qnts.col = 1:4, qnts.alpha = 0.25, qnts.smooth = TRUE,
     ylim = c(0, 3), legend = TRUE)

# Average across simulations at beginning, middle, 
df <- as.data.frame(sim)
df[c(3, 4, 5), ]
df[c(2, 100, 300), ]



# Plot network
par(mar = c(3,3,1,1), mgp = c(2,1,0))
plot(sim,type="network",at=2,sims="mean",
     col.status=TRUE, main="Prevalence at t2")
plot(sim,type="network",at=10,sims="mean",
     col.status=TRUE, main="Prevalence at t10")
plot(sim,type="network",at=100,sims="mean",
     col.status=TRUE, main="Prevalence at t100")