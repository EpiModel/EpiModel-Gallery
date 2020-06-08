
## Modeling Two Competing Strains in an SIS Epidemic
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##
## Author: Steven M. Goodreau (University of Washington)
## Date: September 2018
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
nw <- network.initialize(n = 1000, directed = FALSE)

# Define two different formation models and target stats associated with them
#   Model 1:                            # All ties equally likely, no constraint on degree
formation.mod1 <- ~edges
target.stats.mod1 <- 300
#   Model 2:                            # Concurrent term w/ target stat 0 means nobody
formation.mod2 <- ~edges + concurrent   # allowed to have concurrent ties - but still
target.stats.mod2 <- c(300, 0)          # same overall number of ties as Model 1

# Parameterize the dissolution model (same for Models 1 and 2)
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 50)
coef.diss

# Fit the models
est.mod1 <- netest(nw, formation.mod1, target.stats.mod1, coef.diss)
est.mod2 <- netest(nw, formation.mod2, target.stats.mod2, coef.diss)

# Model diagnostics
dx.mod1 <- netdx(est.mod1, nsims = nsims, ncores = ncores, nsteps = nsteps,
                 set.control.ergm = control.simulate.ergm(MCMC.burnin = 1e5))
print(dx.mod1)
plot(dx.mod1, plots.joined = FALSE)

dx.mod2 <- netdx(est.mod2, nsims = nsims, ncores = ncores, nsteps = nsteps,
                 set.control.ergm = control.simulate.ergm(MCMC.burnin = 1e5))
print(dx.mod2)
plot(dx.mod2, plots.joined = FALSE)


# Epidemic model simulation -----------------------------------------------

# Parameterizing an SIS epidemic
# Note: strain 1 is highly infectious but short-lived;
#       strain 2 has much lower infection but longer duration
param <- param.net(inf.prob = 0.5, inf.prob.st2 = 0.01,
                   rec.rate = 0.05, rec.rate.st2 = 0.005,
                   pct.st2 = 0.5,
                   act.rate = 2)

# Initial conditions
init <- init.net(i.num = 50)

# Read in the module functions
if (interactive()) {
  source("2018-09-CompetingStrains/module-fx.R", echo = TRUE)
} else {
  source("module-fx.R")
}

# Control settings
control <- control.net(type = NULL,
                       nsims = nsims,
                       ncores = ncores,
                       nsteps = nsteps,
                       infection.FUN = infection.2strains,
                       recovery.FUN = recov.2strains)

# Run the network model simulations with netsim
sim.mod1 <- netsim(est.mod1, param, init, control)
sim.mod2 <- netsim(est.mod2, param, init, control)


# Plotting results -----------------------------------------------

## In model 1, strain 1 dominates strain 2.
## In model 2, strain 1 goes extinct while strain 2 persists.
## The only difference between the two models was that one
## enforced a monogamy rule and the other did not.

par(mfrow = c(1, 2), mar = c(3,3,2,1), mgp = c(2,1,0))
plot(sim.mod1, y = c("i.num.st1", "i.num.st2"),
     sim.lines = TRUE, mean.line = TRUE, mean.lwd = 2,
     qnts = FALSE, main = "Model 1")

plot(sim.mod2, y = c("i.num.st1", "i.num.st2"),
     sim.lines = TRUE, mean.line = TRUE, mean.lwd = 2,
     qnts = FALSE, main = "Model 2")


## Probing further  --------------------------------------------------------

# At what level of concurrency does the cross-over point occur?

if (interactive()) {

# Check  how many nodes had concurrent ties on average in model 1
dx.mod1a <- netdx(est.mod1, nsims = 10, ncores = 10, nsteps = 100,
                  set.control.ergm = control.simulate.ergm(MCMC.burnin = 1e5),
                  nwstats.formula = ~edges+concurrent)
dx.mod1a             # Roughly 120

# Define model, and set up a vector of concurrency levels to explore
formation.mod3 <- ~edges + concurrent
concties.mod3 <- seq(0,120,10)
est.mod3 <- list()
sim.mod3 <- list()

# Run models
# Warning: this loop can take 30+ minutes to run
for (i in 1:length(concties.mod3)) {
  target.stats.mod3 <- c(300, concties.mod3[i])
  est.mod3[[i]] <- suppressMessages(netest(nw, formation.mod3, target.stats.mod3, coef.diss))
  sim.mod3[[i]] <- netsim(est.mod3[[i]], param, init, control)
  cat("\n ConcTies =", concties.mod3[i], "complete ...")
}

# Process output
i.num.st1.final.mod3 <- sapply(1:13, function(x) rowMeans(sim.mod3[[x]]$epi$i.num.st1)[nsteps])
i.num.st2.final.mod3 <- sapply(1:13, function(x) rowMeans(sim.mod3[[x]]$epi$i.num.st2)[nsteps])

# Plot results
par(mfrow = c(1, 1))
plot(concties.mod3, i.num.st1.final.mod3, type = "b", col = "purple",
     xlab = "exp. # of persons with concurrent ties",
     ylab = "prevalence of strains at time step 1000", ylim = c(0, 500))
points(concties.mod3, i.num.st2.final.mod3, type = "b", col = "green")
legend(0, 500, legend = c("strain 1", "strain 2"), col = c("purple", "green"), lwd = 2)

}
