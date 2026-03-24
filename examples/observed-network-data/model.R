
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

if (interactive()) {
  nsims <- 5
  ncores <- 5
  nsteps <- 100
} else {
  nsims <- 1
  ncores <- 1
  nsteps <- 50
}


# 1. Import Observed Network Data -------------------------------------------

# The standard approach in EpiModel uses egocentrically observed network data
# to fit a temporal ERGM, then simulates from that model fit. This example
# demonstrates an alternative: using a dynamic network census -- an observed
# network with all nodes and edges recorded over discrete time steps.

# Load observed network from the networkDynamicData package
library(networkDynamicData)
data(concurrencyComparisonNets)
nw <- base
print(nw)

# Remove the existing disease status attribute from the observed data,
# since we will be simulating our own epidemic over this network
nw <- network::delete.vertex.attribute(nw, "status.active")

# Load custom module functions
source("examples/observed-network-data/module-fx.R")


# 2. Example 1: Basic SI Epidemic -------------------------------------------

# Simple SI model: constant per-act transmission probability of 50%
param <- param.net(inf.prob = 0.5, act.rate = 1)

# Initial conditions: 10 infected nodes to seed the epidemic
init <- init.net(i.num = 10)

# Control settings: custom modules for observed network.
# resimulate.network = FALSE because we use the fixed observed network
# rather than resimulating from an ERGM fit each timestep.
# save.nwstats = FALSE because there is no ERGM fit to compute stats from.
# nsteps must link to the number of observed time steps in the network.
control <- control.net(
  type = NULL,
  nsteps = nsteps,
  nsims = nsims,
  ncores = ncores,
  initialize.FUN = init_obsnw,
  infection.FUN = infect_obsnw,
  prevalence.FUN = prevalence.net,
  resimulate.network = FALSE,
  save.nwstats = FALSE,
  verbose = FALSE
)

# Run the network model simulation with netsim
sim <- netsim(nw, param, init, control)
print(sim)

# Examine the data
df <- as.data.frame(sim)
head(df, 25)


# 3. Caveats: Simulating Past the Observation Window -------------------------

# The observed network has edges active over a finite window (~100 timesteps).
# Dyads that are active at the last observed time remain active indefinitely
# (networkDynamic convention), so nothing prevents simulating past the
# observations -- but the results become meaningless.

if (interactive()) {

  # Dyads remain active beyond the observation window
  head(get.dyads.active(nw, at = 1), 10)
  head(get.dyads.active(nw, at = 100), 10)
  head(get.dyads.active(nw, at = 200), 10)
  head(get.dyads.active(nw, at = Inf), 10)

  # Simulate deliberately past the observation window (200 steps on a
  # ~100-step network) to demonstrate the problem
  control_caveat <- control.net(
    type = NULL,
    nsteps = 200,
    nsims = nsims,
    ncores = ncores,
    initialize.FUN = init_obsnw,
    infection.FUN = infect_obsnw,
    prevalence.FUN = prevalence.net,
    resimulate.network = FALSE,
    save.nwstats = FALSE,
    verbose = FALSE
  )

  sim_caveat <- netsim(nw, param, init, control_caveat)

  # Past the observation window, prevalence plateaus and incidence drops to
  # zero because the frozen edge set no longer reflects real dynamics
  par(mfrow = c(1, 2), mar = c(3, 3, 2, 1), mgp = c(2, 1, 0))
  plot(sim_caveat, y = "i.num",
       main = "Prevalence (Past Observations)",
       ylab = "Infected Count", xlab = "Time Step",
       mean.col = "firebrick", mean.lwd = 2, mean.smooth = FALSE,
       qnts.col = "firebrick", qnts.alpha = 0.2, qnts.smooth = FALSE,
       legend = FALSE)
  plot(sim_caveat, y = "si.flow",
       main = "Incidence (Past Observations)",
       ylab = "New Infections", xlab = "Time Step",
       mean.col = "steelblue", mean.lwd = 2, mean.smooth = FALSE,
       qnts.col = "steelblue", qnts.alpha = 0.2, qnts.smooth = FALSE,
       legend = FALSE)

}


# 4. Example 2: Time-Varying Transmission ------------------------------------

# Two-stage disease model where transmission probability depends on
# infection duration. Primary stage (first 5 timesteps) has a low
# transmission probability (5% per act), while secondary stage has a
# higher probability (15% per act).
#
# track.nw.attr = TRUE stores time-varying disease status on the network,
# enabling visualization with plot(sim, type = "network", col.status = TRUE).
param <- param.net(
  inf.prob.stage1 = 0.05,
  inf.prob.stage2 = 0.15,
  dur.stage1 = 5,
  act.rate = 1,
  track.nw.attr = TRUE
)

# Initial conditions
init <- init.net(i.num = 10)

# Control settings (same modules -- behavior changes via parameters)
control <- control.net(
  type = NULL,
  nsteps = nsteps,
  nsims = nsims,
  ncores = ncores,
  initialize.FUN = init_obsnw,
  infection.FUN = infect_obsnw,
  prevalence.FUN = prevalence.net,
  resimulate.network = FALSE,
  save.nwstats = FALSE,
  verbose = FALSE
)

# Run the network model simulation with netsim
sim2 <- netsim(nw, param, init, control)
print(sim2)

# Plot network at early and late time steps, colored by disease status
par(mfrow = c(1, 2), mar = c(1, 1, 1, 1))
plot(sim2, type = "network", col.status = TRUE, at = 2, sims = 1)
plot(sim2, type = "network", col.status = TRUE, at = nsteps, sims = 1)

# Extract individual-level disease status from the network at different times
nwd <- get_network(sim2, 1)
head(get.vertex.attribute.active(nwd, "testatus", at = 1), 25)
head(get.vertex.attribute.active(nwd, "testatus", at = nsteps), 25)


# 5. Analysis ---------------------------------------------------------------

# Compute prevalence for both examples
sim <- mutate_epi(sim, prev = i.num / num)
sim2 <- mutate_epi(sim2, prev = i.num / num)


## --- Plot 1: Prevalence Comparison ---
# Example 1 (constant high transmission) vs. Example 2 (time-varying lower
# transmission) over the same observed network.
par(mfrow = c(1, 2), mar = c(3, 3, 2, 1), mgp = c(2, 1, 0))
plot(sim, y = "i.num",
     main = "Ex 1: Constant Transmission",
     ylab = "Infected Count", xlab = "Time Step",
     mean.col = "firebrick", mean.lwd = 2, mean.smooth = TRUE,
     qnts.col = "firebrick", qnts.alpha = 0.2, qnts.smooth = TRUE,
     legend = FALSE)
plot(sim2, y = "i.num",
     main = "Ex 2: Time-Varying Transmission",
     ylab = "Infected Count", xlab = "Time Step",
     mean.col = "firebrick", mean.lwd = 2, mean.smooth = TRUE,
     qnts.col = "firebrick", qnts.alpha = 0.2, qnts.smooth = TRUE,
     legend = FALSE)


## --- Plot 2: Incidence Comparison ---
par(mfrow = c(1, 2), mar = c(3, 3, 2, 1), mgp = c(2, 1, 0))
plot(sim, y = "si.flow",
     main = "Ex 1: Incidence",
     ylab = "New Infections", xlab = "Time Step",
     mean.col = "steelblue", mean.lwd = 2, mean.smooth = TRUE,
     qnts.col = "steelblue", qnts.alpha = 0.2, qnts.smooth = TRUE,
     legend = FALSE)
plot(sim2, y = "si.flow",
     main = "Ex 2: Incidence",
     ylab = "New Infections", xlab = "Time Step",
     mean.col = "steelblue", mean.lwd = 2, mean.smooth = TRUE,
     qnts.col = "steelblue", qnts.alpha = 0.2, qnts.smooth = TRUE,
     legend = FALSE)


## --- Summary Table ---
df1 <- as.data.frame(sim)
df2 <- as.data.frame(sim2)

data.frame(
  Metric = c("Final prevalence",
             "Cumulative infections",
             "Mean incidence per step"),
  Example1 = c(
    round(mean(df1$prev[df1$time == max(df1$time)], na.rm = TRUE), 3),
    round(mean(tapply(df1$si.flow, df1$sim, sum, na.rm = TRUE))),
    round(mean(df1$si.flow, na.rm = TRUE), 2)
  ),
  Example2 = c(
    round(mean(df2$prev[df2$time == max(df2$time)], na.rm = TRUE), 3),
    round(mean(tapply(df2$si.flow, df2$sim, sum, na.rm = TRUE))),
    round(mean(df2$si.flow, na.rm = TRUE), 2)
  )
)
