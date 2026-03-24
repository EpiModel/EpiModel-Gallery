
##
## Epidemics with Multiple (Multilayer) Networks
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##
## Authors: Chad Klumb, Samuel M. Jenness
## Date: December 2022
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
  nsims <- 2
  ncores <- 2
  nsteps <- 100
}


# 1. Overview ----------------------------------------------------------------

# This example demonstrates how to model an epidemic on a MULTILAYER NETWORK
# in EpiModel. A multilayer network has two (or more) layers that share the
# same node set but have different edge sets -- representing distinct types of
# relationships (e.g., main partnerships vs. casual contacts, or sexual
# contacts vs. needle-sharing contacts).
#
# Key features demonstrated:
#   - Two network layers with different formation and dissolution models
#   - Cross-layer dependency: degree in one layer affects edge formation in the
#     other (via edge-dependent nodal attributes and nodefactor ERGM terms)
#   - The dat.updates callback for maintaining cross-layer attributes during
#     the simulation
#   - The multilayer() helper for specifying per-layer network statistics
#
# We use a simple SI model in a closed population to isolate the multilayer
# network mechanics from disease model complexity. No custom module-fx.R is
# needed -- this example uses EpiModel's built-in SI modules.


# 2. Network Setup -----------------------------------------------------------

# Initialize a network of 500 nodes with a binary demographic attribute
# ("race"). This attribute will be used for assortative mixing (homophily)
# in layer 1's formation model.
n <- 500
nw <- network_initialize(n)
nw <- set_vertex_attribute(nw, "race", rep(0:1, length.out = n))


# --- Generating starting degree values with san() ---
#
# Because the two layers depend on each other's degree distribution (via
# cross-layer nodefactor terms), we need reasonable starting values for
# the degree attributes BEFORE fitting the ERGM models. We use san()
# (simulated annealing) to generate edge sets consistent with our target
# statistics.

# Step 1: Generate layer 1 edges consistent with target mean degree and
# race homophily. Layer 1 has target stats of 90 total edges and 60 that
# are race-homophilous (same-race partnerships).
nw_san <- san(nw ~ edges + nodematch("race"), target.stats = c(90, 60))

# Record each node's degree in layer 1 as a vertex attribute, and create
# a binary indicator (0 vs. 1+) for the cross-layer term.
nw_san <- set_vertex_attribute(nw_san, "deg.net1", get_degree(nw_san))
nw_san <- set_vertex_attribute(nw_san, "deg1+.net1",
                               pmin(get_degree(nw_san), 1))

# Step 2: Generate layer 2 edges, accounting for the cross-layer constraint.
# The nodefactor("deg1+.net1") term ensures that nodes already active in
# layer 1 are less likely to form layer 2 edges (negative correlation).
# The Sum() term provides additional control over the degree distribution
# by linking layer 2 degree to layer 1 degree categories.
nw_san_2 <- san(nw_san ~ edges + degree(1) + nodefactor("deg1+.net1")
                  + Sum(matrix(seq_len(network.size(nw)), nrow = 1) ~
                        degrange(1, by = "deg.net1", levels = -1),
                    label = "Sum"),
                target.stats = c(75, 120, 10, 10))
nw_san <- set_vertex_attribute(nw_san, "deg1+.net2",
                               pmin(get_degree(nw_san_2), 1))

# Verify that both SAN networks approximate their target statistics.
# These don't need to be exact -- they are starting conditions for ERGM.
# SAN verification (Layer 1 targets: 90 edges, 60 homophilous)
summary(nw_san ~ edges + nodematch("race") + nodefactor("deg1+.net2"))

# SAN verification (Layer 2 targets: 75 edges, 120 degree-1)
summary(nw_san_2 ~ edges + degree(1) + nodefactor("deg1+.net1"))

# Set the binary degree attributes on the original network object for netest.
nw <- set_vertex_attribute(nw, "deg1+.net1", pmin(get_degree(nw_san), 1))
nw <- set_vertex_attribute(nw, "deg1+.net2", pmin(get_degree(nw_san_2), 1))


# 3. Network Model Estimation -----------------------------------------------

# Formation models:
#   Layer 1: edges + race homophily + cross-layer degree constraint
#   Layer 2: edges + degree(1) distribution + cross-layer degree constraint
#
# The cross-layer nodefactor terms create NEGATIVE degree correlation across
# layers: having edges in one layer makes edges in the other less likely.
# This models a realistic constraint -- individuals have finite relational
# capacity, so being active in one partnership type reduces availability for
# the other.
formation.1 <- ~edges + nodematch("race") + nodefactor("deg1+.net2")
formation.2 <- ~edges + degree(1) + nodefactor("deg1+.net1")

# Target statistics:
#   Layer 1: 90 edges total, 60 race-homophilous, 10 among those with layer-2
#            edges (out of 90 total -- strong negative correlation)
#   Layer 2: 75 edges total, 120 degree-1 nodes, 10 among those with layer-1
#            edges (out of 75 total -- strong negative correlation)
target.stats.1 <- c(90, 60, 10)
target.stats.2 <- c(75, 120, 10)

# Dissolution models:
#   Layer 1: mean partnership duration of 100 time steps (longer, more stable)
#   Layer 2: mean partnership duration of 75 time steps (shorter, more transient)
# Different durations demonstrate that each layer can have its own dissolution
# dynamics independently.
coef.diss.1 <- dissolution_coefs(dissolution = ~offset(edges), duration = 100)
coef.diss.2 <- dissolution_coefs(dissolution = ~offset(edges), duration = 75)

# Fit the two ERGM models. Although we think of this as ONE multilayer
# network, estimation is done separately for each layer. The layers are
# linked through the shared nodal attributes (deg1+.net1, deg1+.net2).
est.1 <- netest(nw, formation.1, target.stats.1, coef.diss.1)
est.2 <- netest(nw, formation.2, target.stats.2, coef.diss.2)


# 4. Model Diagnostics ------------------------------------------------------

# Run diagnostics for each layer. Because the network is relatively small
# (500 nodes) and some target statistics are small (10 for cross-layer terms),
# there may be more variability than in a single-layer model.

# --- Layer 1 diagnostics ---
dx.1 <- netdx(est.1, nsims = nsims, ncores = ncores, nsteps = nsteps)
print(dx.1)
plot(dx.1)

# --- Layer 2 diagnostics ---
dx.2 <- netdx(est.2, nsims = nsims, ncores = ncores, nsteps = nsteps)
print(dx.2)
plot(dx.2)


# 5. Epidemic Model Setup ----------------------------------------------------

# Simple SI model parameters. High transmission probability and act rate are
# deliberate -- they generate a clearly visible epidemic so we can focus on
# the network mechanics rather than fine-tuning epidemiological realism.
#
# Parameters:
#   inf.prob: per-act transmission probability (0.5 = 50% per act)
#   act.rate: number of acts per partnership per time step
param <- param.net(inf.prob = 0.5, act.rate = 2)

# Initial conditions: 10 infected individuals (2% prevalence in 500 nodes)
init <- init.net(i.num = 10)


# --- Cross-layer update function ---
#
# This is the KEY MECHANISM for multilayer feedback. During each time step,
# EpiModel resimulates each network layer sequentially. Before resimulating
# layer k, it calls this function to update the cross-layer degree attributes
# that layer k's formation model depends on.
#
# The `network` argument tells you which layer is about to be resimulated:
#   0 = before layer 1: look up current degree in layer 2, update deg1+.net2
#   1 = before layer 2: look up current degree in layer 1, update deg1+.net1
#   2 = after layer 2:  recalculate deg1+.net2 for network summary statistics
#
# Without this function, the cross-layer degree attributes would be frozen at
# their initial values and the layers would evolve independently.
network_layer_updates <- function(dat, at, network) {
  if (network == 0) {
    dat <- set_attr(dat, "deg1+.net2", pmin(get_degree(dat$run$el[[2]]), 1))
  } else if (network == 1) {
    dat <- set_attr(dat, "deg1+.net1", pmin(get_degree(dat$run$el[[1]]), 1))
  } else if (network == 2) {
    dat <- set_attr(dat, "deg1+.net2", pmin(get_degree(dat$run$el[[2]]), 1))
  }
  return(dat)
}

# Control settings:
#   type = "SI": use EpiModel's built-in SI modules
#   tergmLite = TRUE: use the lightweight edgelist representation (required
#     for the dat$run$el[[]] access in network_layer_updates)
#   resimulate.network = TRUE: redraw both layers each time step
#   dat.updates: our cross-layer degree update callback
#   nwstats.formula: uses multilayer() to track per-layer network statistics.
#     Layer 1 uses the default "formation" formula; layer 2 additionally
#     tracks the full degree distribution (degree 0 through 4).
control <- control.net(
  type = "SI",
  nsteps = nsteps,
  nsims = nsims,
  ncores = ncores,
  tergmLite = TRUE,
  resimulate.network = TRUE,
  dat.updates = network_layer_updates,
  nwstats.formula = multilayer(
    "formation",
    ~edges + nodefactor("deg1+.net1") + degree(0:4)
  )
)

# Run the simulation. Passing a LIST of netest objects tells EpiModel this
# is a multilayer model. The list order determines the layer numbering.
sim <- netsim(list(est.1, est.2), param, init, control)


# 6. Analysis ---------------------------------------------------------------

## --- Network Formation Diagnostics During Simulation ---
# Verify that network structure is maintained during the epidemic. These
# plots should show the simulated formation statistics tracking close to
# their targets throughout the run.

par(mfrow = c(1, 1), mar = c(3, 3, 2, 1), mgp = c(2, 1, 0))

# --- Layer 1 formation stats during simulation ---
print(sim, network = 1)
plot(sim, network = 1, type = "formation",
     main = "Layer 1: Formation Diagnostics During Simulation")

# --- Layer 2 formation stats during simulation ---
print(sim, network = 2)
plot(sim, network = 2, type = "formation",
     main = "Layer 2: Formation Diagnostics During Simulation")


## --- Plot 1: SI Compartment Counts ---
plot(sim, y = c("s.num", "i.num"),
     main = "SI Epidemic on Multilayer Network",
     ylab = "Count", xlab = "Time Steps",
     mean.col = c("steelblue", "firebrick"),
     mean.lwd = 2, mean.smooth = TRUE,
     qnts.col = c("steelblue", "firebrick"),
     qnts.alpha = 0.2, qnts.smooth = TRUE,
     legend = TRUE)


## --- Plot 2: Prevalence ---
sim <- mutate_epi(sim, prev = i.num / num)
plot(sim, y = "prev",
     main = "Prevalence over Time",
     ylab = "Prevalence", xlab = "Time Steps",
     mean.col = "firebrick", mean.lwd = 2, mean.smooth = TRUE,
     qnts.col = "firebrick", qnts.alpha = 0.2, qnts.smooth = TRUE)


## --- Plot 3: Incidence ---
plot(sim, y = "si.flow",
     main = "New Infections per Time Step",
     ylab = "Incidence", xlab = "Time Steps",
     mean.col = "firebrick", mean.lwd = 2, mean.smooth = TRUE,
     qnts.col = "firebrick", qnts.alpha = 0.2, qnts.smooth = TRUE)


## --- Summary Table ---
df <- as.data.frame(sim)
last_t <- max(df$time)

# Summary statistics
data.frame(
  Metric = c("Population size", "Final susceptible", "Final infected",
             "Final prevalence", "Cumulative new infections",
             "Peak incidence (per time step)"),
  Value = c(n,
            mean(df$s.num[df$time == last_t]),
            mean(df$i.num[df$time == last_t]),
            round(mean(df$prev[df$time == last_t]), 3),
            sum(df$si.flow, na.rm = TRUE),
            max(df$si.flow, na.rm = TRUE))
)
