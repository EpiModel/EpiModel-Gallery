
##
## Epidemics with Multiple Networks
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##
## Author: Chad Klumb (borrowing some code and comments from the SEIR example)
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

# Simulate data consistent with the target statistics we will use
nw <- network_initialize(500)
nw <- set_vertex_attribute(nw, "race", rep(c(0,1), length.out = 500))
nw_san <- san(nw ~ edges + nodematch("race"), target.stats = c(90, 60))
nw_san <- set_vertex_attribute(nw_san, "deg.net1", get_degree(nw_san))
nw_san <- set_vertex_attribute(nw_san, "deg1+.net1",
                               pmin(get_degree(nw_san), 1))
nw_san_2 <- san(nw_san ~ edges + degree(1) + nodefactor("deg1+.net1")
                  + Sum(matrix(seq_len(network.size(nw)), nrow = 1) ~
                        degrange(1, by = "deg.net1", levels = -1),label="Sum"),
                target.stats = c(75, 120, 10, 10))
nw_san <- set_vertex_attribute(nw_san, "deg1+.net2", pmin(get_degree(nw_san_2), 1))

# confirm we have hit our intended targets for both edge types
summary(nw_san ~ edges + nodematch("race") + nodefactor("deg1+.net2"))
summary(nw_san_2 ~ edges + degree(1) + nodefactor("deg1+.net1"))

# set edge-dependent nodal attributes on the network object for use in netest
nw <- set_vertex_attribute(nw, "deg1+.net1", pmin(get_degree(nw_san), 1))
nw <- set_vertex_attribute(nw, "deg1+.net2", pmin(get_degree(nw_san_2), 1))

# Network model estimation ------------------------------------------------

# Define the formation models
formation.1 <- ~edges + nodematch("race") + nodefactor("deg1+.net2")
formation.2 <- ~edges + degree(1) + nodefactor("deg1+.net1")

# Input the appropriate target statistics for each term
target.stats.1 <- c(90, 60, 10)
target.stats.2 <- c(75, 120, 10)

# Parameterize the dissolution models
coef.diss.1 <- dissolution_coefs(dissolution = ~offset(edges), duration = 100)
coef.diss.2 <- dissolution_coefs(dissolution = ~offset(edges), duration = 75)

# Fit the models
est.1 <- netest(nw, formation.1, target.stats.1, coef.diss.1)
est.2 <- netest(nw, formation.2, target.stats.2, coef.diss.2)

# Model diagnostics
dx.1 <- netdx(est.1, nsims = nsims, ncores = ncores, nsteps = nsteps)
print(dx.1)
plot(dx.1)

dx.2 <- netdx(est.2, nsims = nsims, ncores = ncores, nsteps = nsteps)
print(dx.2)
plot(dx.2)

# Epidemic model simulation -----------------------------------------------

# Model parameters
param <- param.net(inf.prob = 0.5, act.rate = 2,
                   ei.rate = 0.01, ir.rate = 0.01)

# Initial conditions
init <- init.net(i.num = 10)

# Control settings

# This function defines activities that should occur between resimulation of individual
#   network layers
network_layer_updates <- function(dat, at, network) {
  if (network == 0) {
    # update deg1+.net2 prior to network 1 simulation,
    # as network 1's formation model depends on this
    # nodal attribute
    dat <- set_attr(dat, "deg1+.net2", pmin(get_degree(dat$run$el[[2]]), 1))
  } else if (network == 1) {
    # update deg1+.net1 prior to network 2 simulation,
    # as network 2's formation model depends on this
    # nodal attribute
    dat <- set_attr(dat, "deg1+.net1", pmin(get_degree(dat$run$el[[1]]), 1))
  } else if (network == 2) {
    # update deg1+.net2 after network 2 simulation,
    # for summary statistics
    dat <- set_attr(dat, "deg1+.net2", pmin(get_degree(dat$run$el[[2]]), 1))
  }
  return(dat)
}

# Note the nwstats.formula uses the multilayer function to define the stats
#   specific to each layer
control <- control.net(type = "SI",
                       nsteps = nsteps,
                       nsims = nsims,
                       ncores = ncores,
                       tergmLite = TRUE,
                       resimulate.network = TRUE,
                       dat.updates = network_layer_updates,
                       nwstats.formula = multilayer("formation",
                                                    ~edges + nodefactor("deg1+.net1") + degree(0:4)))

# Run the network model simulation with netsim
sim <- netsim(list(est.1, est.2), param, init, control)

# Examine results
print(sim, network = 1)
print(sim, network = 2)
plot(sim, network = 1, type = "formation")
plot(sim, network = 2, type = "formation")
