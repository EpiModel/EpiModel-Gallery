
##
## Multilayer Networks: Independent Layers (the simplest multilayer model)
## EpiModel Gallery (https://github.com/EpiModel/EpiModel-Gallery)
##
## Authors: Samuel M. Jenness
## Date: June 2026
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

# This example is the entry point to MULTILAYER NETWORK modeling in EpiModel.
# A multilayer network has two (or more) layers that share the same node set
# but have different edge sets, representing distinct kinds of relationships
# (e.g., main vs. casual partnerships). An epidemic can spread over edges in
# EITHER layer.
#
# The layers here are INDEPENDENT: a person's number of partners in one layer
# does not constrain the other. That keeps the mechanics minimal. The only
# multilayer-specific step is passing netsim() a LIST of fitted network models.
# There is no simulated annealing, no cross-layer nodal attribute, and no
# dat.updates callback. Those are needed only when the layers depend on each
# other, which is the companion cross-layer dependency example (multinets).
#
# If you have built a single-layer edges model with netest() and netsim(),
# you already know almost everything here.
#
# We use a simple SI model in a closed population so the focus stays on the
# network mechanics. No custom module-fx.R is needed (built-in SI modules).


# 2. Network Setup -----------------------------------------------------------

# One shared node set carries both layers. A binary "race" attribute is used
# for assortative mixing (homophily) in layer 1.
n <- 500
nw <- network_initialize(n)
nw <- set_vertex_attribute(nw, "race", rep(0:1, length.out = n))


# 3. Network Model Estimation ------------------------------------------------

# Each layer is an ordinary single-layer TERGM, fit exactly as in the earlier
# Gallery examples. Nothing here is multilayer-specific yet.
#
#   Layer 1 (main):   steady, long-duration ties with race homophily.
#   Layer 2 (casual): higher-turnover, short-duration ties, degree-1 driven.
#
# Layer 2's formation references ONLY its own structure (degree(1)); it does
# not reference layer 1. That independence is what removes the need for any
# cross-layer machinery.
est1 <- netest(
  nw,
  formation = ~edges + nodematch("race"),
  target.stats = c(90, 60),                # 90 edges total, 60 same-race
  coef.diss = dissolution_coefs(~offset(edges), duration = 200) # long-lasting
)

est2 <- netest(
  nw,
  formation = ~edges + degree(1),
  target.stats = c(75, 120),               # 75 edges total, 120 degree-1 nodes
  coef.diss = dissolution_coefs(~offset(edges), duration = 20)  # transient
)


# 4. Model Diagnostics ------------------------------------------------------

# Diagnose each layer separately, as you would any single-layer model.
dx1 <- netdx(est1, nsims = nsims, ncores = ncores, nsteps = nsteps)
print(dx1)
plot(dx1)

dx2 <- netdx(est2, nsims = nsims, ncores = ncores, nsteps = nsteps)
print(dx2)
plot(dx2)


# 5. Epidemic Model Setup ----------------------------------------------------

# Simple SI parameters. A high transmission probability and act rate are
# deliberate, generating a clearly visible epidemic so the focus stays on the
# network mechanics rather than epidemiological realism.
#   inf.prob: per-act transmission probability
#   act.rate: acts per partnership per time step
param <- param.net(inf.prob = 0.5, act.rate = 2)
init <- init.net(i.num = 10)

# The control settings differ from a single-layer model in only two ways, and
# one of them is optional:
#   resimulate.network = TRUE: redraw both dynamic layers each time step.
#   nwstats.formula = multilayer(...): OPTIONAL per-layer diagnostics, mapping
#     one formula to each layer by position. Omit it and EpiModel tracks each
#     layer's own formation formula by default.
# Note what is ABSENT relative to the cross-layer example: no dat.updates
# callback, because neither layer's formation reads from the other.
control <- control.net(
  type = "SI",
  nsteps = nsteps,
  nsims = nsims,
  ncores = ncores,
  tergmLite = TRUE,
  resimulate.network = TRUE,
  nwstats.formula = multilayer(
    ~edges + nodematch("race"),            # layer 1 diagnostics
    ~edges + degree(0:3)                   # layer 2 diagnostics
  )
)

# Passing a LIST of netest objects is the one step that makes this multilayer.
# The list order sets the layer numbering (est1 = layer 1, est2 = layer 2).
sim <- netsim(list(est1, est2), param, init, control)


# 6. Analysis ---------------------------------------------------------------

## --- Per-layer formation diagnostics during the simulation ---
par(mfrow = c(1, 1), mar = c(3, 3, 2, 1), mgp = c(2, 1, 0))

print(sim, network = 1)
plot(sim, network = 1, type = "formation", main = "Layer 1 (main) formation")

print(sim, network = 2)
plot(sim, network = 2, type = "formation", main = "Layer 2 (casual) formation")


## --- Plot 1: SI Compartment Counts ---
plot(sim, y = c("s.num", "i.num"),
     main = "SI Epidemic on Two Independent Layers",
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

# The short-duration casual layer tends to drive early spread, while the
# long-duration main layer sustains transmission over time. Because
# transmission can occur on either layer, the epidemic reflects both.
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
