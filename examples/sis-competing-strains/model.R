
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


# 1. Network Model Estimation ------------------------------------------------

# Initialize a 1000-node network
nw <- network_initialize(n = 1000)

# Two network formation models, identical in mean degree but differing
# in whether concurrency (simultaneous partnerships) is allowed:
#
#   Model 1 (concurrency allowed): edges-only model with no degree
#     constraints. At 300 edges in a 1000-node network, mean degree = 0.6.
#     Some nodes will naturally have 2+ concurrent ties.
#
#   Model 2 (strict monogamy): adds a concurrent term constrained to 0,
#     prohibiting any node from having more than 1 active tie at a time.
#     Same mean degree (0.6), but no concurrency.
formation.mod1 <- ~edges
target.stats.mod1 <- 300

formation.mod2 <- ~edges + concurrent
target.stats.mod2 <- c(300, 0)

# Partnership duration: mean of 50 weeks (~1 year), same for both models
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 50)
coef.diss

# Fit both network models
est.mod1 <- netest(nw, formation.mod1, target.stats.mod1, coef.diss)
est.mod2 <- netest(nw, formation.mod2, target.stats.mod2, coef.diss)

# Diagnostics
dx.mod1 <- netdx(est.mod1, nsims = nsims, ncores = ncores, nsteps = nsteps,
                 nwstats.formula = ~edges + concurrent + degree(0:5),
                 set.control.ergm = control.simulate.formula(MCMC.burnin = 1e5))
print(dx.mod1)
plot(dx.mod1, plots.joined = FALSE)

dx.mod2 <- netdx(est.mod2, nsims = nsims, ncores = ncores, nsteps = nsteps,
                 nwstats.formula = ~edges + concurrent + degree(0:5),
                 set.control.ergm = control.simulate.formula(MCMC.burnin = 1e5))
print(dx.mod2)
plot(dx.mod2, plots.joined = FALSE)


# 2. Epidemic Model Parameters ------------------------------------------------

# Load custom module functions (strain initialization, infection, recovery)
source("examples/sis-competing-strains/module-fx.R")

# Two competing strains with an infectiousness-duration trade-off:
#
#   Strain 1 ("acute"): highly infectious but short-lived
#     inf.prob = 0.5: 50% per-act transmission probability
#     rec.rate = 0.05: mean infectious duration ~20 weeks
#
#   Strain 2 ("chronic"): low infectiousness but long-lived
#     inf.prob.st2 = 0.01: 1% per-act transmission probability
#     rec.rate.st2 = 0.005: mean infectious duration ~200 weeks
#
#   act.rate = 2: two acts per partnership per week
#   pct.st2 = 0.5: initial infections split 50/50 between strains
#
# The key insight: both strains have very different transmission
# strategies (fast-and-furious vs. slow-and-steady). Which strategy
# "wins" depends critically on the network structure -- specifically,
# whether concurrency is allowed.
param <- param.net(
  inf.prob = 0.5, inf.prob.st2 = 0.01,
  rec.rate = 0.05, rec.rate.st2 = 0.005,
  pct.st2 = 0.5,
  act.rate = 2
)

# Initial conditions: 50 infected nodes (5% prevalence) to seed the epidemic
init <- init.net(i.num = 50)

# Control settings: custom modules for strain initialization, infection,
# and recovery. Module order: initStrain -> infection -> recovery
control <- control.net(
  type = NULL,
  nsims = nsims,
  ncores = ncores,
  nsteps = nsteps,
  initStrain.FUN = init_strain,
  recovery.FUN = recov.2strains,
  infection.FUN = infection.2strains,
  verbose = FALSE
)


# 3. Run Simulations ----------------------------------------------------------

# Model 1: Concurrency allowed (random network)
sim.mod1 <- netsim(est.mod1, param, init, control)
print(sim.mod1)

# Model 2: Strict monogamy (no concurrency)
sim.mod2 <- netsim(est.mod2, param, init, control)
print(sim.mod2)


# 4. Analysis -----------------------------------------------------------------

# Compute total and strain-specific prevalence
sim.mod1 <- mutate_epi(sim.mod1, prev = i.num / num,
                       prev.st1 = i.num.st1 / num,
                       prev.st2 = i.num.st2 / num)
sim.mod2 <- mutate_epi(sim.mod2, prev = i.num / num,
                       prev.st1 = i.num.st1 / num,
                       prev.st2 = i.num.st2 / num)


## --- Plot 1: Strain Competition -- Concurrency vs. Monogamy ---
# The headline result: concurrency reverses which strain dominates.
# With concurrency, the fast-spreading strain 1 dominates.
# Under monogamy, the slow-and-steady strain 2 wins and strain 1 goes extinct.
par(mfrow = c(1, 2), mar = c(3, 3, 2, 1), mgp = c(2, 1, 0))
plot(sim.mod1, y = c("i.num.st1", "i.num.st2"),
     main = "Concurrency Allowed",
     ylab = "Infected (count)", xlab = "Week",
     mean.col = c("firebrick", "steelblue"), mean.lwd = 2, mean.smooth = TRUE,
     qnts.col = c("firebrick", "steelblue"), qnts.alpha = 0.2,
     qnts.smooth = TRUE, legend = TRUE)
plot(sim.mod2, y = c("i.num.st1", "i.num.st2"),
     main = "Strict Monogamy",
     ylab = "Infected (count)", xlab = "Week",
     mean.col = c("firebrick", "steelblue"), mean.lwd = 2, mean.smooth = TRUE,
     qnts.col = c("firebrick", "steelblue"), qnts.alpha = 0.2,
     qnts.smooth = TRUE, legend = TRUE)


## --- Plot 2: Total Prevalence Comparison ---
# How does total disease burden compare across network structures?
par(mfrow = c(1, 1), mar = c(3, 3, 2, 1), mgp = c(2, 1, 0))
plot(sim.mod1, y = "prev",
     main = "Total Prevalence by Network Structure",
     ylab = "Prevalence (I / N)", xlab = "Week",
     mean.col = "firebrick", mean.lwd = 2, mean.smooth = TRUE,
     qnts.col = "firebrick", qnts.alpha = 0.2, qnts.smooth = TRUE,
     legend = FALSE)
plot(sim.mod2, y = "prev",
     mean.col = "steelblue", mean.lwd = 2, mean.smooth = TRUE,
     qnts.col = "steelblue", qnts.alpha = 0.2, qnts.smooth = TRUE,
     legend = FALSE, add = TRUE)
legend("topright",
       legend = c("Concurrency Allowed", "Strict Monogamy"),
       col = c("firebrick", "steelblue"), lwd = 2, bty = "n")


## --- Plot 3: Incidence by Strain ---
# New infections per week, broken out by strain and network model.
# Shows the transmission "velocity" of each strain.
par(mfrow = c(1, 2), mar = c(3, 3, 2, 1), mgp = c(2, 1, 0))
plot(sim.mod1, y = c("si.flow", "si.flow.st2"),
     main = "Incidence: Concurrency",
     ylab = "New Infections / Week", xlab = "Week",
     mean.col = c("firebrick", "steelblue"),
     mean.lwd = 2, mean.smooth = TRUE,
     qnts = FALSE, legend = TRUE)
plot(sim.mod2, y = c("si.flow", "si.flow.st2"),
     main = "Incidence: Monogamy",
     ylab = "New Infections / Week", xlab = "Week",
     mean.col = c("firebrick", "steelblue"),
     mean.lwd = 2, mean.smooth = TRUE,
     qnts = FALSE, legend = TRUE)


## --- Summary Table ---
df1 <- as.data.frame(sim.mod1)
df2 <- as.data.frame(sim.mod2)

# Use last 20% of timesteps to estimate equilibrium
late1 <- df1$time > nsteps * 0.8
late2 <- df2$time > nsteps * 0.8

data.frame(
  Metric = c("Equilibrium prevalence (total)",
             "Strain 1 prevalence",
             "Strain 2 prevalence",
             "Dominant strain",
             "Cumulative infections (St1)",
             "Cumulative infections (St2)"),
  Concurrency = c(
    round(mean(df1$prev[late1], na.rm = TRUE), 3),
    round(mean(df1$prev.st1[late1], na.rm = TRUE), 3),
    round(mean(df1$prev.st2[late1], na.rm = TRUE), 3),
    ifelse(mean(df1$i.num.st1[late1], na.rm = TRUE) >
             mean(df1$i.num.st2[late1], na.rm = TRUE),
           "Strain 1", "Strain 2"),
    round(mean(tapply(df1$si.flow, df1$sim, sum, na.rm = TRUE))),
    round(mean(tapply(df1$si.flow.st2, df1$sim, sum, na.rm = TRUE)))
  ),
  Monogamy = c(
    round(mean(df2$prev[late2], na.rm = TRUE), 3),
    round(mean(df2$prev.st1[late2], na.rm = TRUE), 3),
    round(mean(df2$prev.st2[late2], na.rm = TRUE), 3),
    ifelse(mean(df2$i.num.st1[late2], na.rm = TRUE) >
             mean(df2$i.num.st2[late2], na.rm = TRUE),
           "Strain 1", "Strain 2"),
    round(mean(tapply(df2$si.flow, df2$sim, sum, na.rm = TRUE))),
    round(mean(tapply(df2$si.flow.st2, df2$sim, sum, na.rm = TRUE)))
  )
)


# 5. Sensitivity Analysis: Concurrency Gradient -------------------------------

# At what level of concurrency does the crossover occur? We sweep from
# strict monogamy (concurrent = 0) to the natural concurrency level
# (~120 nodes), tracking equilibrium strain prevalence at each level.

if (interactive()) {

  # Check natural concurrency level in Model 1
  dx.mod1a <- netdx(est.mod1, nsims = 10, ncores = ncores, nsteps = 100,
                    set.control.ergm = control.simulate.formula(MCMC.burnin = 1e5),
                    nwstats.formula = ~edges + concurrent, verbose = FALSE)
  dx.mod1a

  # Sweep concurrency from 0 to 120 in steps of 10
  formation.mod3 <- ~edges + concurrent
  concties.mod3 <- seq(0, 120, 10)
  est.mod3 <- list()
  sim.mod3 <- list()

  for (i in seq_along(concties.mod3)) {
    cat("\rConcurrency sweep:", i, "/", length(concties.mod3),
        "(concurrent ties =", concties.mod3[i], ")    ")
    target.stats.mod3 <- c(300, concties.mod3[i])
    est.mod3[[i]] <- suppressMessages(
      netest(nw, formation.mod3, target.stats.mod3, coef.diss)
    )
    sim.mod3[[i]] <- netsim(est.mod3[[i]], param, init, control)
  }
  cat("\n")

  # Extract equilibrium strain prevalence (last timestep mean)
  i.num.st1.final <- sapply(seq_along(concties.mod3), function(x) {
    tail(as.data.frame(sim.mod3[[x]], out = "mean")[["i.num.st1"]], 1)
  })
  i.num.st2.final <- sapply(seq_along(concties.mod3), function(x) {
    tail(as.data.frame(sim.mod3[[x]], out = "mean")[["i.num.st2"]], 1)
  })

  ## --- Plot 4: Concurrency Gradient ---
  # The dose-response of concurrency on strain competition. Strain 1
  # prevalence rises monotonically with concurrency. Strain 2 initially
  # increases then declines as strain 1 crowds it out. The crossover
  # occurs at approximately 70 concurrent nodes.
  par(mfrow = c(1, 1), mar = c(3, 3, 2, 1), mgp = c(2, 1, 0))
  plot(concties.mod3, i.num.st1.final, type = "b", pch = 19,
       col = "firebrick", lwd = 2,
       xlab = "Nodes with Concurrent Ties",
       ylab = "Equilibrium Infected (count)",
       main = "Strain Competition vs. Concurrency",
       ylim = c(0, max(c(i.num.st1.final, i.num.st2.final), na.rm = TRUE)))
  points(concties.mod3, i.num.st2.final, type = "b", pch = 19,
         col = "steelblue", lwd = 2)
  legend("topleft",
         legend = c("Strain 1 (acute)", "Strain 2 (chronic)"),
         col = c("firebrick", "steelblue"), lwd = 2, pch = 19, bty = "n")

  # Crossover summary
  data.frame(
    concurrent = concties.mod3,
    strain1 = round(i.num.st1.final),
    strain2 = round(i.num.st2.final)
  )
}
