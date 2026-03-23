
##
## Social Diffusion Model
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##
## Authors: Samuel M. Jenness
## Date: September 2018
##

## Load EpiModel
suppressMessages(library(EpiModel))

# Standard Gallery unit test lines
rm(list = ls())
eval(parse(text = print(commandArgs(TRUE)[1])))

if (interactive()) {
  nsims <- 5
  ncores <- 5
  nsteps <- 300
} else {
  nsims <- 2
  ncores <- 2
  nsteps <- 200
}


# Overview -----------------------------------------------------------------
#
# This example demonstrates how EpiModel's SI framework can be repurposed
# to model SOCIAL DIFFUSION -- the spread of ideas, behaviors, or technologies
# through a social network. Unlike infectious disease, social diffusion often
# exhibits "complex contagion": adoption requires social reinforcement from
# multiple contacts, not just exposure to a single carrier.
#
# We compare three diffusion mechanisms on the same network:
#
#   1. Simple contagion (baseline): constant per-contact adoption probability,
#      using EpiModel's built-in SI model. This is equivalent to infectious
#      disease transmission and serves as the reference case.
#
#   2. Threshold diffusion: adoption occurs ONLY when a non-adopter has at
#      least `min.degree` contacts who have already adopted. Below this
#      threshold, adoption probability is exactly 0. This models phenomena
#      like protest participation (need a critical mass of committed peers).
#
#   3. Dose-response diffusion: adoption probability is a smooth logistic
#      function of the number of adopter contacts. More exposure -> higher
#      probability, but there is no hard cutoff. This is an intermediate
#      model between simple and threshold contagion.
#
# The key insight: the SAME network and SAME initial conditions produce
# dramatically different diffusion dynamics depending on the mechanism.
# Complex contagion (scenarios 2-3) is slower to start but can produce
# sudden "tipping point" cascades once enough of the network is seeded.


# 1. Network Model Estimation -----------------------------------------------

# Initialize a 500-node network representing a social contact network.
# This could represent friendships, coworker relationships, or any social
# ties through which information or behavior can spread.
n <- 500
nw <- network_initialize(n)

# Formation model: edges + isolates.
# The isolates term controls how many nodes have degree 0 (no connections),
# shaping the degree distribution. With 600 edges and 20 isolates:
#   - Mean degree = 2 * 600 / 500 = 2.4
#   - 4% of nodes are isolated (unreachable by network diffusion)
#   - The remaining 96% form the "diffusion-accessible" population
#
# A mean degree of 2.4 is important for the threshold model: with lower
# connectivity (e.g., mean degree 1.8), requiring 2 adopter contacts makes
# diffusion nearly impossible because most nodes don't have enough contacts.
# With mean degree 2.4, the threshold model produces a clear delayed-onset
# S-curve rather than stalling entirely.
formation <- ~edges + isolates

target.stats <- c(600, 20)

# Partnership duration of 50 time steps. In a social network context, this
# represents the average duration of an active social tie (e.g., a friendship
# or regular social contact). Longer durations create a more stable network
# where people maintain the same connections over time. This is realistic for
# social relationships and avoids bias in the network diagnostics that arises
# from very short partnership durations (the "edges dissolution approximation"
# used by netest becomes less accurate when durations are short relative to
# the network size).
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 50)
coef.diss

# Fit the ERGM
est <- netest(nw, formation, target.stats, coef.diss)

# Diagnostics: verify that the simulated network maintains target statistics
# and inspect the degree distribution. The degree distribution is important
# because complex contagion depends on having MULTIPLE adopter contacts,
# which requires sufficient connectivity (nodes with degree >= 2).
dx <- netdx(est, nsims = nsims, ncores = ncores, nsteps = nsteps,
            nwstats.formula = ~edges + isolates + degree(0:5))
print(dx)
plot(dx)


# 2. Load Custom Diffusion Modules ------------------------------------------

source("2018-09-SocialDiffusion/module-fx.R")

# Shared initial conditions: 50 initial adopters (10% seed prevalence).
# A meaningful seed is needed because the threshold diffusion model requires
# local clusters of adopters to form before diffusion can begin. With too
# few initial adopters (e.g., 2% seed) on a stable network (duration = 50),
# the threshold model can stall entirely because non-adopters rarely
# encounter enough adopter contacts simultaneously.
init <- init.net(i.num = 50)


# 3. Scenario 1: Simple Contagion (Baseline) --------------------------------
#
# Uses EpiModel's built-in SI model. Each contact between an adopter and
# a non-adopter has a fixed probability of causing adoption, regardless of
# how many OTHER adopter contacts the non-adopter has. This is the standard
# "simple contagion" assumption used in most epidemic models.
#
# With inf.prob = 0.1 and act.rate = 1, the per-partnership per-timestep
# adoption probability is 0.1. A non-adopter with 1 adopter contact has
# the same per-contact probability as one with 5 adopter contacts.

param_simple <- param.net(inf.prob = 0.1, act.rate = 1)

control_simple <- control.net(
  type = "SI",
  nsteps = nsteps,
  nsims = nsims,
  ncores = ncores,
  verbose = FALSE
)

sim_simple <- netsim(est, param_simple, init, control_simple)
print(sim_simple)


# 4. Scenario 2: Threshold Diffusion ----------------------------------------
#
# A non-adopter only considers adoption when they have at least `min.degree`
# contacts who are already adopters. Below this threshold, adoption probability
# is exactly 0 -- no matter how persuasive any single contact is.
#
# This models "complex contagion" with a hard threshold. Examples:
#   - Technology adoption: "I'll only switch to a new app if at least 2 of
#     my friends already use it"
#   - Behavior change: "I need to see multiple role models before I change"
#   - Collective action: "I'll join the protest only if enough peers commit"
#
# Parameters:
#   inf.prob = 0.5:   per-act adoption probability ONCE the threshold is met.
#                     This is deliberately higher than simple contagion's 0.1
#                     to demonstrate a key insight: even with 5x higher per-act
#                     probability, threshold diffusion is STILL much slower
#                     because the threshold requirement blocks most adoption
#                     attempts entirely.
#   act.rate = 1:     acts per partnership per timestep (same as simple)
#   min.degree = 2:   minimum adopter contacts required before adoption is
#                     possible. With min.degree = 2, a non-adopter must have
#                     at least 2 current partners who have adopted.
#
# Expected dynamics: dramatically slower initial spread than simple contagion.
# With a stable network (duration = 50), non-adopters keep the same partners
# for extended periods. If those partners haven't adopted, the non-adopter
# is locked out of the diffusion process until network turnover brings in
# new (adopter) contacts. Once enough local clusters form, diffusion
# accelerates rapidly through the rest of the network.

param_threshold <- param.net(inf.prob = 0.5, act.rate = 1, min.degree = 2)

control_threshold <- control.net(
  type = NULL,
  nsteps = nsteps,
  nsims = nsims,
  ncores = ncores,
  infection.FUN = diffuse_threshold,
  verbose = FALSE
)

sim_threshold <- netsim(est, param_threshold, init, control_threshold)
print(sim_threshold)


# 5. Scenario 3: Dose-Response Diffusion -------------------------------------
#
# Adoption probability is a smooth logistic function of the number of
# adopter contacts ("exposure"):
#
#   P(adopt per act) = plogis(beta0 + beta1 * exposure)
#                    = 1 / (1 + exp(-(beta0 + beta1 * exposure)))
#
# This is a "soft threshold": no hard cutoff, but adoption probability
# increases smoothly with more adopter contacts. It's intermediate between
# simple contagion (flat probability) and threshold contagion (step function).
#
# Parameters:
#   beta0 = -5.0:  intercept. At exposure = 0 (no adopter contacts), the
#                  adoption probability per act is plogis(-5) = 0.007.
#                  This is very low but nonzero (unlike the threshold model).
#   beta1 = 1.5:   slope. Each additional adopter contact increases the
#                  log-odds by 1.5. The adoption probabilities by exposure:
#                    exposure = 0: plogis(-5.0)       = 0.007
#                    exposure = 1: plogis(-5.0 + 1.5) = 0.029
#                    exposure = 2: plogis(-5.0 + 3.0) = 0.119
#                    exposure = 3: plogis(-5.0 + 4.5) = 0.378
#                    exposure = 4: plogis(-5.0 + 6.0) = 0.731
#   act.rate = 1:  acts per partnership per timestep (same as other scenarios)
#
# Expected dynamics: intermediate between simple and threshold. With a single
# adopter contact, adoption probability is very low (0.029), so diffusion
# requires SOME social reinforcement to proceed quickly. But unlike the
# threshold model, there is no hard cutoff -- even isolated exposure can
# (rarely) cause adoption, preventing the complete stalling that occurs
# with threshold diffusion when clusters are sparse.

param_dose <- param.net(beta0 = -5.0, beta1 = 1.5, act.rate = 1)

control_dose <- control.net(
  type = NULL,
  nsteps = nsteps,
  nsims = nsims,
  ncores = ncores,
  infection.FUN = diffuse_dose_response,
  verbose = FALSE
)

sim_dose <- netsim(est, param_dose, init, control_dose)
print(sim_dose)


# 6. Comparative Analysis ---------------------------------------------------

# Compute derived measures for all three scenarios
sim_simple <- mutate_epi(sim_simple, prev = i.num / num)
sim_threshold <- mutate_epi(sim_threshold, prev = i.num / num)
sim_dose <- mutate_epi(sim_dose, prev = i.num / num)


## --- Plot 1: Adoption Prevalence Comparison (3 Scenarios) ---
# This is the central comparison. All three scenarios use the same network,
# the same initial seed, and the same number of acts per partnership.
# The only difference is the MECHANISM by which exposure translates into
# adoption probability.
#
# Expected pattern:
#   - Simple contagion (green): fastest initial spread, smooth S-curve
#   - Dose-response (blue): intermediate speed, slightly delayed S-curve
#   - Threshold (red): slowest start, with a dramatic delayed onset followed
#     by rapid acceleration once critical mass is reached

par(mfrow = c(1, 1), mar = c(3, 3, 2, 1), mgp = c(2, 1, 0))

plot(sim_simple, y = "prev",
     main = "Adoption Prevalence: Three Diffusion Mechanisms",
     ylab = "Prevalence (Fraction Adopted)", xlab = "Time Steps",
     mean.col = "forestgreen", mean.lwd = 2, mean.smooth = TRUE,
     qnts.col = "forestgreen", qnts.alpha = 0.2, qnts.smooth = TRUE,
     legend = FALSE, ylim = c(0, 1))
plot(sim_threshold, y = "prev", add = TRUE,
     mean.col = "firebrick", mean.lwd = 2, mean.smooth = TRUE,
     qnts.col = "firebrick", qnts.alpha = 0.2, qnts.smooth = TRUE,
     legend = FALSE)
plot(sim_dose, y = "prev", add = TRUE,
     mean.col = "steelblue", mean.lwd = 2, mean.smooth = TRUE,
     qnts.col = "steelblue", qnts.alpha = 0.2, qnts.smooth = TRUE,
     legend = FALSE)
legend("bottomright",
       legend = c("Simple Contagion", "Threshold (min = 2)", "Dose-Response"),
       col = c("forestgreen", "firebrick", "steelblue"),
       lwd = 2, bty = "n")


## --- Plot 2: Adoption Incidence Comparison (3 Scenarios) ---
# New adoptions per time step. Simple contagion shows a single smooth peak.
# Threshold diffusion may show a delayed peak or a more abrupt onset.
# Dose-response should be intermediate.

plot(sim_simple, y = "si.flow",
     main = "New Adoptions per Time Step",
     ylab = "New Adoptions", xlab = "Time Steps",
     mean.col = "forestgreen", mean.lwd = 2, mean.smooth = TRUE,
     qnts.col = "forestgreen", qnts.alpha = 0.2, qnts.smooth = TRUE,
     legend = FALSE)
plot(sim_threshold, y = "si.flow", add = TRUE,
     mean.col = "firebrick", mean.lwd = 2, mean.smooth = TRUE,
     qnts.col = "firebrick", qnts.alpha = 0.2, qnts.smooth = TRUE,
     legend = FALSE)
plot(sim_dose, y = "si.flow", add = TRUE,
     mean.col = "steelblue", mean.lwd = 2, mean.smooth = TRUE,
     qnts.col = "steelblue", qnts.alpha = 0.2, qnts.smooth = TRUE,
     legend = FALSE)
legend("topright",
       legend = c("Simple Contagion", "Threshold (min = 2)", "Dose-Response"),
       col = c("forestgreen", "firebrick", "steelblue"),
       lwd = 2, bty = "n")


## --- Plot 3: Individual Scenario Detail (3-Panel) ---
# Each panel shows the S (non-adopter) and I (adopter) compartment counts
# for one scenario, with quantile bands across simulations. This lets us
# see the stochastic variability within each scenario.

par(mfrow = c(1, 3), mar = c(3, 3, 2, 1), mgp = c(2, 1, 0))

plot(sim_simple, y = c("s.num", "i.num"),
     main = "Simple Contagion",
     ylab = "Count", xlab = "Time Steps",
     mean.col = c("steelblue", "firebrick"),
     mean.lwd = 2, mean.smooth = TRUE,
     qnts.col = c("steelblue", "firebrick"),
     qnts.alpha = 0.2, qnts.smooth = TRUE,
     legend = TRUE)

plot(sim_threshold, y = c("s.num", "i.num"),
     main = "Threshold (min = 2)",
     ylab = "Count", xlab = "Time Steps",
     mean.col = c("steelblue", "firebrick"),
     mean.lwd = 2, mean.smooth = TRUE,
     qnts.col = c("steelblue", "firebrick"),
     qnts.alpha = 0.2, qnts.smooth = TRUE,
     legend = TRUE)

plot(sim_dose, y = c("s.num", "i.num"),
     main = "Dose-Response",
     ylab = "Count", xlab = "Time Steps",
     mean.col = c("steelblue", "firebrick"),
     mean.lwd = 2, mean.smooth = TRUE,
     qnts.col = c("steelblue", "firebrick"),
     qnts.alpha = 0.2, qnts.smooth = TRUE,
     legend = TRUE)


## --- Summary Table ---

df_simple <- as.data.frame(sim_simple)
df_threshold <- as.data.frame(sim_threshold)
df_dose <- as.data.frame(sim_dose)
last_t <- max(df_simple$time)

# Time to 50% adoption: the first time step where mean prevalence >= 0.5.
# For scenarios that never reach 50% within the simulation, this returns NA.
calc_t50 <- function(df) {
  mean_prev <- tapply(df$prev, df$time, mean, na.rm = TRUE)
  reached <- which(mean_prev >= 0.5)
  if (length(reached) == 0) NA else min(reached)
}
t50_simple <- calc_t50(df_simple)
t50_threshold <- calc_t50(df_threshold)
t50_dose <- calc_t50(df_dose)

# Summary statistics
data.frame(
  Metric = c("Final prevalence",
             "Time to 50% adoption",
             "Peak new adoptions/step",
             "Cumulative adoptions",
             "Mean prevalence"),
  Simple = c(round(mean(df_simple$prev[df_simple$time == last_t],
                        na.rm = TRUE), 3),
             t50_simple,
             round(max(tapply(df_simple$si.flow, df_simple$time,
                              mean, na.rm = TRUE), na.rm = TRUE), 1),
             round(sum(tapply(df_simple$si.flow, df_simple$time,
                              mean, na.rm = TRUE), na.rm = TRUE)),
             round(mean(df_simple$prev, na.rm = TRUE), 3)),
  Threshold = c(round(mean(df_threshold$prev[df_threshold$time == last_t],
                           na.rm = TRUE), 3),
                t50_threshold,
                round(max(tapply(df_threshold$si.flow, df_threshold$time,
                                 mean, na.rm = TRUE), na.rm = TRUE), 1),
                round(sum(tapply(df_threshold$si.flow, df_threshold$time,
                                 mean, na.rm = TRUE), na.rm = TRUE)),
                round(mean(df_threshold$prev, na.rm = TRUE), 3)),
  Dose_Response = c(round(mean(df_dose$prev[df_dose$time == last_t],
                               na.rm = TRUE), 3),
                    t50_dose,
                    round(max(tapply(df_dose$si.flow, df_dose$time,
                                     mean, na.rm = TRUE), na.rm = TRUE), 1),
                    round(sum(tapply(df_dose$si.flow, df_dose$time,
                                     mean, na.rm = TRUE), na.rm = TRUE)),
                    round(mean(df_dose$prev, na.rm = TRUE), 3))
)

# Key takeaway: all three scenarios eventually reach full adoption (the
# network is well-connected enough), but the SPEED and SHAPE of diffusion
# differ dramatically. Simple contagion produces the classic smooth S-curve.
# Threshold contagion delays diffusion until critical mass is reached, then
# accelerates sharply. Dose-response is a smooth middle ground.
#
# Note that the threshold model has 5x HIGHER per-act adoption probability
# (0.5 vs. 0.1) yet is still ~4x SLOWER to reach 50% adoption. This
# powerfully demonstrates why complex contagion produces fundamentally
# different dynamics: the MECHANISM of social influence (requiring multiple
# reinforcing contacts) matters more than the probability per contact.
#
# This has implications for intervention design: seeding adoption in highly
# connected clusters is much more important for complex contagion than for
# simple contagion, where even random seeding produces rapid diffusion.
