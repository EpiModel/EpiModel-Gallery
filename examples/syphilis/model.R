##
## Syphilis Progression Model
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##
## Authors: Samuel M. Jenness, Yuan Zhao
## Date: November 2018
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
  nsims <- 1
  ncores <- 1
  nsteps <- 50
}


# 1. Network Model Estimation ------------------------------------------------

# Initialize a 1000-node network representing a sexual contact network.
n <- 1000
nw <- network_initialize(n)

# Formation model: edges + isolates. The isolates term controls the number
# of nodes with degree 0, shaping the degree distribution beyond what
# edges alone can achieve.
formation <- ~edges + isolates

# Target: 300 edges (mean degree 0.6) with 480 isolates (48% with no current
# partner). This creates a network where about half the population is not
# sexually active at any given time.
target.stats <- c(300, 480)

# Partnership duration of 10 weeks. No departure rate adjustment is needed
# because this model has no vital dynamics (closed population).
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 10)
coef.diss

# Fit the ERGM
est <- netest(nw, formation, target.stats, coef.diss)

# Diagnostics: verify formation targets and degree distribution
dx <- netdx(est, nsims = nsims, ncores = ncores, nsteps = nsteps,
            nwstats.formula = ~edges + isolates + degree(0:3))
print(dx)
plot(dx)

# Note: the diagnostics are roughly ~5% off here due to the use of the "edges
# dissolution approximation" for quickly fitting the TERGM in the context of a
# short overall partner duration (10 weeks). There are various ways to handle
# this (see ?netest) but ignore this issue for current purposes.


# 2. Epidemic Model Setup ----------------------------------------------------

# Load custom module functions (infection, progression, treatment/screening)
source("examples/syphilis/module-fx.R")

# Syphilis natural history parameters:
#
# Stage durations (expressed as 1 / mean weeks in stage):
#   Incubation -> Primary:    ~4 weeks   (ipr.rate; CDC range 10-90 days,
#                                          median ~21 days)
#   Primary -> Secondary:     ~9 weeks   (prse.rate)
#   Secondary -> Early Lat:   ~17 weeks  (seel.rate)
#   Early Lat -> Late Lat:    ~22 weeks  (elll.rate; CDC defines early
#                                          latent as < 1 year since
#                                          acquisition, so ~52 wks total)
#   Late Lat -> Tertiary:     ~29 years  (llter.rate) -- most never progress
#
# Stage-specific transmission probability (per act). Grounded in CDC's
# clinical framing that sexual transmission occurs primarily via contact
# with mucocutaneous syphilitic lesions, and the surveillance definition
# of "infectious syphilis" = primary + secondary + early latent.
# References:
#   - CDC, About Syphilis (https://www.cdc.gov/syphilis/about/index.html)
#   - CDC STI Treatment Guidelines, Syphilis (2021)
#     (https://www.cdc.gov/std/treatment-guidelines/syphilis.htm)
#   - Garnett GP, Aral SO, Hoyle DV, Cates W, Anderson RM (1997). The
#     natural history of syphilis. Sex Transm Dis 24(4):185-200.
#
#   Primary, Secondary:   0.18 per act (inf.prob.early). Mucocutaneous
#                         lesions (chancre, mucous patches, condylomata
#                         lata) are present -- the conventional route of
#                         sexual transmission.
#   Early latent:         0.09 per act (inf.prob.latent). Half the early-
#                         stage rate, representing residual transmission
#                         during clinical relapses of secondary lesions
#                         within the first year after infection. Matches
#                         CDC's "infectious syphilis" classification.
#   Incubating:           0 per act (inf.prob.incubating). Default: no
#                         transmission before the chancre appears, since
#                         the mucocutaneous-lesion route is by definition
#                         absent. Users can set > 0 to explore pre-chancre
#                         transmission as a sensitivity analysis.
#   Late latent, Tertiary: 0 per act (hard-coded in module-fx.R). Not
#                         sexually transmissible.
#
# Symptom recognition (probability of clinical detection at stage entry):
#   Primary:   20.5% -- many chancres are internal or go unnoticed
#   Secondary: 10.6% -- rash may be attributed to other causes
#   Tertiary:  100%  -- always symptomatic (set in module code)
#
#   This means ~80% of primary and ~89% of secondary infections are
#   clinically silent, detectable only through screening.
#
# Treatment:
#   Symptomatic treatment rate: 80% per week (early stages), 100% (tertiary)
#   Screening rate: 1/52 per week (~annual) for asymptomatic infected

param_scr <- param.net(
  inf.prob.incubating = 0,
  inf.prob.early = 0.18,
  inf.prob.latent = 0.09,
  act.rate = 2,
  ipr.rate = 1 / 4,
  prse.rate = 1 / 9,
  seel.rate = 1 / 17,
  elll.rate = 1 / 22,
  llter.rate = 1 / 1508,
  pri.sym = 0.205,
  sec.sym = 0.106,
  early.trt = 0.8,
  late.trt = 1.0,
  scr.rate = 1 / 52
)

# Initial conditions: 10 infected (1% prevalence)
init <- init.net(i.num = 10)

# Control settings:
#   type = NULL: fully custom module set
#   resimulate.network = FALSE: no vital dynamics, so the network is static
#     and does not need to be resimulated each step
#   Module order: infection -> progression -> treatment/screening -> prevalence
control <- control.net(
  type = NULL,
  nsteps = nsteps,
  nsims = nsims,
  ncores = ncores,
  infection.FUN = infect,
  progress.FUN = progress,
  tnt.FUN = tnt,
  resimulate.network = FALSE,
  module.order = c("resim_nets.FUN",
                   "infection.FUN",
                   "progress.FUN",
                   "tnt.FUN",
                   "prevalence.FUN")
)


# 3. Scenario 1: With Annual Screening ---------------------------------------

sim_scr <- netsim(est, param_scr, init, control)
print(sim_scr)


# 4. Scenario 2: No Screening (Baseline) -------------------------------------

# Same parameters but screening rate = 0. Without screening, only the ~20%
# of primary and ~11% of secondary cases that show recognizable symptoms
# can be detected and treated. The vast majority of infections progress
# silently to late latent, creating a large untreatable reservoir.
param_noscr <- param.net(
  inf.prob.incubating = 0,
  inf.prob.early = 0.18,
  inf.prob.latent = 0.09,
  act.rate = 2,
  ipr.rate = 1 / 4,
  prse.rate = 1 / 9,
  seel.rate = 1 / 17,
  elll.rate = 1 / 22,
  llter.rate = 1 / 1508,
  pri.sym = 0.205,
  sec.sym = 0.106,
  early.trt = 0.8,
  late.trt = 1.0,
  scr.rate = 0
)

sim_noscr <- netsim(est, param_noscr, init, control)
print(sim_noscr)


# 5. Analysis ----------------------------------------------------------------

# Compute derived measures.
#
# "infectious.num" follows CDC's surveillance definition of infectious
# syphilis: primary + secondary + early latent. With the default
# `inf.prob.incubating = 0`, incubating cases are not actively transmitting
# and so are excluded from this aggregate. "latent.num" covers the
# non-transmitting reservoir (incubating + late latent).
sim_scr <- mutate_epi(sim_scr,
  prev = i.num / num,
  infectious.num = pr.num + se.num + el.num,
  latent.num = inc.num + ll.num
)
sim_noscr <- mutate_epi(sim_noscr,
  prev = i.num / num,
  infectious.num = pr.num + se.num + el.num,
  latent.num = inc.num + ll.num
)


## --- Plot 1: Syphilis Prevalence Comparison ---
# Without screening, prevalence climbs toward saturation as the asymptomatic
# majority progresses to late latent and can never be detected or treated.
# With screening, prevalence stabilizes at a much lower endemic level.
par(mfrow = c(1, 1), mar = c(3, 3, 2, 1), mgp = c(2, 1, 0))
plot(sim_noscr, y = "prev",
     main = "Syphilis Prevalence: Screening vs. No Screening",
     ylab = "Prevalence", xlab = "Weeks",
     mean.col = "firebrick", mean.lwd = 2, mean.smooth = TRUE,
     qnts.col = "firebrick", qnts.alpha = 0.2, qnts.smooth = TRUE,
     legend = FALSE)
plot(sim_scr, y = "prev", add = TRUE,
     mean.col = "steelblue", mean.lwd = 2, mean.smooth = TRUE,
     qnts.col = "steelblue", qnts.alpha = 0.2, qnts.smooth = TRUE,
     legend = FALSE)
legend("topleft", legend = c("No Screening", "With Screening"),
       col = c("firebrick", "steelblue"), lwd = 2, bty = "n")


## --- Plot 2: Stage Composition (2-Panel) ---
# Left: Infectious stages (primary + secondary + early latent), matching
#   CDC's surveillance definition of "infectious syphilis". These cases
#   transmit at the rates governed by `inf.prob.early` (primary/secondary)
#   and `inf.prob.latent` (early latent).
# Right: Non-transmitting reservoir (incubating + late latent). With
#   default `inf.prob.incubating = 0`, incubating individuals are infected
#   but pre-chancre and therefore not yet transmitting. Late latent is
#   the larger long-term reservoir of asymptomatic, untreated infection.
par(mfrow = c(1, 2), mar = c(3, 3, 2, 1), mgp = c(2, 1, 0))

plot(sim_noscr, y = "infectious.num",
     main = "Infectious Stages (primary + secondary + early latent)",
     ylab = "Count", xlab = "Weeks",
     mean.col = "firebrick", mean.lwd = 2, mean.smooth = TRUE,
     qnts.col = "firebrick", qnts.alpha = 0.2, qnts.smooth = TRUE,
     legend = FALSE)
plot(sim_scr, y = "infectious.num", add = TRUE,
     mean.col = "steelblue", mean.lwd = 2, mean.smooth = TRUE,
     qnts.col = "steelblue", qnts.alpha = 0.2, qnts.smooth = TRUE,
     legend = FALSE)
legend("topleft", legend = c("No Screening", "With Screening"),
       col = c("firebrick", "steelblue"), lwd = 2, cex = 0.8, bty = "n")

plot(sim_noscr, y = "latent.num",
     main = "Non-transmitting (incubating + late latent)",
     ylab = "Count", xlab = "Weeks",
     mean.col = "firebrick", mean.lwd = 2, mean.smooth = TRUE,
     qnts.col = "firebrick", qnts.alpha = 0.2, qnts.smooth = TRUE,
     legend = FALSE)
plot(sim_scr, y = "latent.num", add = TRUE,
     mean.col = "steelblue", mean.lwd = 2, mean.smooth = TRUE,
     qnts.col = "steelblue", qnts.alpha = 0.2, qnts.smooth = TRUE,
     legend = FALSE)
legend("topleft", legend = c("No Screening", "With Screening"),
       col = c("firebrick", "steelblue"), lwd = 2, cex = 0.8, bty = "n")


## --- Summary Table ---

df_scr <- as.data.frame(sim_scr)
df_noscr <- as.data.frame(sim_noscr)
last_t <- max(df_scr$time)

# Summary statistics
data.frame(
  Metric = c("Mean prevalence",
             "Peak prevalence",
             "Cumulative new infections",
             "Total recoveries",
             "Peak latent reservoir"),
  Screening = c(round(mean(df_scr$prev, na.rm = TRUE), 3),
                round(max(df_scr$prev, na.rm = TRUE), 3),
                sum(df_scr$si.flow, na.rm = TRUE),
                sum(df_scr$rec.flow, na.rm = TRUE),
                max(df_scr$latent.num, na.rm = TRUE)),
  No_Screening = c(round(mean(df_noscr$prev, na.rm = TRUE), 3),
                   round(max(df_noscr$prev, na.rm = TRUE), 3),
                   sum(df_noscr$si.flow, na.rm = TRUE),
                   sum(df_noscr$rec.flow, na.rm = TRUE),
                   max(df_noscr$latent.num, na.rm = TRUE))
)

# Note: cumulative infections may be HIGHER with screening because recovered
# individuals return to the susceptible pool and can be reinfected. This is
# a feature of SIS-type dynamics. The relevant public health metric is
# prevalence (the disease burden at any point), not cumulative incidence.
