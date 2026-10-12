
##
## Rabies on an Observed Raccoon Contact Network
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##
## Author: Samuel M. Jenness (Emory University)
## Date: September 2026
##

suppressMessages(library(EpiModel))

# Standard Gallery unit test lines
rm(list = ls())
eval(parse(text = print(commandArgs(TRUE)[1])))

# Two run modes:
#   interactive(): all four introduction weeks and three vaccination coverages,
#     200 simulations per scenario.
#   non-interactive (CI): two introduction weeks, two coverages, 10 simulations
#     per scenario, so the script completes in well under a minute.
if (interactive() || exists("run_full")) {
  nsims <- 200
  ncores <- 5
  intro.weeks <- c(4, 17, 30, 43)
  vax.covs <- c(0, 0.3, 0.6)
} else {
  nsims <- 10
  ncores <- 1
  intro.weeks <- c(4, 30)
  vax.covs <- c(0, 0.6)
}
# The observed year is played twice (see below), so a chain started late in
# year one can run into year two.
nsteps <- 104


# 1. The Observed Contact Network -------------------------------------------

# Reynolds et al. (2015) fitted proximity-logging collars to raccoons in a
# suburban forest preserve in Illinois and recorded every close contact
# (within 1 to 1.5 m, contacts under 1 s excluded) for a year. The data
# deposited with the paper are 52 weekly matrices of contact seconds between
# raccoon dyads; the weeks the authors mark as the breeding season are 24 to
# 39. raccoon_contacts.csv is those matrices in long form: one row per dyad
# and week with at least one contact, from the Animal Social Network
# Repository copy of the deposit. Raccoon ids run 1 to 24; a raccoon with no
# row in a week was not in contact with any other collared raccoon that week
# (in most cases its collar was not active).
contacts <- read.csv(if (file.exists("examples/rabies-raccoon-network/raccoon_contacts.csv"))
                       "examples/rabies-raccoon-network/raccoon_contacts.csv"
                     else "raccoon_contacts.csv")
n <- 24

cat(sprintf("%d dyad-weeks of contact among %d raccoons over %d weeks\n",
            nrow(contacts), n, max(contacts$week)))
cat("Contacts per week (raccoon pairs in contact):\n")
print(table(factor(contacts$week, levels = 1:52)))
cat("\nContact hours per dyad-week:\n")
print(round(quantile(contacts$seconds / 3600, c(0.1, 0.25, 0.5, 0.75, 0.9, 1)), 2))

# The year is observed once, but a rabies chain started in the autumn of that
# year has not run its course by week 52: with an incubation of a month or
# more, a year holds only a few generations. We therefore play the observed
# year twice, on the assumption that the seasonal contact pattern repeats,
# and follow every introduction for at least 52 weeks. This is the one
# modeling assumption layered on top of the data.
contacts2 <- rbind(contacts, transform(contacts, week = week + 52))

# A dynamic census. Each dyad-week becomes an edge spell covering that week,
# so that netsim's time step `at` reads the contacts of week `at`.
nw <- network_initialize(n)
nw <- set_vertex_attribute(nw, "id", 1:n)
nwd <- networkDynamic::networkDynamic(
  base.net = nw,
  edge.spells = data.frame(onset = contacts2$week,
                           terminus = contacts2$week + 1,
                           tail = contacts2$raccoon1,
                           head = contacts2$raccoon2),
  verbose = FALSE
)
nwd %n% "net.obs.period" <- list(observations = list(c(1, nsteps + 1)),
                                 mode = "discrete", time.increment = 1,
                                 time.unit = "week")

# The census layer. There is no ERGM here: the whole network was observed,
# so there is nothing to estimate and nothing for netdx to check.
obs <- netcensus(nwd)
obs

# Contact intensity. The layer is binary (in contact this week or not); the
# hours of contact for each pair and week go to the infection module as a
# parameter, one named vector per week.
contact.hours <- lapply(1:nsteps, function(w) {
  cw <- contacts2[contacts2$week == w, ]
  setNames(cw$seconds / 3600, paste(cw$raccoon1, cw$raccoon2, sep = "_"))
})


# 2. Parameters ---------------------------------------------------------------

# Time step: one week.
#
# inf.hazard: transmission hazard per hour of proximity between an infectious
#   and a susceptible raccoon, so the per-pair weekly transmission probability
#   is 1 - exp(-inf.hazard * hours). At 2 per hour, a pair with the median
#   contact time (about three minutes) transmits with probability 0.09, a
#   pair in contact for an hour with probability 0.86, and a pair that shares
#   a den for a day with probability 1. The value is illustrative; the
#   sensitivity section varies it.
# ei.rate: weekly rate of leaving the incubating state; 1/5 gives a mean
#   incubation of five weeks, in the range of a few weeks to a few months
#   reported for rabies in raccoons.
# Infectious raccoons transmit for one weekly step and die at its end.
# intro.week: the week in which one infectious raccoon is introduced.
# vax.cov: share of raccoons immune at baseline through oral rabies
#   vaccination.
param <- param.net(
  inf.hazard = 2,
  ei.rate = 1 / 5,
  intro.week = 4,
  vax.cov = 0,
  contact.hours = contact.hours
)

# Every raccoon starts susceptible; the introduction module seeds the
# infection at intro.week.
init <- init.net(i.num = 0)

source(if (file.exists("examples/rabies-raccoon-network/module-fx.R"))
         "examples/rabies-raccoon-network/module-fx.R"
       else "module-fx.R")

# The census layer needs no resimulation of its own, but under tergmLite the
# built-in resim_nets module is what reads the week's contacts from the
# observed network, so it stays in the pipeline and runs before infection.
control <- control.net(
  type = NULL,
  nsims = nsims,
  ncores = ncores,
  nsteps = nsteps,
  tergmLite = TRUE,
  vaccinate.FUN = vaccinate,
  introduce.FUN = introduce,
  infection.FUN = infect,
  progress.FUN = progress,
  prevalence.FUN = prevalence,
  module.order = c("vaccinate.FUN", "resim_nets.FUN", "summary_nets.FUN",
                   "introduce.FUN", "infection.FUN", "progress.FUN",
                   "nwupdate.FUN", "prevalence.FUN"),
  save.transmat = TRUE,
  save.nwstats = FALSE,
  verbose = FALSE
)


# 3. Scenarios: Introduction Week by Vaccination Coverage ---------------------

# Each scenario introduces one infectious raccoon in a given week into a
# population with a given vaccination coverage. The outcome is the number of
# raccoons infected beyond the introduced one by the end of the second year.
scenarios <- expand.grid(intro.week = intro.weeks, vax.cov = vax.covs)
sims <- vector("list", nrow(scenarios))
for (k in seq_len(nrow(scenarios))) {
  p <- param
  p$intro.week <- scenarios$intro.week[k]
  p$vax.cov <- scenarios$vax.cov[k]
  sims[[k]] <- netsim(obs, p, init, control)
}


# 4. Analysis -------------------------------------------------------------------

# Secondary cases per simulation: everyone who ended the run exposed,
# infectious, or dead, minus the introduced raccoon.
outbreak_size <- function(sim) {
  df <- as.data.frame(sim)
  last <- df[df$time == nsteps, ]
  with(last, e.num + i.num + d.num - 1)
}
scenarios$size <- lapply(sims, outbreak_size)
scenarios$mean.size <- sapply(scenarios$size, mean)
scenarios$p.any <- sapply(scenarios$size, function(x) mean(x >= 1))
scenarios$p.five <- sapply(scenarios$size, function(x) mean(x >= 5))

cat("\nSecondary cases by introduction week and vaccination coverage\n")
print(data.frame(intro.week = scenarios$intro.week,
                 vax.cov = scenarios$vax.cov,
                 mean.size = round(scenarios$mean.size, 2),
                 p.any = round(scenarios$p.any, 2),
                 p.five = round(scenarios$p.five, 2)))

# Unvaccinated scenarios: distribution of outbreak size by introduction week
base <- which(scenarios$vax.cov == 0)
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
boxplot(scenarios$size[base], names = paste("wk", scenarios$intro.week[base]),
        xlab = "Week of introduction", ylab = "Secondary cases by week 104",
        main = "No vaccination", col = "gray90")

# Vaccination: probability of five or more secondary cases by coverage, for
# each introduction week
p5 <- with(scenarios, tapply(p.five, list(vax.cov, intro.week), identity))
matplot(as.numeric(rownames(p5)), p5, type = "b", pch = 16, lty = 1,
        xlab = "Vaccination coverage", ylab = "P(5 or more secondary cases)",
        main = "By introduction week", ylim = c(0, 1))
legend("topright", legend = paste("wk", colnames(p5)), col = seq_len(ncol(p5)),
       pch = 16, lty = 1, bty = "n")

# Epidemic curves for the unvaccinated introduction weeks
par(mfrow = c(1, 1), mar = c(4, 4, 2, 1))
cols <- seq_along(base)
for (j in seq_along(base)) {
  df <- as.data.frame(sims[[base[j]]], out = "mean")
  if (j == 1) {
    plot(df$time, df$d.num, type = "l", col = cols[j], lwd = 2,
         ylim = c(0, n), xlab = "Week", ylab = "Cumulative rabies deaths",
         main = "Mean over simulations, no vaccination")
  } else {
    lines(df$time, df$d.num, col = cols[j], lwd = 2)
  }
}
legend("topleft", legend = paste("introduced wk", scenarios$intro.week[base]),
       col = cols, lwd = 2, bty = "n")

# Transmission chains: contact hours of the pairs that transmitted, pooled
# over the unvaccinated scenarios
tm <- do.call(rbind, lapply(base, function(k) {
  do.call(rbind, Filter(function(x) NROW(x) > 0, sims[[k]]$stats$transmat))
}))
if (!is.null(tm) && nrow(tm) > 0) {
  cat("\nContact hours of transmitting pairs (unvaccinated scenarios):\n")
  print(round(quantile(tm$hours, c(0.1, 0.25, 0.5, 0.75, 0.9)), 2))
  cat("Share of transmissions over pairs with at least one hour of contact:",
      round(mean(tm$hours >= 1), 2), "\n")
}
