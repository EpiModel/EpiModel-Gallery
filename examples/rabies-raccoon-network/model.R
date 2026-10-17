##
## Rabies on an Observed Raccoon Contact Network
## EpiModel Gallery (https://github.com/EpiModel/EpiModel-Gallery)
##
## Author: Samuel M. Jenness (Emory University)
## Date: September 2026
##

suppressMessages(library(EpiModel))

# Standard Gallery unit test lines
rm(list = ls())
eval(parse(text = print(commandArgs(TRUE)[1])))

# Two run modes:
#   interactive(): an introduction in every week of the year, 100 simulations
#     per week, at four levels of baseline immunity.
#   non-interactive (CI): four introduction weeks, about five simulations per
#     week, two levels of immunity, so the script completes in well under a
#     minute.
# Each simulation draws its introduction week at random from intro.weeks, so
# the number of simulations per week varies around nsims.week.
if (interactive() || exists("run_full")) {
  nsims.week <- 100
  ncores <- 5
  intro.weeks <- 2:52
  vax.covs <- c(0, 0.25, 0.5, 0.75)
} else {
  nsims.week <- 5
  ncores <- 1
  intro.weeks <- c(6, 19, 32, 45)
  vax.covs <- c(0, 0.5)
}
# The observed year is played twice (see below), so a chain started late in
# year one can run into year two.
nsteps <- 104
set.seed(2015)


# 1. The Observed Contact Network -------------------------------------------

# Hirsch et al. (2013) fitted proximity-logging collars to nearly all the
# resident adult raccoons in a 20-ha core of a suburban forest preserve in
# Illinois. Reynolds et al. (2015) deposited 52 weekly matrices of contact
# seconds (within 1 to 1.5 m, contacts under 1 s excluded) between 24 of
# them, starting in early July 2004; the weeks they mark as the breeding
# season are 24 to 39. raccoon_contacts.csv is those matrices in long form:
# one row per pair and week with at least one contact, from the Animal
# Social Network Repository copy of the deposit. A raccoon with no row in a
# week has no recorded contact that week, most often because its collar was
# not working.
contacts.raw <- read.csv(if (file.exists("examples/rabies-raccoon-network/raccoon_contacts.csv"))
                           "examples/rabies-raccoon-network/raccoon_contacts.csv"
                         else "raccoon_contacts.csv")
stopifnot(all(contacts.raw$raccoon1 < contacts.raw$raccoon2))
n <- 24

cat(sprintf("%d pair-weeks of contact among %d raccoons over %d weeks\n",
            nrow(contacts.raw), n, max(contacts.raw$week)))

# Weekly summaries: raccoons with any recorded contact, pairs in contact,
# total contact hours, mean hours per pair in contact, and the share of the
# possible pairs among the recorded raccoons that were in contact
wk <- split(contacts.raw, factor(contacts.raw$week, levels = 1:52))
weekly <- data.frame(week = 1:52,
                     animals = sapply(wk, function(d) length(unique(c(d$raccoon1, d$raccoon2)))),
                     pairs = sapply(wk, nrow),
                     hours = sapply(wk, function(d) sum(d$seconds) / 3600))
weekly$hours.per.pair <- weekly$hours / weekly$pairs
weekly$density <- weekly$pairs / choose(weekly$animals, 2)

# Four 13-week seasons, with week 1 in early July
season.of <- function(week) {
  cut(week, c(0, 13, 26, 39, 52),
      labels = c("summer (Jul-Sep)", "autumn (Oct-Dec)",
                 "winter (Jan-Mar)", "spring (Apr-Jun)"))
}
data.season <- round(sapply(weekly[, c("animals", "pairs", "density", "hours.per.pair")],
                            function(x) tapply(x, season.of(weekly$week), mean)), 2)
print(data.season)
cat("Correlation of pairs in contact with raccoons recorded:",
    round(cor(weekly$pairs, weekly$animals), 2), "\n")

par(mfrow = c(1, 3), mar = c(4, 4, 2, 1))
labs <- c(animals = "Raccoons recorded", pairs = "Pairs in contact",
          hours.per.pair = "Hours per pair in contact")
for (v in names(labs)) {
  y <- weekly[[v]]
  plot(1:52, y, type = "h", lwd = 3, col = "gray40", xlab = "Week",
       ylab = labs[v], main = labs[v], ylim = c(0, max(y)))
  rect(24, 0, 39, max(y) * 1.1, col = adjustcolor("orange", 0.2), border = NA)
}

# Three data problems. Weeks 3 and 4 record as many pairs as the weeks around
# them but almost no contact time, and week 4 is the first week of one
# collar and the last of two others. Weeks 19 and 20 are identical row for
# row. And three pair-weeks record exactly a full week (604,800 s) together.
print(round(weekly[1:8, c("week", "animals", "pairs", "hours")], 2))
w19 <- contacts.raw[contacts.raw$week == 19, -1]
w20 <- contacts.raw[contacts.raw$week == 20, -1]
cat("Weeks 19 and 20 identical:",
    isTRUE(all.equal(w19[order(w19$raccoon1, w19$raccoon2), ],
                     w20[order(w20$raccoon1, w20$raccoon2), ],
                     check.attributes = FALSE)), "\n")
print(contacts.raw[contacts.raw$seconds == 7 * 24 * 3600, ], row.names = FALSE)

# Weeks 3 and 4 are treated as missing and filled from the nearest observed
# week on either side (week 2 for week 3, week 5 for week 4). The duplicated
# weeks 19 and 20 are kept: they look like real contact data, repeated.
contacts <- rbind(contacts.raw[!contacts.raw$week %in% c(3, 4), ],
                  transform(contacts.raw[contacts.raw$week == 2, ], week = 3),
                  transform(contacts.raw[contacts.raw$week == 5, ], week = 4))
contacts <- contacts[order(contacts$week, contacts$raccoon1, contacts$raccoon2), ]

# The year is observed once, but a rabies chain started late in that year
# has not run its course by week 52: with an incubation period of a month or
# more, a year holds only a few generations. We therefore play the observed
# year again as year two, on the assumption that the seasonal contact
# pattern repeats. Year two includes only the raccoons still recorded at the
# end of year one; the others (collars that failed for good, or animals that
# died or left) do not come back.
last.week <- sapply(1:n, function(i) {
  max(contacts$week[contacts$raccoon1 == i | contacts$raccoon2 == i])
})
stayers <- which(last.week == 52)
cat("Raccoons in year two:", length(stayers), "of", n, "\n")
year2 <- contacts[contacts$raccoon1 %in% stayers & contacts$raccoon2 %in% stayers, ]
contacts2 <- rbind(contacts, transform(year2, week = week + 52))

# A dynamic census. Each pair-week becomes an edge spell covering that week,
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

# The census layer. There is no ERGM here: the network was observed, so
# there is nothing to estimate and nothing for netdx to check.
obs <- netcensus(nwd)
obs

# The census layer works with EpiModel's built-in models as it is: a quick
# check with the built-in SI model and a fixed probability per contact.
sim_si <- netsim(obs, param.net(inf.prob = 0.1, act.rate = 1), init.net(i.num = 2),
                 control.net(type = "SI", nsims = 1, nsteps = 52, tergmLite = TRUE,
                             resimulate.network = TRUE, verbose = FALSE))
print(tail(as.data.frame(sim_si)[, c("time", "s.num", "i.num")], 1), row.names = FALSE)

# Contact intensity. The layer is binary (in contact this week or not); the
# hours of contact for each pair and week go to the infection module as a
# parameter, one named vector per week, with the smaller id first in each
# name.
contact.hours <- lapply(1:nsteps, function(w) {
  cw <- contacts2[contacts2$week == w, ]
  setNames(cw$seconds / 3600, paste(pmin(cw$raccoon1, cw$raccoon2),
                                    pmax(cw$raccoon1, cw$raccoon2), sep = "_"))
})


# 2. Parameters ---------------------------------------------------------------

# Time step: one week.
#
# inf.hazard: transmission hazard per hour of proximity between an infectious
#   and a susceptible raccoon, so the per-pair weekly transmission
#   probability is 1 - exp(-inf.hazard * hours). Rabies is transmitted by
#   bites, and proximity time stands in for the opportunity; at 0.5 per hour
#   a three-minute contact transmits about one time in forty, an hour about
#   two times in five, and a day spent resting together almost always. The
#   value is illustrative; the sensitivity section varies it.
# inc.mean, inc.shape: the incubation period, weeks from infection to the
#   week in which a raccoon is infectious, is drawn at infection from a gamma
#   distribution with mean 6 weeks and shape 6, rounded to whole weeks. Its
#   mode is 5 weeks (Tinline et al. 2002) and 94% of draws fall between 3 and
#   12 weeks, the range reported for raccoons (Hirsch et al. 2013).
# Infectious raccoons transmit for one weekly step and die at its end.
# intro.week: the candidate weeks for the one introduced infectious raccoon;
#   each simulation draws one.
# vax.cov: share of raccoons immune at baseline through oral rabies
#   vaccination.
param <- param.net(
  inf.hazard = 0.5,
  inc.mean = 6,
  inc.shape = 6,
  intro.week = intro.weeks,
  vax.cov = 0,
  contact.hours = contact.hours
)

# Every raccoon starts susceptible; the introduction module seeds the
# infection at the simulation's introduction week.
init <- init.net(i.num = 0)

# Transmission potential by week, from the data alone. An infectious raccoon
# infects each of its contacts that week with probability
# 1 - exp(-inf.hazard * hours). For each raccoon in contact that week, the
# sum of those probabilities is the expected number it infects directly, and
# one minus the product of their complements is the probability that it
# infects at least one. Averaged over the raccoons in contact that week, as
# the introduction module picks the introduced animal, these are the
# expected direct infections and the chance of any onward transmission from
# an introduction in that week.
potential <- t(sapply(1:52, function(w) {
  cw <- contacts[contacts$week == w, ]
  p <- 1 - exp(-param$inf.hazard * cw$seconds / 3600)
  ids <- unique(c(cw$raccoon1, cw$raccoon2))
  mine <- lapply(ids, function(i) p[cw$raccoon1 == i | cw$raccoon2 == i])
  c(direct = mean(sapply(mine, sum)),
    any = mean(sapply(mine, function(q) 1 - prod(1 - q))))
}))
par(mfrow = c(1, 1), mar = c(4, 4, 2, 1))
plot(1:52, potential[, "direct"], type = "h", lwd = 4, col = "gray40",
     xlab = "Week", ylab = "Expected raccoons infected directly",
     ylim = c(0, max(potential[, "direct"]) * 1.1))
rect(24, 0, 39, max(potential[, "direct"]) * 1.1, col = adjustcolor("orange", 0.2),
     border = NA)
abline(h = 1, lty = 3, col = "gray50")
abline(v = c(13.5, 26.5, 39.5), lty = 2, col = "gray70")
cat("\nExpected direct infections from an introduction, by week:\n")
print(round(setNames(potential[, "direct"], 1:52), 2))

source(if (file.exists("examples/rabies-raccoon-network/module-fx.R"))
         "examples/rabies-raccoon-network/module-fx.R"
       else "module-fx.R")

# The census layer needs no resimulation of its own, but under tergmLite the
# built-in resim_nets module is what reads the week's contacts from the
# observed network, so it stays in the pipeline and runs before infection.
control <- control.net(
  type = NULL,
  nsims = nsims.week * length(intro.weeks),
  ncores = ncores,
  nsteps = nsteps,
  tergmLite = TRUE,
  resimulate.network = TRUE,
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


# 3. Simulations ---------------------------------------------------------------

# One netsim() call per level of baseline immunity. Each simulation draws
# its introduction week from intro.weeks, so one call covers every week. The
# outcome is the number of additional cases: raccoons exposed, infectious,
# or dead by week 104, minus the introduced one.
sims <- lapply(vax.covs, function(v) {
  p <- param
  p$vax.cov <- v
  netsim(obs, p, init, control)
})
names(sims) <- vax.covs


# 4. Analysis -------------------------------------------------------------------

# One row per simulation: introduction week, season, additional cases, and
# the attack rate among the other unvaccinated raccoons
outcomes <- function(sim) {
  df <- as.data.frame(sim)
  last <- df[df$time == nsteps, ]
  out <- data.frame(week = last$intro.week,
                    size = last$e.num + last$i.num + last$d.num - 1,
                    others = last$num - last$r.num - 1)
  out$season <- season.of(out$week)
  out$attack <- out$size / out$others
  out
}
res <- lapply(sims, outcomes)

# Means and proportions with 95% intervals over simulations
ci_mean <- function(x) {
  m <- mean(x)
  se <- sd(x) / sqrt(length(x))
  sprintf("%.2f (%.2f, %.2f)", m, m - 1.96 * se, m + 1.96 * se)
}
ci_prop <- function(x) {
  p <- mean(x)
  se <- sqrt(p * (1 - p) / length(x))
  sprintf("%.2f (%.2f, %.2f)", p, max(0, p - 1.96 * se), min(1, p + 1.96 * se))
}
summarize_runs <- function(d) {
  data.frame(Simulations = nrow(d),
             `P(any additional case)` = ci_prop(d$size >= 1),
             `Mean additional cases` = ci_mean(d$size),
             `Mean if any` = round(mean(d$size[d$size >= 1]), 2),
             `P(5 or more)` = ci_prop(d$size >= 5),
             check.names = FALSE)
}

## Season of introduction, no vaccination
r0 <- res[["0"]]
tab_season <- do.call(rbind, lapply(split(r0, r0$season, drop = TRUE), summarize_runs))
tab_season <- cbind(Season = rownames(tab_season), tab_season)
cat("\nAdditional cases by season of introduction, no vaccination\n")
print(tab_season, row.names = FALSE)
season.stats <- sapply(split(r0$size, r0$season), function(x) {
  c(p.any = mean(x >= 1), mean = mean(x), p.five = mean(x >= 5))
})

by.week <- data.frame(
  week = intro.weeks,
  p.any = sapply(intro.weeks, function(w) mean(r0$size[r0$week == w] >= 1)),
  mean = sapply(intro.weeks, function(w) mean(r0$size[r0$week == w])),
  se = sapply(intro.weeks, function(w) sd(r0$size[r0$week == w]) / sqrt(sum(r0$week == w)))
)
cat("Correlation of simulated P(any) by week with the expectation from the data:",
    round(cor(by.week$p.any, potential[intro.weeks, "any"], use = "complete.obs"), 2), "\n")

par(mfrow = c(2, 1), mar = c(4, 4, 1, 1))
plot(intro.weeks, by.week$p.any, pch = 16, col = "#2c3e50", ylim = c(0, 1),
     xlim = c(1, 52), xlab = "Week of introduction", ylab = "P(any additional case)")
lines(intro.weeks, potential[intro.weeks, "any"], col = "#e74c3c", lwd = 2)
abline(v = c(13.5, 26.5, 39.5), lty = 2, col = "gray70")
legend("topleft", c("simulated", "expected from the data"), pch = c(16, NA),
       lty = c(NA, 1), lwd = c(NA, 2), col = c("#2c3e50", "#e74c3c"), bty = "n")
plot(intro.weeks, by.week$mean, pch = 16, col = "#2c3e50", xlim = c(1, 52),
     ylim = c(0, max(by.week$mean + 1.96 * by.week$se, na.rm = TRUE)),
     xlab = "Week of introduction", ylab = "Mean additional cases")
ok <- which(by.week$se > 0)
arrows(intro.weeks[ok], by.week$mean[ok] - 1.96 * by.week$se[ok], intro.weeks[ok],
       by.week$mean[ok] + 1.96 * by.week$se[ok], angle = 90, code = 3, length = 0.02,
       col = "gray50")
segments(c(1, 14, 27, 40), season.stats["mean", ], c(13, 26, 39, 52),
         season.stats["mean", ], col = "#e74c3c", lwd = 3)
abline(v = c(13.5, 26.5, 39.5), lty = 2, col = "gray70")

cols <- c("#3498db", "#e74c3c", "#f39c12", "#27ae60")
df0 <- as.data.frame(sims[["0"]])
curves <- tapply(df0$d.num, list(df0$time, season.of(df0$intro.week)), mean)
par(mfrow = c(1, 1), mar = c(4, 4, 2, 1))
matplot(as.numeric(rownames(curves)), curves, type = "l", lty = 1, lwd = 2, col = cols,
        xlab = "Week", ylab = "Cumulative rabies deaths")
abline(v = 52, lty = 3, col = "gray50")
legend("topleft", legend = paste("introduced in", colnames(curves)), col = cols,
       lwd = 2, bty = "n")

## What the transmission assumption implies: contact hours of the pairs that
## transmitted, pooled over the unvaccinated simulations
tm <- do.call(rbind, Filter(function(x) NROW(x) > 0, sims[["0"]]$stats$transmat))
q <- round(quantile(tm$hours, c(0.1, 0.25, 0.5, 0.75, 0.9)), 2)
cat("\nContact hours of transmitting pairs (no vaccination):\n")
print(q)
cat("Share of transmissions over pairs with at least one hour of contact:",
    round(mean(tm$hours >= 1), 2), "\n")
cat("Share of pair-weeks in the data with at least one hour of contact:",
    round(mean(contacts$seconds >= 3600), 2), "\n")
cat("Share of transmissions in year two:", round(mean(tm$at > 52), 2), "\n")

## Baseline immunity, introduction week drawn over the year
tab_vax <- do.call(rbind, lapply(names(res), function(v) {
  d <- res[[v]]
  k <- round(as.numeric(v) * n)
  cbind(`Raccoons immune` = sprintf("%d of %d (%.0f%%)", k, n, 100 * k / n),
        summarize_runs(d)[, -1],
        `Attack rate, unvaccinated others` = ci_mean(d$attack))
}))
cat("\nAdditional cases by baseline immunity, introduction week drawn over the year\n")
print(tab_vax, row.names = FALSE)
vax.stats <- sapply(res, function(d) {
  c(p.any = mean(d$size >= 1), mean = mean(d$size), p.five = mean(d$size >= 5),
    attack = mean(d$attack))
})

p5 <- sapply(res, function(d) tapply(d$size >= 5, d$season, mean))
par(mfrow = c(1, 1), mar = c(4, 4, 2, 1))
matplot(round(as.numeric(colnames(p5)) * n) / n, t(p5), type = "b", pch = 16, lty = 1,
        lwd = 2, col = cols, xlab = "Share of raccoons immune at baseline",
        ylab = "P(5 or more additional cases)", ylim = c(0, max(p5, na.rm = TRUE) * 1.1))
legend("topright", legend = paste("introduced in", rownames(p5)), col = cols,
       pch = 16, lty = 1, bty = "n")


# 5. Sensitivity: the Transmission Hazard ------------------------------------

# Mean additional cases by season of introduction for other hazards per
# contact hour, with no vaccination and a quarter of the simulations
hazards <- c(0.25, 1, 2)
control_sens <- control
control_sens$nsims <- max(length(intro.weeks), control$nsims / 4)
sens <- lapply(hazards, function(h) {
  p <- param
  p$inf.hazard <- h
  d <- outcomes(netsim(obs, p, init, control_sens))
  c(tapply(d$size, d$season, mean), `P(any), all weeks` = mean(d$size >= 1))
})
tab_sens <- rbind(do.call(rbind, sens),
                  c(tapply(r0$size, r0$season, mean), mean(r0$size >= 1)))
tab_sens <- cbind(inf.hazard = c(hazards, param$inf.hazard), tab_sens)
tab_sens <- tab_sens[order(tab_sens[, "inf.hazard"]), ]
cat("\nMean additional cases by season of introduction and transmission hazard\n")
print(round(tab_sens, 2))
