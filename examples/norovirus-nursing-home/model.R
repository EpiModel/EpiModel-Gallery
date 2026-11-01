
##
## Norovirus in a Nursing Home: Person-to-Person and Environmental Transmission
## EpiModel Gallery Example
##
## Description: An individual-based model of norovirus in a 96-bed US nursing
##   home over one winter season, in 8-hour shifts. Residents live in rooms on
##   three units; staff work shifts on a weekly roster; contact happens between
##   roommates, in each unit's dining and common areas, between staff on the
##   same team, and in care, and each unit has an environmental reservoir that
##   infectious people contaminate and everyone on the unit picks up from.
##   Norovirus enters with staff infected in the community and with new
##   admissions. When the facility declares an outbreak, control measures
##   switch on. Part 1 describes a season; Part 2 compares control measures;
##   Part 3 recalibrates the model at five environmental shares of
##   transmission that reproduce the same outbreaks, and asks whether the
##   ranking of the measures depends on the share.
##
## Author: Samuel M. Jenness
## Date: November 2026
##

# Load EpiModel
suppressMessages(library(EpiModel))

# Standard Gallery unit test lines
rm(list = ls())
eval(parse(text = print(commandArgs(TRUE)[1])))

if (interactive() || exists("run_full")) {
  nsims <- 1200
  ncores <- 10
  ndays <- 180
} else {
  nsims <- 3
  ncores <- 1
  ndays <- 60
}
nsteps <- 3 * ndays

source("examples/norovirus-nursing-home/module-fx.R")


# Facility ------------------------------------------------------------------------

# Three units of 32 beds: 13 double rooms (beds 1 to 26) and 6 single rooms
# (beds 27 to 32). Unit 1 is the short-stay rehabilitation unit.
n.units <- 3
beds.per.unit <- 32
res_unit <- rep(seq_len(n.units), each = beds.per.unit)
res_bed <- rep(seq_len(beds.per.unit), n.units)
room_in_unit <- ifelse(res_bed <= 26, (res_bed + 1) %/% 2, 13 + (res_bed - 26))
res_room <- (res_unit - 1) * 19 + room_in_unit

# Staff pools. Each unit has a pool of nursing assistants (CNAs) and of
# licensed nurses for each shift, sized so that the pool covers the shift's
# care groups every day of the week when each member takes two days off.
# Float staff (activities, therapy, dietary, housekeeping) work days across
# all units.
cna_pool <- c(day = 6, evening = 5, night = 3)
nurse_pool <- c(day = 3, evening = 2, night = 2)
n_float <- 14
staff <- do.call(rbind, lapply(seq_len(n.units), function(u) {
  do.call(rbind, lapply(1:3, function(sh) {
    rbind(data.frame(role = "cna", unit = u, staff_shift = sh,
                     pool_k = seq_len(cna_pool[sh])),
          data.frame(role = "nurse", unit = u, staff_shift = sh,
                     pool_k = seq_len(nurse_pool[sh])))
  }))
}))
staff <- rbind(staff, data.frame(role = "float", unit = 0, staff_shift = 1,
                                 pool_k = seq_len(n_float)))
# Two days off a week, staggered across each pool.
staff$off1 <- ((2 * staff$pool_k - 2) %% 7) + 1
staff$off2 <- ((2 * staff$pool_k - 1) %% 7) + 1
staff$team <- ifelse(staff$role == "float", "float",
                     paste0("u", staff$unit, "s", staff$staff_shift))
staff$float_unit <- ifelse(staff$role == "float",
                           (staff$pool_k - 1) %% n.units + 1, NA)

n_res <- length(res_unit)
n_staff <- nrow(staff)
n <- n_res + n_staff
table(staff$role, staff$staff_shift)


# Network with clique layers ------------------------------------------------------

# Every attribute any module reads is set on the base network, so netsim()
# copies it into the simulation. Residents have no staff attributes and staff
# no room or bed.
nw <- network_initialize(n)
attrs <- list(
  type = c(rep("res", n_res), rep("staff", n_staff)),
  role = c(rep("resident", n_res), staff$role),
  unit = c(res_unit, staff$unit),
  res_unit = c(res_unit, rep(NA, n_staff)),
  room = c(res_room, rep(NA, n_staff)),
  bed = c(res_bed, rep(NA, n_staff)),
  team = c(rep(NA, n_res), staff$team),
  staff_id = c(rep(NA, n_res), seq_len(n_staff)),
  staff_shift = c(rep(NA, n_res), staff$staff_shift),
  pool_k = c(rep(NA, n_res), staff$pool_k),
  off1 = c(rep(NA, n_res), staff$off1),
  off2 = c(rep(NA, n_res), staff$off2),
  float_unit = c(rep(NA, n_res), staff$float_unit),
  infRoute = rep(NA_character_, n),
  stage_end = rep(NA_real_, n),
  carer = rep(NA_real_, n),
  nurse_of = rep(NA_real_, n),
  onsetTime = rep(NA_real_, n),
  symEnd = rep(NA_real_, n),
  hold_until = rep(NA_real_, n),
  on_duty = rep(0, n),
  vomit = rep(0, n),
  presentee = rep(0, n),
  isolated = rep(0, n),
  excluded = rep(0, n)
)
for (a in names(attrs)) nw <- set_vertex_attribute(nw, a, attrs[[a]])

# Three clique layers: roommates, residents of the same unit, and staff on
# the same team. A new resident joins the room and unit of the bed they are
# admitted to, which the admission module writes to their attributes.
room_back <- function(dat, at, new_ids, network) get_attr(dat, "room")[new_ids]
unit_back <- function(dat, at, new_ids, network) get_attr(dat, "res_unit")[new_ids]
est_room <- netclique(nw, group.attr = "room", arrivals = "join",
                      arrivals.FUN = room_back)
est_unit <- netclique(nw, group.attr = "res_unit", arrivals = "join",
                      arrivals.FUN = unit_back)
est_team <- netclique(nw, group.attr = "team")
layers <- list(est_room, est_unit, est_team)


# Parameters ----------------------------------------------------------------------

# Vectors of three are values for the day, evening, and night shifts, and
# rel.inf.a holds values for residents and staff. Contact weights are hours of
# close contact per pair per shift, with hands-on care weighted up for contact
# with stool and vomit. beta and env.dep, the two
# transmission parameters, are set in Part 3 by the calibration chosen there.
param_base <- list(
  # facility and staffing
  n.units = n.units, beds.per.unit = beds.per.unit, cna.groups = c(4, 3, 2),
  los = c(40, 730, 730), fill.prob = 0.5,
  # natural history (days)
  inc.meanlog = 0.182, inc.sdlog = 0.495,
  sym.mean = 2, sym.shape = 4, shed.mean = 10, shed.shape = 2,
  p.asym = 0.3, p.vomit = 0.6, p.immune = 0.45,
  # person-to-person transmission
  beta = NA, rel.inf.a = c(0.05, 0.01), vomit.mult = 1.5,
  w.room = c(1, 1, 1), w.common = c(0.03, 0.02, 0), w.team = c(0.1, 0.1, 0.1),
  w.care = c(0.5, 0.35, 0.2), w.nurse = c(0.1, 0.05, 0.03),
  w.float = c(0.2, 0, 0), n.float.res = 6,
  # environment
  env.dep = NA, env.decay = 0.18, env.w.res = c(1, 1, 0.25),
  clean.routine = 0.3, clean.bleach = 0.6,
  # introductions (per staff member per day, and per admission, at the peak)
  intro.staff = 7.9e-5, intro.adm = 1.86e-3, season.peak = 76, season.sd = 45,
  # staff who work while ill
  p.present = 0.3,
  # outbreak declaration
  trig.n = 3, trig.window = 9, trig.delay = 3, end.window = 12, ob.gap = 42,
  # control measures (all off) and their settings
  ctl.isolate = 0, ctl.exclude = 0, ctl.cohort = 0, ctl.close = 0,
  ctl.clean = 0, ctl.handhyg = 0,
  iso.hold = 6, excl.hold = 6, prec.mult = 0.5, hh.mult = 0.8
)

# Five calibrations that reproduce the same outbreaks with different shares
# of transmission through the environment (Part 3). The page's main analysis
# uses the middle one.
calib <- data.frame(
  share = c(0.10, 0.20, 0.30, 0.45, 0.60),
  beta = c(0.935, 0.796, 0.670, 0.490, 0.317),
  env.dep = c(0.0520, 0.1010, 0.1485, 0.2205, 0.2894)
)
central <- 3
param_base$beta <- calib$beta[central]
param_base$env.dep <- calib$env.dep[central]


# Introductions -------------------------------------------------------------------

# Norovirus enters through staff exposed in the community and through new
# admissions, at rates that follow a winter season peaking in mid-January.
# Introductions are rare: at these rates most seasons have none, and a season
# without one has no cases under any policy. Each simulated season therefore
# draws its introductions before it runs, conditional on at least one, and
# season-level results are weighted by the probability of at least one,
# 1 - exp(-Lambda), where Lambda is the expected number per season.
adm_per_day <- sum(beds.per.unit * (1 - exp(-1 / param_base$los)))
draw_intros <- function(p, ndays) {
  season <- exp(-(seq_len(ndays) - p$season.peak)^2 / (2 * p$season.sd^2))
  rate <- c(p$intro.staff * n_staff * season, p$intro.adm * adm_per_day * season)
  Lambda <- sum(rate)
  k <- qpois(runif(1, ppois(0, Lambda), 1), Lambda)
  ev <- sample.int(length(rate), k, replace = TRUE, prob = rate)
  list(staff = ev[ev <= ndays], adm = ev[ev > ndays] - ndays,
       weight = 1 - exp(-Lambda))
}


# Control settings ----------------------------------------------------------------

control <- control.net(
  type = NULL, nsteps = nsteps, nsims = 1, ncores = 1,
  tergmLite = TRUE, resimulate.network = TRUE, save.nwstats = FALSE,
  staffing.FUN = staffing,
  infection.FUN = infect,
  environ.FUN = environ,
  progress.FUN = progress,
  outbreak.FUN = outbreak,
  departures.FUN = discharge,
  arrivals.FUN = admit,
  module.order = c("resim_nets.FUN", "summary_nets.FUN", "staffing.FUN",
                   "infection.FUN", "environ.FUN", "progress.FUN",
                   "outbreak.FUN", "departures.FUN", "arrivals.FUN",
                   "nwupdate.FUN"),
  save.other = "ob_log",
  verbose = FALSE
)


# Running a scenario --------------------------------------------------------------

# Each season is its own single-simulation netsim() call, run in parallel. A
# fixed future.seed gives season s the same random number stream in every
# scenario, so scenarios share their introductions and differ only once
# their control measures first act. Each season returns its totals, its
# outbreaks, and the daily onsets of each outbreak from its first case (NA
# past the end of the season).
routes <- c("room", "common", "team", "care", "nurse", "float")
run_scenario <- function(changes = list(), nsims, keep_epi = 0, seed = 2026) {
  p <- modifyList(param_base, changes)
  one <- function(s) {
    intro <- draw_intros(p, ndays)
    param <- do.call(param.net, c(p, list(intro.staff.days = intro$staff,
                                          intro.adm.days = intro$adm)))
    status0 <- ifelse(runif(n) < p$p.immune, "r", "s")
    sim <- netsim(layers, param, init.net(status.vector = status0), control)
    d <- as.data.frame(sim)
    d[is.na(d)] <- 0
    res_on <- d$onset.res.u1 + d$onset.res.u2 + d$onset.res.u3
    tot <- c(res_cases = sum(res_on), staff_cases = sum(d$onset.staff),
             colSums(d[, paste0("inf.", routes)]),
             inf.env.res = sum(d$inf.env.res), inf.env.staff = sum(d$inf.env.staff),
             intro_staff = sum(d$intro.staff), intro_adm = sum(d$intro.adm),
             staff_shifts_lost = sum(d$staff.absent),
             ob_days = sum(d$ob.active) / 3, admissions = sum(d$admissions),
             outbreaks = NROW(sim$ob_log[[1]]))
    # Person-to-person transmissions by route and by who infected whom.
    tm <- as.data.frame(get_transmat(sim))
    pairs <- c("res>res", "staff>res", "res>staff", "staff>staff")
    by_route <- paste0(rep(routes, each = 2), ".", c("res", "staff"))
    who <- setNames(numeric(length(pairs) + length(by_route)), c(pairs, by_route))
    if (nrow(tm) > 0) {
      who[pairs] <- table(factor(paste(tm$infType, tm$susType, sep = ">"),
                                 levels = pairs))
      who[by_route] <- table(factor(paste(tm$route, tm$susType, sep = "."),
                                    levels = by_route))
    }
    ob <- sim$ob_log[[1]]
    curves <- NULL
    if (NROW(ob) > 0) {
      day_on <- tapply(res_on + d$onset.staff, (d$time - 1) %/% 3 + 1, sum)
      curves <- t(sapply((ob$first - 1) %/% 3 + 1, function(d0) day_on[d0 + 0:39]))
      ob$season <- s
    }
    list(tot = c(tot, who, weight = intro$weight), ob = ob, curves = curves,
         epi = if (s <= keep_epi) d else NULL)
  }
  runs <- future.apply::future_lapply(seq_len(nsims), one, future.seed = seed)
  list(tot = as.data.frame(do.call(rbind, lapply(runs, `[[`, "tot"))),
       ob = do.call(rbind, lapply(runs, `[[`, "ob")),
       curves = do.call(rbind, lapply(runs, `[[`, "curves")),
       epi = lapply(runs, `[[`, "epi"))
}
if (ncores > 1) {
  future::plan("multisession", workers = ncores)
} else {
  future::plan("sequential")
}


# Analysis helpers ------------------------------------------------------------------

# A season total per 100 facility-seasons, with a 95% Monte Carlo interval.
per100 <- function(r, col) {
  x <- 100 * r$tot$weight * r$tot[[col]]
  m <- mean(x)
  h <- 1.96 * sd(x) / sqrt(length(x))
  c(mean = m, lo = m - h, hi = m + h)
}

# Percent change in a season total against a scenario run on the same random
# number streams, with a bootstrap interval that resamples seasons in pairs.
pct_change <- function(r, ref, col = "res_cases", B = 2000) {
  x <- r$tot[[col]]
  y <- ref$tot[[col]]
  est <- 100 * (sum(x) / sum(y) - 1)
  bs <- replicate(B, {
    i <- sample.int(length(x), replace = TRUE)
    100 * (sum(x[i]) / sum(y[i]) - 1)
  })
  c(est = est, lo = quantile(bs, 0.025, names = FALSE, na.rm = TRUE),
    hi = quantile(bs, 0.975, names = FALSE, na.rm = TRUE))
}

# One row per outbreak: size, attack rates among residents and staff of the
# units it touched (plus float staff), and duration from first to last onset.
staff_per_unit <- sum(staff$role != "float" & staff$unit == 1)
ob_stats <- function(r) {
  ob <- r$ob
  if (NROW(ob) == 0) {
    return(data.frame(size = numeric(0), ar_res = numeric(0), ar_staff = numeric(0),
                      days = numeric(0), units = numeric(0), resurged = logical(0)))
  }
  u <- pmax(ob$units, 1)
  data.frame(size = ob$cases.res + ob$cases.staff,
             ar_res = ob$cases.res / (beds.per.unit * u),
             ar_staff = ob$cases.staff / (staff_per_unit * u + n_float),
             days = (ob$last - ob$first) / 3,
             units = ob$units,
             resurged = ob$declarations > 1)
}

# Share of infections acquired in the facility that came from the environment.
env_share <- function(r) {
  env <- sum(r$tot$inf.env.res + r$tot$inf.env.staff)
  env / (env + sum(r$tot[, paste0("inf.", routes)]))
}

med_iqr <- function(x, digits = 0, pct = FALSE) {
  if (length(x) == 0) return("none")
  q <- quantile(x, c(0.5, 0.25, 0.75), names = FALSE) * if (pct) 100 else 1
  fmt <- paste0("%.", digits, "f [%.", digits, "f, %.", digits, "f]")
  sprintf(fmt, q[1], q[2], q[3])
}

# Severity among residents, per case (US long-term care outbreaks: 21.6
# hospitalizations and 2.3 deaths per 1,000 cases, 63% of cases residents).
p.hosp <- 0.034
p.death <- 0.0037


# Part 1: A Season with No Control Measures -------------------------------------

base <- run_scenario(nsims = nsims, keep_epi = min(nsims, 100))
ob_base <- ob_stats(base)

# The model against the targets.
targets <- data.frame(
  Quantity = c("Outbreaks per 100 facility-seasons",
               "Cases per outbreak, median [IQR]",
               "Resident attack rate on affected units (%)",
               "Staff attack rate (%)",
               "Residents as a share of cases (%)",
               "Duration, first to last onset (days)",
               "Units affected, mean"),
  Model = c(sprintf("%.1f", 100 * base$tot$weight[1] * mean(base$tot$outbreaks)),
            med_iqr(ob_base$size),
            med_iqr(ob_base$ar_res, pct = TRUE),
            med_iqr(ob_base$ar_staff, pct = TRUE),
            sprintf("%.0f", 100 * sum(base$tot$res_cases) /
                      sum(base$tot$res_cases + base$tot$staff_cases)),
            med_iqr(ob_base$days),
            sprintf("%.1f", mean(ob_base$units))),
  Reference = c("13 (Oregon, active surveillance); 4.1 (US reports, a lower bound)",
                "21 [11, 35] (US reports, 2015 to 2025)",
                "35.5 (Dutch wards, median)",
                "23.7 [11.6, 38.3] (published outbreaks); 10.9 (US reports)",
                "63 (107 US outbreaks)",
                "14.5 [9, 24.5] without measures (published); 8 to 13 (US reports)",
                "2.4 wards (Netherlands)")
)
print(targets)

# Burden per 100 facility-seasons.
burden <- rbind(
  "Resident cases" = per100(base, "res_cases"),
  "Resident hospitalizations" = p.hosp * per100(base, "res_cases"),
  "Resident deaths" = p.death * per100(base, "res_cases"),
  "Staff cases" = per100(base, "staff_cases"),
  "Staff shifts lost" = per100(base, "staff_shifts_lost")
)
round(burden, 1)

# Where infections happen: by route and by who was infected.
route_lab <- c(room = "Roommates", common = "Dining and common areas",
               team = "Staff team", care = "CNA care", nurse = "Nurse rounds",
               float = "Float staff", env = "Environment")
by_route <- rbind(
  Residents = c(colSums(base$tot[, paste0(routes, ".res")]), sum(base$tot$inf.env.res)),
  Staff = c(colSums(base$tot[, paste0(routes, ".staff")]), sum(base$tot$inf.env.staff))
)
colnames(by_route) <- route_lab
route_pct <- 100 * by_route / sum(by_route)
round(route_pct, 1)

who <- colSums(base$tot[, c("res>res", "staff>res", "res>staff", "staff>staff")])
round(100 * who / sum(who), 1)

par(mfrow = c(1, 1), mar = c(4, 11, 1, 1))
ord <- order(colSums(route_pct))
barplot(route_pct[, ord], horiz = TRUE, las = 1, col = c("#2a78d6", "#d9822b"),
        border = NA, xlab = "Share of infections acquired in the facility (%)")
legend("bottomright", fill = c("#2a78d6", "#d9822b"), border = NA, bty = "n",
       legend = c("Residents infected", "Staff infected"))

# One season with an outbreak of about median size and duration, day by day,
# from the seasons whose full output was kept.
kept <- which(!vapply(base$epi, is.null, logical(1)))
cand <- which(base$ob$season %in% kept)
if (length(cand) > 0) {
  dist <- abs(ob_base$size[cand] / median(ob_base$size) - 1) +
    abs(ob_base$days[cand] / median(ob_base$days) - 1)
  pick <- cand[which.min(dist)]
  ex <- base$epi[[base$ob$season[pick]]]
  ex$day <- (ex$time - 1) %/% 3 + 1
  daily <- aggregate(cbind(onset.res.u1, onset.res.u2, onset.res.u3, onset.staff,
                           ob.active) ~ day, data = ex, FUN = sum)
  env_day <- aggregate(cbind(env.u1, env.u2, env.u3) ~ day,
                       data = ex[ex$time %% 3 == 1, ], FUN = mean)
  first <- (base$ob$first[pick] - 1) %/% 3 + 1
  win <- daily$day >= first - 5 & daily$day <= first + 45
  x <- daily$day[win] - first
  bars <- as.matrix(daily[win, c("onset.res.u1", "onset.res.u2", "onset.res.u3",
                                 "onset.staff")])
  on <- daily$ob.active[win] > 0
  cols <- c("#2a78d6", "#1baf7a", "#6250d6", "#e34948")
  par(mfrow = c(2, 1), mar = c(1, 4.5, 1, 1), oma = c(3, 0, 0, 0))
  plot(NA, xlim = range(x), ylim = c(0, max(rowSums(bars), 1)), xaxt = "n",
       xlab = "", ylab = "Symptom onsets per day")
  if (any(on)) {
    rect(min(x[on]) - 0.5, 0, max(x[on]) + 0.5, par("usr")[4],
         col = adjustcolor("#e34948", 0.1), border = NA)
  }
  base_y <- rep(0, length(x))
  for (j in 1:4) {
    rect(x - 0.4, base_y, x + 0.4, base_y + bars[, j], col = cols[j], border = NA)
    base_y <- base_y + bars[, j]
  }
  legend("topright", fill = cols, border = NA, bty = "n", ncol = 2,
         legend = c("Unit 1 residents", "Unit 2 residents", "Unit 3 residents",
                    "Staff"))
  ew <- env_day$day >= first - 5 & env_day$day <= first + 45
  matplot(env_day$day[ew] - first, as.matrix(env_day[ew, -1]), type = "l", lty = 1,
          lwd = 2, col = cols[1:3], xlim = range(x), xlab = "",
          ylab = "Contamination (E)")
  mtext("Days from the first case", side = 1, line = 1.5, outer = TRUE)
}


# Part 2: Control Measures --------------------------------------------------------

bundle <- list(ctl.isolate = 1, ctl.exclude = 1, ctl.cohort = 1, ctl.close = 1,
               ctl.clean = 1)
arms <- list(
  isolate = list(ctl.isolate = 1),
  exclude = list(ctl.exclude = 1),
  cohort = list(ctl.cohort = 1),
  close = list(ctl.close = 1),
  clean = list(ctl.clean = 1),
  handhyg = list(ctl.handhyg = 1),
  bundle = bundle,
  early = c(bundle, list(trig.n = 2, trig.window = 6, trig.delay = 0)),
  sickleave = list(p.present = 0.1)
)
arm_lab <- c(isolate = "Isolate ill residents", exclude = "Exclude ill staff",
             cohort = "Cohort float staff", close = "Close units to admissions",
             clean = "Bleach twice a day", handhyg = "Soap and water after care",
             bundle = "All five measures", early = "All five, declared earlier",
             sickleave = "Paid sick leave")
runs <- lapply(arms, function(a) run_scenario(a, nsims))

# The hand hygiene effect is an assumption; two more settings bracket it.
hh_runs <- list(
  "0.95" = run_scenario(list(ctl.handhyg = 1, hh.mult = 0.95), nsims),
  "0.8" = runs$handhyg,
  "0.5" = run_scenario(list(ctl.handhyg = 1, hh.mult = 0.5), nsims)
)

ctl_table <- do.call(rbind, lapply(names(arms), function(a) {
  r <- runs[[a]]
  pc <- pct_change(r, base)
  ps <- pct_change(r, base, "staff_cases")
  data.frame(
    Scenario = arm_lab[[a]],
    res = per100(r, "res_cases")[["mean"]],
    res_pct = pc[["est"]], res_lo = pc[["lo"]], res_hi = pc[["hi"]],
    staff_pct = ps[["est"]],
    shifts = per100(r, "staff_shifts_lost")[["mean"]],
    outbreaks = 100 * r$tot$weight[1] * mean(r$tot$outbreaks),
    size = if (NROW(r$ob) > 0) median(ob_stats(r)$size) else NA
  )
}))
print(ctl_table, digits = 3)

par(mfrow = c(1, 1), mar = c(4, 13, 1, 1))
o <- order(ctl_table$res_pct, decreasing = TRUE)
y <- seq_along(o)
plot(ctl_table$res_pct[o], y, xlim = range(c(0, ctl_table$res_lo, ctl_table$res_hi),
                                         na.rm = TRUE), yaxt = "n", pch = 19,
     col = "#1a6b54", xlab = "Change in resident cases per season (%)", ylab = "")
segments(ctl_table$res_lo[o], y, ctl_table$res_hi[o], y, col = "#1a6b54", lwd = 2)
abline(v = 0, lty = 2, col = "gray50")
axis(2, at = y, labels = ctl_table$Scenario[o], las = 1, tick = FALSE)

hh_pct <- sapply(hh_runs, function(r) pct_change(r, base))
round(hh_pct, 1)


# Part 3: The Same Outbreaks, Different Routes -----------------------------------

ridge_arms <- c("isolate", "exclude", "cohort", "clean", "bundle")
ridge <- lapply(seq_len(nrow(calib)), function(k) {
  if (k == central) return(c(list(none = base), runs[ridge_arms]))
  cal <- list(beta = calib$beta[k], env.dep = calib$env.dep[k])
  out <- list(none = run_scenario(cal, nsims))
  for (a in ridge_arms) out[[a]] <- run_scenario(c(cal, arms[[a]]), nsims)
  out
})

# Each calibration against the outbreak targets.
ridge_fit <- do.call(rbind, lapply(seq_along(ridge), function(k) {
  r <- ridge[[k]]$none
  s <- ob_stats(r)
  data.frame(
    calibration = k, beta = calib$beta[k], env.dep = calib$env.dep[k],
    env_share = 100 * env_share(r),
    outbreaks = 100 * r$tot$weight[1] * mean(r$tot$outbreaks),
    res_cases = per100(r, "res_cases")[["mean"]],
    size = med_iqr(s$size), ar_res = med_iqr(s$ar_res, pct = TRUE),
    days = med_iqr(s$days)
  )
}))
ridge_fit

# The average outbreak, day by day from its first case, in each calibration.
ridge_cols <- c("#2a78d6", "#3b72a3", "#7f8e9c", "#c98a45", "#a8692a")
par(mfrow = c(1, 1), mar = c(4, 4.5, 1, 1))
curves <- lapply(ridge, function(x) {
  if (is.null(x$none$curves)) rep(NA, 40) else colMeans(x$none$curves, na.rm = TRUE)
})
ymax <- max(unlist(curves), 1, na.rm = TRUE)
plot(NA, xlim = c(0, 39), ylim = c(0, ymax * 1.05), xlab = "Days from the first case",
     ylab = "Symptom onsets per day, mean per outbreak")
for (k in seq_along(curves)) {
  if (all(is.na(curves[[k]]))) next
  lines(0:39, curves[[k]], col = ridge_cols[k], lwd = 2.5)
}
legend("topright", col = ridge_cols, lwd = 2.5, bty = "n",
       legend = sprintf("%.0f%% environmental", ridge_fit$env_share))

# How much each measure prevents, against the environmental share.
ridge_pct <- lapply(ridge_arms, function(a) {
  t(sapply(ridge, function(x) pct_change(x[[a]], x$none)))
})
names(ridge_pct) <- ridge_arms
arm_cols <- c(isolate = "#1a6b54", exclude = "#2a78d6", cohort = "#d9822b",
              clean = "#a8692a", bundle = "#46515b")
xs <- ridge_fit$env_share
par(mfrow = c(1, 1), mar = c(4, 4.5, 1, 1))
plot(NA, xlim = range(c(0, xs, 65), na.rm = TRUE),
     ylim = c(min(c(0, unlist(ridge_pct)), na.rm = TRUE) - 12, 0),
     xlab = "Share of transmission through the environment (%)",
     ylab = "Change in resident cases per season (%)")
abline(h = 0, lty = 2, col = "gray50")
for (a in ridge_arms) {
  m <- ridge_pct[[a]]
  polygon(c(xs, rev(xs)), c(m[, "lo"], rev(m[, "hi"])),
          col = adjustcolor(arm_cols[[a]], 0.15), border = NA)
  lines(xs, m[, "est"], col = arm_cols[[a]], lwd = 2.5, type = "b", pch = 19)
}
legend("bottom", col = arm_cols, lwd = 2.5, pch = 19, bty = "n", ncol = 3,
       legend = arm_lab[ridge_arms])
