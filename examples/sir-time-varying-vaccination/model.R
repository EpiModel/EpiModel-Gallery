
##
## SIR Model with Time-Varying (Phased) Vaccination
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##
## Author: Samuel M. Jenness (Emory University)
## Date: May 2026
##

# Load EpiModel
suppressMessages(library(EpiModel))

# Standard Gallery unit test lines
rm(list = ls())
eval(parse(text = print(commandArgs(TRUE)[1])))

if (interactive()) {
  nsims <- 5
  ncores <- 5
  nsteps <- 400
  nsteps_endemic <- 1000
} else {
  nsims <- 1
  ncores <- 1
  nsteps <- 50
  nsteps_endemic <- 80
}


# 1. Network Model Estimation ------------------------------------------------

# A simple 500-node contact network. The intervention timing is the
# pedagogical focus, so the network structure is kept intentionally simple.
n <- 500
nw <- network_initialize(n)

# Mean degree = 1.5 (375 edges in a 500-node population).
formation <- ~edges
target.stats <- c(round(1.5 * n / 2))

# Closed population, partnership duration ~60 timesteps.
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 60)
coef.diss

# Fit the ERGM
est <- netest(nw, formation, target.stats, coef.diss)

# Diagnostics
dx <- netdx(est, nsims = nsims, ncores = ncores, nsteps = nsteps,
            nwstats.formula = ~edges + degree(0:3))
print(dx)
plot(dx)


# 2. Epidemic Model Parameters -----------------------------------------------

source("examples/sir-time-varying-vaccination/module-fx.R")

# Disease parameters tuned for a clean single-wave epidemic on the closed
# SIR scenarios, peaking around timestep 60-80 at ~10-15% prevalence.
#
#   inf.prob = 0.10, act.rate = 1, rec.rate = 0.04 (mean inf. dur. ~25)
#   R0 (network-effective) is roughly 2-3.
#
# The endemic scenario adds:
#   rs.rate  -- waning of natural immunity (R -> S)
#   vax.wane -- waning of vaccine immunity (V -> S)
#   These keep the susceptible pool replenished, allowing the reactive
#   program to demonstrate sustained on/off cycling over a long horizon.

init <- init.net(i.num = 10)

# Module order is left to the default. EpiModel inserts our custom
# infection.FUN, recovery.FUN, and vaccinate.FUN slots alongside the
# built-in resim_nets / summary_nets / nwupdate / prevalence modules.
control <- control.net(
  type = NULL,
  nsims = nsims,
  ncores = ncores,
  nsteps = nsteps,
  infection.FUN = infect,
  recovery.FUN = recov,
  vaccinate.FUN = vaccinate,
  verbose = FALSE
)

# Base parameter set. vax.starts and vax.ends are length-2 vectors so
# the pulse scenario (two activation windows) can be expressed via the
# scenarios API's `_N` vector-position suffix. Single-window scenarios
# fill the second position with -1, which the vaccinate module
# interprets as "unused" (since at <= -1 is always FALSE).
param_base <- param.net(
  inf.prob = 0.10,
  act.rate = 1,
  rec.rate = 0.04,
  rs.rate = 0,
  exog.inf.prob = 0,
  vax.rate = 0.05,
  vax.efficacy = 0.9,
  vax.wane = 0,
  vax.starts = c(-1, -1),
  vax.ends = c(-1, -1),
  vax.prev.on = -1,
  vax.prev.off = -1
)


# 3. Closed-SIR Scenarios ----------------------------------------------------

# Five scenarios share the same disease parameters; only the
# vaccination schedule differs.
#   none    no vaccination
#   early   one window starting at step 30
#   late    one window starting at step 100
#   pulse   two windows (40-70 and 120-150)
#   react   prevalence-triggered with hysteresis (5% on / 2% off)
scenarios_sir.df <- data.frame(
  .scenario.id   = c("none", "early", "late",  "pulse", "react"),
  .at            = 0,
  vax.rate       = c(0,      0.05,    0.05,    0.05,    0.05),
  vax.starts_1   = c(-1,     30,      100,     40,      -1),
  vax.starts_2   = c(-1,     -1,      -1,      120,     -1),
  vax.ends_1     = c(-1,     nsteps,  nsteps,  70,      -1),
  vax.ends_2     = c(-1,     -1,      -1,      150,     -1),
  vax.prev.on    = c(-1,     -1,      -1,      -1,      0.05),
  vax.prev.off   = c(-1,     -1,      -1,      -1,      0.02)
)
scenarios_sir.list <- create_scenario_list(scenarios_sir.df)

sims <- list()
for (scn in scenarios_sir.list) {
  cat(sprintf("Scenario: %s\n", scn$id))
  sims[[scn$id]] <- netsim(est, use_scenario(param_base, scn),
                           init, control)
  print(sims[[scn$id]])
}
sim_none <- sims[["none"]]
sim_early <- sims[["early"]]
sim_late <- sims[["late"]]
sim_pulse <- sims[["pulse"]]
sim_react <- sims[["react"]]


# 4. Endemic Scenarios (SIRS + Waning Vaccine) -------------------------------

# Switching on waning of natural immunity (rs.rate) and vaccine immunity
# (vax.wane) converts the model to a setting that sustains an endemic
# equilibrium rather than burning out after a single wave. Run for a
# longer horizon to make the long-run behaviour visible.

control_endemic <- control.net(
  type = NULL,
  nsims = nsims,
  ncores = ncores,
  nsteps = nsteps_endemic,
  infection.FUN = infect,
  recovery.FUN = recov,
  vaccinate.FUN = vaccinate,
  verbose = FALSE
)

scenarios_endemic.df <- data.frame(
  .scenario.id = c("endemic_none", "endemic_react"),
  .at          = 0,
  vax.rate     = c(0,      0.07),
  rs.rate      = c(0.015,  0.015),
  vax.wane     = c(0.03,   0.03),
  vax.prev.on  = c(-1,     0.05),
  vax.prev.off = c(-1,     0.02)
)
scenarios_endemic.list <- create_scenario_list(scenarios_endemic.df)

sims_endemic <- list()
for (scn in scenarios_endemic.list) {
  cat(sprintf("Scenario: %s\n", scn$id))
  sims_endemic[[scn$id]] <- netsim(est, use_scenario(param_base, scn),
                                   init, control_endemic)
  print(sims_endemic[[scn$id]])
}
sim_endemic_none <- sims_endemic[["endemic_none"]]
sim_endemic_react <- sims_endemic[["endemic_react"]]


# 5. Analysis ----------------------------------------------------------------

labels <- c(none = "No vaccination",
            early = "Early (start 30)",
            late = "Late (start 100)",
            pulse = "Pulses (40-70, 120-150)",
            react = "Reactive (>5% / <2%)")
cols <- c(none = "gray40", early = "seagreen", late = "firebrick",
          pulse = "#8e44ad", react = "steelblue")

sims <- lapply(sims, function(s) mutate_epi(s, prev = i.num / num))


## --- Plot 1: Prevalence Overlay (Closed-SIR Scenarios) ---
# Same disease, same per-step vaccination rate, different activation
# schedules. The Early and Reactive curves clip the epidemic peak; Late
# only mops up after the peak; Pulses produce a knock-down/rebound shape.
par(mfrow = c(1, 1), mar = c(3, 3, 2, 1), mgp = c(2, 1, 0))
plot(sims$none, y = "prev",
     main = "Prevalence by Intervention Timing",
     ylab = "Prevalence (I / N)", xlab = "Time Step",
     mean.col = cols["none"], mean.lwd = 2, mean.smooth = TRUE,
     qnts = FALSE, legend = FALSE)
for (k in c("early", "late", "pulse", "react")) {
  plot(sims[[k]], y = "prev",
       mean.col = cols[k], mean.lwd = 2, mean.smooth = TRUE,
       qnts = FALSE, legend = FALSE, add = TRUE)
}
legend("topright", legend = labels, col = cols, lwd = 2, bty = "n",
       cex = 0.8)


## --- Plot 2: Disease Burden + Per-Dose Efficiency ---
cum_inf <- sapply(sims, function(s) {
  df <- as.data.frame(s)
  mean(tapply(df$si.flow, df$sim, sum, na.rm = TRUE))
})
total_vax <- sapply(sims, function(s) {
  df <- as.data.frame(s)
  mean(tapply(df$vax.flow, df$sim, sum, na.rm = TRUE))
})
averted <- cum_inf["none"] - cum_inf
avert_per_dose <- ifelse(total_vax > 0, averted / total_vax, NA)

par(mfrow = c(1, 2), mar = c(7, 4, 2, 1), mgp = c(2.5, 1, 0))
bp <- barplot(cum_inf, col = cols, las = 2,
              ylab = "Cumulative infections",
              main = "Disease Burden", names.arg = labels,
              cex.names = 0.8, ylim = c(0, max(cum_inf) * 1.15))
text(bp, cum_inf + max(cum_inf) * 0.03, round(cum_inf),
     cex = 0.8, font = 2)

apd_plot <- avert_per_dose
apd_plot[is.na(apd_plot)] <- 0
bp2 <- barplot(apd_plot, col = cols, las = 2,
               ylab = "Infections averted per dose",
               main = "Per-Dose Efficiency", names.arg = labels,
               cex.names = 0.8,
               ylim = c(0, max(apd_plot, na.rm = TRUE) * 1.15 + 0.01))
text(bp2, apd_plot + max(apd_plot, na.rm = TRUE) * 0.03,
     ifelse(is.na(avert_per_dose), "n/a", sprintf("%.2f", apd_plot)),
     cex = 0.8, font = 2)


## --- Plot 3: Endemic Counterfactual vs. Reactive Intervention ---
# The headline phenomenon. Top: same SIRS dynamics, no vaccination, the
# disease settles into a sustained endemic equilibrium that oscillates
# around the activation threshold. Bottom: with reactive vaccination, the
# program toggles on and off many times across the horizon, suppressing
# each rising wave below the threshold and standing down when prevalence
# falls. A single representative simulation is plotted so the binary
# on/off activity flag is visible (averaging across sims smooths the flag
# into a fractional value and hides the cycling).
sim_endemic_none <- mutate_epi(sim_endemic_none, prev = i.num / num)
sim_endemic_react <- mutate_epi(sim_endemic_react, prev = i.num / num)
df_n1 <- subset(as.data.frame(sim_endemic_none), sim == 1)
df_r1 <- subset(as.data.frame(sim_endemic_react), sim == 1)
ymax <- max(df_n1$i.num / df_n1$num, na.rm = TRUE) * 1.2

par(mfrow = c(2, 1), mar = c(3, 4, 2, 4), mgp = c(2.5, 1, 0))
plot(df_n1$time, df_n1$prev, type = "l", lwd = 2, col = "firebrick",
     xlab = "Time Step", ylab = "Prevalence (I / N)",
     main = "Endemic SIRS without Intervention (sim 1)",
     ylim = c(0, ymax))
abline(h = 0.05, lty = 2, col = "firebrick", lwd = 0.6)

plot(df_r1$time, df_r1$prev, type = "l", lwd = 2, col = "steelblue",
     xlab = "Time Step", ylab = "Prevalence (I / N)",
     main = "Endemic SIRS + Reactive Vaccination (sim 1)",
     ylim = c(0, ymax))
abline(h = 0.05, lty = 2, col = "firebrick")
abline(h = 0.02, lty = 2, col = "seagreen")
par(new = TRUE)
plot(df_r1$time, df_r1$vax.active, type = "s", lwd = 1.5,
     col = "gray40", axes = FALSE, xlab = "", ylab = "",
     ylim = c(0, 1.05))
axis(4)
mtext("Vaccination active (0/1)", side = 4, line = 2.5)
legend("topright",
       legend = c("Prevalence", "Activation 5%",
                  "Deactivation 2%", "Program active"),
       col = c("steelblue", "firebrick", "seagreen", "gray40"),
       lty = c(1, 2, 2, 1), lwd = c(2, 1, 1, 1.5), bty = "n",
       cex = 0.75)


## --- Summary Table ---
final_v <- sapply(sims, function(s) {
  df <- as.data.frame(s)
  round(mean(df$v.num[df$time == max(df$time)], na.rm = TRUE))
})
pct_averted <- ifelse(names(averted) == "none", "--",
                      paste0(round(averted / cum_inf["none"] * 100), "%"))
apd_str <- ifelse(is.na(avert_per_dose), "--",
                  sprintf("%.2f", avert_per_dose))

summary_tbl <- data.frame(
  Scenario = labels,
  Cum_inf = round(cum_inf),
  Averted = ifelse(names(averted) == "none", 0, round(averted)),
  Pct_averted = pct_averted,
  Doses = round(total_vax),
  Averted_per_dose = apd_str,
  Vacc_immune_end = final_v,
  row.names = NULL
)
print(summary_tbl)


# 6. Timing Sensitivity Sweep (Interactive Only) -----------------------------

# How does the cumulative burden change as the campaign start time slides
# across the full simulation horizon? The dose-response curve reveals the
# "deadline" for intervention effectiveness -- the start time beyond which
# the campaign's marginal benefit drops off sharply.
if (interactive()) {
  sweep_starts <- seq(0, 200, by = 20)
  sweep.df <- data.frame(
    .scenario.id = paste0("start_", sweep_starts),
    .at          = 0,
    vax.starts_1 = sweep_starts,
    vax.starts_2 = -1,
    vax.ends_1   = nsteps,
    vax.ends_2   = -1
  )
  sweep.list <- create_scenario_list(sweep.df)
  sweep_cum <- numeric(length(sweep.list))
  for (i in seq_along(sweep.list)) {
    cat("Sweep:", sweep_starts[i],
        "(", i, "/", length(sweep.list), ")\n")
    s <- netsim(est, use_scenario(param_base, sweep.list[[i]]),
                init, control)
    df <- as.data.frame(s)
    sweep_cum[i] <- mean(tapply(df$si.flow, df$sim, sum, na.rm = TRUE))
  }
  sweep_pct <- 100 * (cum_inf["none"] - sweep_cum) / cum_inf["none"]

  par(mfrow = c(1, 2), mar = c(3.5, 4, 2, 1), mgp = c(2.4, 1, 0))
  plot(sweep_starts, sweep_cum, type = "b", pch = 19, lwd = 2,
       col = "steelblue", xlab = "Campaign start timestep",
       ylab = "Cumulative infections",
       main = "Burden vs. Start Time")
  abline(h = cum_inf["none"], lty = 2, col = "gray50")
  plot(sweep_starts, sweep_pct, type = "b", pch = 19, lwd = 2,
       col = "seagreen", xlab = "Campaign start timestep",
       ylab = "% infections averted",
       main = "Benefit by Start Time", ylim = c(0, 100))

  print(data.frame(start = sweep_starts,
                   cum_inf = round(sweep_cum),
                   pct_averted = round(sweep_pct, 1)))
}
