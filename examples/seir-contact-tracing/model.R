
##
## SEIR with Contact Tracing for an Acute, Immunizing Infection
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##
## Author: Samuel M. Jenness (Emory University)
## Date: May 2026
##

suppressMessages(library(EpiModel))

# Standard Gallery unit test lines
rm(list = ls())
eval(parse(text = print(commandArgs(TRUE)[1])))

# Two run modes:
#   interactive(): full network, 5 sims, 200-step horizon, all 4 scenarios.
#   non-interactive (CI): small network, 2 sims, 80 steps. CI mode is
#     calibrated to complete in well under a minute.
if (interactive()) {
  nsims <- 5
  ncores <- 5
  nsteps <- 200
  n <- 500
} else {
  nsims <- 2
  ncores <- 2
  nsteps <- 80
  n <- 200
}


# 1. Network Model Estimation ----------------------------------------------

# A single-layer dynamic contact network of n nodes. Mean degree 3
# (target.stats = 1.5 * n edges) with short partnership duration so the
# tracing lookback window has meaningful turnover.
#
# Each timestep is interpreted as 1 day for parameter readability.
nw <- network_initialize(n)
formation <- ~edges
target.stats <- c(round(1.5 * n))

coef.diss <- dissolution_coefs(~offset(edges), duration = 10)
est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)

dx <- netdx(est, nsims = nsims, ncores = ncores, nsteps = nsteps,
            nwstats.formula = ~edges + degree(0:4),
            verbose = FALSE)
print(dx)
if (interactive()) plot(dx)


# 2. Parameters and Modules -----------------------------------------------

source("examples/seir-contact-tracing/module-fx.R")

# Disease parameters reflect a COVID-like acute, immunizing infection
# with a presymptomatic infectious window. Each timestep = 1 day.
#
#   inf.prob     = 0.05   per-act transmission probability
#   act.rate     = 2      acts per partnership per day
#   ei.rate      = 1/3    mean 3 days latent
#   ips.rate     = 1/2    mean 2 days presymptomatic
#   isr.rate     = 1/6    mean 6 days symptomatic
#   dx.rate.symp = 0.5    per-step probability an undiagnosed symptomatic
#                         case is detected
#   iso.duration = 10     days an index isolates after diagnosis
#
# Tracing parameters (the scenario knobs):
#   trace.reach.prob -- per-partner probability that a contact is
#                       successfully reached and advised to quarantine
#   trace.delay      -- days from diagnosis to contact reach
#   trace.lookback   -- days of partner history traced
#   quar.duration    -- days a reached contact stays in quarantine
#   quar.act.mult    -- act.rate multiplier when either endpoint is
#                       currently quarantined (0 = perfect isolation)

make_param <- function(trace.reach.prob = 0.0,
                       trace.delay = 0,
                       trace.lookback = 3,
                       quar.duration = 10,
                       quar.act.mult = 0.1) {
  param.net(
    inf.prob = 0.05,
    act.rate = 2,
    ei.rate = 1 / 3,
    ips.rate = 1 / 2,
    isr.rate = 1 / 6,
    dx.rate.symp = 0.5,
    iso.duration = 10,
    trace.reach.prob = trace.reach.prob,
    trace.delay = trace.delay,
    trace.lookback = trace.lookback,
    quar.duration = quar.duration,
    quar.act.mult = quar.act.mult
  )
}

init <- init.net(i.num = 10)

# control.net() switches that activate the cumulative-edgelist machinery:
#   cumulative.edgelist     = TRUE   build the running edge history
#   truncate.el.cuml        = N      drop edges older than N steps
#   save.cumulative.edgelist = TRUE  attach the final history to the sim
#
# truncate.el.cuml is the destructive trim: edges older than the lookback
# window are discarded from `dat` to keep memory bounded. We pick the
# lookback from the parameter set so the same control object works for
# every scenario.
trace.lookback.default <- 3

control <- control.net(
  type = NULL,
  nsims = nsims,
  ncores = ncores,
  nsteps = nsteps,
  cumulative.edgelist = TRUE,
  truncate.el.cuml = trace.lookback.default,
  save.cumulative.edgelist = TRUE,
  initialize.FUN = initialize.net,
  initAttr.FUN = init_attrs,
  infection.FUN = infect,
  progress.FUN = progress,
  trace.FUN = trace,
  prevalence.FUN = prev,
  verbose = FALSE
)


# 3. Scenarios -------------------------------------------------------------

# Four scenarios on the same fitted network, same disease parameters,
# same 10 seed infections. They differ only in the tracing configuration.
#
#   none      no tracing, the counterfactual
#   fast_high fast (1 day) trace, high (80%) reach
#   slow_high slow (4 day) trace, same 80% reach
#   fast_low  fast (1 day) trace, low (30%) reach

scenarios <- list(
  none      = list(reach = 0.0, delay = 0),
  fast_high = list(reach = 0.8, delay = 1),
  slow_high = list(reach = 0.8, delay = 4),
  fast_low  = list(reach = 0.3, delay = 1)
)

sims <- list()
for (s in names(scenarios)) {
  cat("\n--- Running scenario:", s, "---\n")
  cfg <- scenarios[[s]]
  p <- make_param(trace.reach.prob = cfg$reach,
                  trace.delay = cfg$delay)
  sims[[s]] <- netsim(est, p, init, control)
  print(sims[[s]])
}


# 4. Analysis --------------------------------------------------------------

labels <- c(none = "No tracing",
            fast_high = "Fast + high (delay 1, 80%)",
            slow_high = "Slow + high (delay 4, 80%)",
            fast_low  = "Fast + low (delay 1, 30%)")
cols <- c(none = "gray40", fast_high = "seagreen",
          slow_high = "firebrick", fast_low = "steelblue")

# Per-scenario summary metrics from the sim time series. The model
# tracks Ip and Is separately, so total infectious prevalence is
# computed as ip.num + is.num. EpiModel's built-in prevalence module
# overwrites i.num after our progress module runs, so we recompute it
# here from the parallel ip.num and is.num counts.
summarise_sim <- function(sim, npop) {
  df <- as.data.frame(sim)
  df$infectious <- df$ip.num + df$is.num
  df$infectious[is.na(df$infectious)] <- 0
  # Cumulative incidence (mean across sims) at the final step.
  cum_inc <- mean(tapply(df$se.flow, df$sim, sum, na.rm = TRUE))
  # Peak prevalence and peak day, computed per sim then averaged.
  peak_per_sim <- tapply(df$infectious, df$sim, max, na.rm = TRUE)
  peak_day_per_sim <- tapply(seq_len(nrow(df)), df$sim, function(idx) {
    df$time[idx][which.max(df$infectious[idx])]
  })
  peak_prev <- mean(peak_per_sim, na.rm = TRUE) / npop
  peak_day <- mean(peak_day_per_sim, na.rm = TRUE)
  total_dx <- mean(tapply(df$dx.flow, df$sim, sum, na.rm = TRUE))
  total_trace_idx <- mean(tapply(df$trace.idx.flow, df$sim, sum,
                                 na.rm = TRUE))
  total_reach <- mean(tapply(df$trace.reach.flow, df$sim, sum,
                             na.rm = TRUE))
  total_quar <- mean(tapply(df$trace.quar.flow, df$sim, sum,
                            na.rm = TRUE))
  data.frame(cum_inc = cum_inc,
             peak_prev = peak_prev,
             peak_day = peak_day,
             total_dx = total_dx,
             total_trace_idx = total_trace_idx,
             total_reach = total_reach,
             total_quar = total_quar)
}

summary_tbl <- do.call(rbind, lapply(sims, summarise_sim, npop = n))
summary_tbl$scenario <- labels[rownames(summary_tbl)]
summary_tbl$reach_per_idx <- ifelse(summary_tbl$total_dx > 0,
                                    summary_tbl$total_reach /
                                      summary_tbl$total_dx,
                                    NA)
summary_tbl <- summary_tbl[, c("scenario", "cum_inc", "peak_prev",
                               "peak_day", "total_dx",
                               "total_trace_idx", "total_reach",
                               "total_quar", "reach_per_idx")]
rownames(summary_tbl) <- NULL
print(summary_tbl)


## --- Plot 1: Cumulative incidence over time, all scenarios --------------
# Cumulative SE flow per sim, averaged across sims at each step. Captures
# the "final size" signal that policy makers actually care about. The
# flow columns are NA at t = 1 (no module has run yet) so we treat that
# as zero before cumulating.
par(mfrow = c(1, 1), mar = c(3.5, 4, 2.5, 1), mgp = c(2.4, 1, 0))
ymax <- 0
cum_inc_list <- list()
for (s in names(sims)) {
  df <- as.data.frame(sims[[s]])
  df$se.flow[is.na(df$se.flow)] <- 0
  cum_by_sim <- by(df, df$sim, function(d) cumsum(d$se.flow))
  M <- do.call(cbind, lapply(cum_by_sim, function(v) {
    out <- rep(NA_real_, nsteps)
    out[seq_along(v)] <- v
    out
  }))
  cum_mean <- rowMeans(M, na.rm = TRUE)
  cum_inc_list[[s]] <- cum_mean
  if (any(is.finite(cum_mean))) {
    ymax <- max(ymax, max(cum_mean, na.rm = TRUE))
  }
}
plot(seq_len(nsteps), cum_inc_list[["none"]], type = "n",
     ylim = c(0, max(ymax, 1) * 1.05),
     xlab = "Time step (days)",
     ylab = "Cumulative new infections",
     main = "Cumulative Incidence by Tracing Configuration")
for (s in names(sims)) {
  lines(seq_len(nsteps), cum_inc_list[[s]], lwd = 2, col = cols[s])
}
legend("topleft", legend = labels, col = cols, lwd = 2, bty = "n",
       cex = 0.85)


## --- Plot 2: Daily new infections, all scenarios -----------------------
# Peak suppression is visually clearer on the incidence curve than on the
# cumulative curve. The fast + high scenario should clip the peak hardest.
# Daily means are noisy at this sim count, so we overlay a 7-day centered
# moving average and draw the raw curves at low alpha for context.
par(mfrow = c(1, 1), mar = c(3.5, 4, 2.5, 1), mgp = c(2.4, 1, 0))
smooth_ma <- function(x, k = 7) {
  if (length(x) < k) return(x)
  out <- as.numeric(stats::filter(x, rep(1 / k, k), sides = 2))
  names(out) <- names(x)
  out
}
ymax2 <- 0
inc_list <- list()
inc_smooth <- list()
for (s in names(sims)) {
  df <- as.data.frame(sims[[s]])
  inc_mean <- tapply(df$se.flow, df$time, mean, na.rm = TRUE)
  inc_list[[s]] <- inc_mean
  inc_smooth[[s]] <- smooth_ma(inc_mean, k = 7)
  ymax2 <- max(ymax2, max(inc_smooth[[s]], na.rm = TRUE))
}
plot(as.numeric(names(inc_list[["none"]])), inc_list[["none"]],
     type = "n", ylim = c(0, ymax2 * 1.2),
     xlab = "Time step (days)",
     ylab = "New infections (daily mean)",
     main = "Daily New Infections by Tracing Configuration")
for (s in names(sims)) {
  lines(as.numeric(names(inc_list[[s]])), inc_list[[s]],
        lwd = 1, col = adjustcolor(cols[s], alpha.f = 0.25))
}
for (s in names(sims)) {
  lines(as.numeric(names(inc_smooth[[s]])), inc_smooth[[s]],
        lwd = 2.5, col = cols[s])
}
legend("topright", legend = labels, col = cols, lwd = 2.5, bty = "n",
       cex = 0.85)


## --- Plot 3: Tracing cascade ----------------------------------------------
# Per scenario, three bars: total diagnoses, total partners reached,
# total partner-quarantines initiated. The ratio annotations help readers
# see that "partners-per-index" and "quarantines-per-partner" are
# independent properties of the tracing program.
cascade <- t(as.matrix(summary_tbl[, c("total_dx", "total_reach",
                                       "total_quar")]))
cascade_labels <- c(none      = "No tracing",
                    fast_high = "Fast + high\n(d=1, 80%)",
                    slow_high = "Slow + high\n(d=4, 80%)",
                    fast_low  = "Fast + low\n(d=1, 30%)")
colnames(cascade) <- cascade_labels[names(sims)]
rownames(cascade) <- c("Diagnoses", "Partners reached",
                       "Quarantines initiated")

par(mfrow = c(1, 1), mar = c(4.5, 4, 3, 1), mgp = c(2.5, 1, 0))
bp <- barplot(cascade, beside = TRUE,
              col = c("#34495e", "#f39c12", "#e74c3c"),
              names.arg = colnames(cascade),
              las = 1, cex.names = 0.85,
              ylab = "Count (mean per sim)",
              main = "Tracing Cascade by Scenario",
              ylim = c(0, max(cascade, na.rm = TRUE) * 1.25))
legend("topright", legend = rownames(cascade),
       fill = c("#34495e", "#f39c12", "#e74c3c"), bty = "n",
       cex = 0.85)

# Annotate per-scenario ratios beneath the group of bars.
ratio_str <- ifelse(summary_tbl$total_dx > 0,
                    sprintf("%.1f reached/idx",
                            summary_tbl$reach_per_idx),
                    "")
mtext(side = 3, at = colMeans(bp), line = -0.5,
      text = ratio_str, cex = 0.8, font = 3, col = "gray30")


## --- Summary table -------------------------------------------------------
print_tbl <- data.frame(
  Scenario = summary_tbl$scenario,
  Cum_inc = round(summary_tbl$cum_inc),
  Peak_prev = round(summary_tbl$peak_prev, 3),
  Peak_day = round(summary_tbl$peak_day),
  Total_dx = round(summary_tbl$total_dx),
  Total_reach = round(summary_tbl$total_reach),
  Reach_per_idx = round(summary_tbl$reach_per_idx, 2)
)
print(print_tbl)
