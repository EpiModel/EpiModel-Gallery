
##
## SEIR with Contact Tracing for an Acute, Immunizing Infection
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##
## Author: Samuel M. Jenness (Emory University)
## Date: May 2026
##


# Attribute initializer ----------------------------------------------------

init_attrs <- function(dat, at) {
  # One-shot setup of the auxiliary node attributes used by the tracing
  # and quarantine machinery. Runs on the first module pass per sim. The
  # standard `active`, `status`, and `infTime` attributes are managed by
  # EpiModel itself; we add four bookkeeping attributes:
  #
  #   dx.time       -- integer time step a node was diagnosed (NA otherwise)
  #   dx.this.step  -- 0/1 flag set in the step of diagnosis (the tracing
  #                    trigger); reset to 0 in the next progress() call
  #   quar.until    -- last time step on which a node is still quarantined.
  #                    Transmission across an edge is reduced when either
  #                    endpoint has at <= quar.until.
  #   traced.count  -- cumulative number of times a node has been traced

  if (is.null(get_attr(dat, "dx.time", override.null.error = TRUE))) {
    active <- get_attr(dat, "active")
    n <- length(active)
    dat <- set_attr(dat, "dx.time", rep(NA_integer_, n))
    dat <- set_attr(dat, "dx.this.step", rep(0L, n))
    dat <- set_attr(dat, "quar.until", rep(NA_integer_, n))
    dat <- set_attr(dat, "traced.count", rep(0L, n))
  }

  # Seed infections enter as status == "i" (EpiModel's only initial
  # infectious value). Place them all in the presymptomatic substage so
  # they progress through the same Ip -> Is -> R pipeline as later
  # secondary cases.
  if (is.null(get_attr(dat, "inf.stage", override.null.error = TRUE))) {
    active <- get_attr(dat, "active")
    status <- get_attr(dat, "status")
    inf.stage <- rep(NA_character_, length(active))
    inf.stage[status == "i"] <- "ip"
    # Treat the initial seeds as living in the Ip substage. Promote them
    # to a canonical SEIR status by overwriting "i" with "ip"; downstream
    # modules and the prevalence module read inf.stage to count Ip / Is.
    status[status == "i"] <- "ip"
    dat <- set_attr(dat, "status", status)
    dat <- set_attr(dat, "inf.stage", inf.stage)
  }

  return(dat)
}


# Infection module ---------------------------------------------------------

infect <- function(dat, at) {
  # S -> E transmission along discordant edges, with one twist: the per-
  # edge act count is multiplied by quar.act.mult when either endpoint
  # has at <= quar.until. Both index isolation (after diagnosis) and
  # traced-contact quarantine route through that same act-rate cut.
  #
  # We do not use discord_edgelist() here, because in this model the
  # infectious side is identified by status %in% c("ip", "is") rather
  # than the literal status == "i". Walking the edgelist directly is the
  # same pattern as the multilayer rsv example.

  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  infTime <- get_attr(dat, "infTime")
  quar.until <- get_attr(dat, "quar.until")

  inf.prob <- get_param(dat, "inf.prob")
  act.rate <- get_param(dat, "act.rate")
  quar.act.mult <- get_param(dat, "quar.act.mult")

  # Pull the current cross-sectional edgelist (network 1 by default).
  el <- get_edgelist(dat, network = 1)

  nInf <- 0
  if (!is.null(el) && nrow(el) > 0) {
    head <- el[, 1]
    tail <- el[, 2]

    h_status <- status[head]
    t_status <- status[tail]
    h_active <- active[head] == 1
    t_active <- active[tail] == 1

    inf_states <- c("ip", "is")

    # Discordant edges in either orientation (one S, one Ip/Is).
    sus_in_h <- h_status == "s" & t_status %in% inf_states &
                h_active & t_active
    sus_in_t <- t_status == "s" & h_status %in% inf_states &
                h_active & t_active
    disc <- sus_in_h | sus_in_t

    if (any(disc, na.rm = TRUE)) {
      sus_ids <- ifelse(sus_in_h[disc], head[disc], tail[disc])
      inf_ids <- ifelse(sus_in_h[disc], tail[disc], head[disc])

      # Quarantine multiplier: if EITHER endpoint of the partnership is
      # currently quarantined (at <= quar.until), shrink act.rate. The
      # mechanism captures both directions of within-partnership avoidance.
      sus_quar <- !is.na(quar.until[sus_ids]) & at <= quar.until[sus_ids]
      inf_quar <- !is.na(quar.until[inf_ids]) & at <= quar.until[inf_ids]
      eff_acts <- act.rate * ifelse(sus_quar | inf_quar,
                                    quar.act.mult, 1)

      finalProb <- 1 - (1 - inf.prob)^eff_acts
      transmit <- rbinom(length(finalProb), 1, finalProb) == 1
      newInf <- unique(sus_ids[transmit])

      if (length(newInf) > 0) {
        nInf <- length(newInf)
        status[newInf] <- "e"
        infTime[newInf] <- at
        dat <- set_attr(dat, "status", status)
        dat <- set_attr(dat, "infTime", infTime)
      }
    }
  }

  dat <- set_epi(dat, "se.flow", at, nInf)
  return(dat)
}


# Progression module -------------------------------------------------------

progress <- function(dat, at) {
  # State machine:
  #   E   -> Ip  at ei.rate
  #   Ip  -> Is  at ips.rate
  #   Is  -> R   at isr.rate
  #   Is symptomatic -> diagnosed (status unchanged) at dx.rate.symp
  #
  # All transition candidates are taken from a snapshot of the state at
  # entry, so a node cannot cascade through multiple stages in a single
  # step. The diagnosis check also requires inf.stage == "is" (only
  # symptomatic cases are detectable through symptom-based testing); the
  # presymptomatic Ip stage is by construction undetectable by this
  # surveillance pathway.

  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  inf.stage <- get_attr(dat, "inf.stage")
  infTime <- get_attr(dat, "infTime")
  dx.time <- get_attr(dat, "dx.time")
  dx.this.step <- get_attr(dat, "dx.this.step")
  quar.until <- get_attr(dat, "quar.until")

  ei.rate <- get_param(dat, "ei.rate")
  ips.rate <- get_param(dat, "ips.rate")
  isr.rate <- get_param(dat, "isr.rate")
  dx.rate.symp <- get_param(dat, "dx.rate.symp")
  iso.duration <- get_param(dat, "iso.duration")

  # Entry-state snapshots
  status0 <- status
  stage0 <- inf.stage

  ## E -> Ip (require infTime < at so the latent period is at least one step)
  ids_e <- which(active == 1 & status0 == "e" &
                 !is.na(infTime) & infTime < at)
  n_ei <- 0
  if (length(ids_e) > 0) {
    hit <- which(rbinom(length(ids_e), 1, ei.rate) == 1)
    if (length(hit) > 0) {
      new_ip <- ids_e[hit]
      n_ei <- length(new_ip)
      status[new_ip] <- "ip"
      inf.stage[new_ip] <- "ip"
    }
  }

  ## Ip -> Is
  ids_ip <- which(active == 1 & status0 == "ip" &
                  !is.na(stage0) & stage0 == "ip")
  n_ips <- 0
  if (length(ids_ip) > 0) {
    hit <- which(rbinom(length(ids_ip), 1, ips.rate) == 1)
    if (length(hit) > 0) {
      new_is <- ids_ip[hit]
      n_ips <- length(new_is)
      status[new_is] <- "is"
      inf.stage[new_is] <- "is"
    }
  }

  ## Is -> R
  ids_is <- which(active == 1 & status0 == "is" &
                  !is.na(stage0) & stage0 == "is")
  n_ir <- 0
  if (length(ids_is) > 0) {
    hit <- which(rbinom(length(ids_is), 1, isr.rate) == 1)
    if (length(hit) > 0) {
      new_r <- ids_is[hit]
      n_ir <- length(new_r)
      status[new_r] <- "r"
      inf.stage[new_r] <- NA_character_
    }
  }

  # Reset the per-step diagnosis trigger from the previous step before
  # writing this step's diagnoses. The trace module reads dx.this.step
  # later in the same step, so the reset has to happen here.
  dx.this.step <- rep(0L, length(dx.this.step))

  ## Diagnosis of symptomatic, undiagnosed nodes.
  # Eligibility uses the SNAPSHOT inf.stage so a node that just
  # transitioned Ip -> Is in this same step is not diagnosed until the
  # following step (one full step of symptomatic transmission before
  # detection).
  ids_dx <- which(active == 1 & status0 == "is" &
                  !is.na(stage0) & stage0 == "is" &
                  is.na(dx.time))
  n_dx <- 0
  if (length(ids_dx) > 0) {
    hit <- which(rbinom(length(ids_dx), 1, dx.rate.symp) == 1)
    if (length(hit) > 0) {
      new_dx <- ids_dx[hit]
      n_dx <- length(new_dx)
      dx.time[new_dx] <- at
      dx.this.step[new_dx] <- 1L
      # Isolation: cut the index's act rate for the next iso.duration
      # steps. The infection module honours this through quar.until.
      quar.until[new_dx] <- pmax(quar.until[new_dx], at + iso.duration,
                                 na.rm = TRUE)
    }
  }

  dat <- set_attr(dat, "status", status)
  dat <- set_attr(dat, "inf.stage", inf.stage)
  dat <- set_attr(dat, "dx.time", dx.time)
  dat <- set_attr(dat, "dx.this.step", dx.this.step)
  dat <- set_attr(dat, "quar.until", quar.until)

  ## Summary statistics
  is_active <- active == 1
  dat <- set_epi(dat, "ei.flow", at, n_ei)
  dat <- set_epi(dat, "ips.flow", at, n_ips)
  dat <- set_epi(dat, "ir.flow", at, n_ir)
  dat <- set_epi(dat, "dx.flow", at, n_dx)
  dat <- set_epi(dat, "e.num", at, sum(is_active & status == "e"))
  dat <- set_epi(dat, "ip.num", at, sum(is_active & status == "ip"))
  dat <- set_epi(dat, "is.num", at, sum(is_active & status == "is"))
  dat <- set_epi(dat, "r.num", at, sum(is_active & status == "r"))
  # i.num combines Ip + Is so the standard EpiModel plots still work.
  dat <- set_epi(dat, "i.num", at,
                 sum(is_active & status %in% c("ip", "is")))

  return(dat)
}


# Prevalence module --------------------------------------------------------

prev <- function(dat, at) {
  # Replace EpiModel's default prevalence.net so that i.num reflects the
  # Ip + Is split used in this model. The default reads status == "i",
  # which is never true here (we use "ip" and "is" labels), so we
  # overwrite i.num with the corrected count and pass through s.num and
  # num the same way the default does.

  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  is_active <- active == 1

  dat <- set_epi(dat, "s.num", at, sum(is_active & status == "s"))
  dat <- set_epi(dat, "i.num", at,
                 sum(is_active & status %in% c("ip", "is")))
  dat <- set_epi(dat, "num", at, sum(is_active))

  return(dat)
}


# Contact tracing module ---------------------------------------------------

trace <- function(dat, at) {
  # Headline module. For each newly diagnosed index, traverse the
  # cumulative edgelist to find recent partners, Bernoulli-thin by the
  # reach probability, and quarantine the ones we reach.
  #
  # The simplest way to model trace.delay is to fire on indices whose
  # diagnosis happened `trace.delay` steps ago. So a node diagnosed on
  # step 30 triggers contact tracing on step 30 + trace.delay. This
  # avoids a per-node queue while preserving the diagnosis-to-reach gap.
  #
  # The headline EpiModel API showcased here:
  #
  #   1. control.net(cumulative.edgelist = TRUE) makes the simulation
  #      attach a running history of dissolved (and active) edges to
  #      `dat`. The history is truncated by control$truncate.el.cuml,
  #      which is set to trace.lookback in the model script.
  #
  #   2. get_partners(dat, index_pid, truncate = trace.lookback,
  #                   only.active.nodes = TRUE)
  #      returns one row per (index, partner) pair, with partner ids in
  #      the UNIQUE id space because a partner may have already left the
  #      simulation.
  #
  #   3. get_posit_ids(dat, unique_id) converts those unique ids back to
  #      positional ids so we can index attribute vectors. Closed
  #      populations make this round-trip seem cosmetic, but the
  #      get_partners contract guarantees unique ids so future
  #      vital-dynamics extensions of this code do not break silently.

  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  dx.time <- get_attr(dat, "dx.time")
  quar.until <- get_attr(dat, "quar.until")
  traced.count <- get_attr(dat, "traced.count")

  trace.reach.prob <- get_param(dat, "trace.reach.prob")
  trace.delay <- get_param(dat, "trace.delay")
  trace.lookback <- get_param(dat, "trace.lookback")
  quar.duration <- get_param(dat, "quar.duration")

  n_traced <- 0
  n_reached <- 0
  n_quar <- 0

  # Short-circuit when tracing is disabled (the no-tracing scenario).
  if (trace.reach.prob > 0) {

    # Indices whose diagnosis is `trace.delay` steps old this step.
    idsIndex <- which(active == 1 &
                      !is.na(dx.time) &
                      (at - dx.time) == trace.delay)

    if (length(idsIndex) > 0) {

      # --- The cumulative-edgelist round trip ---
      # get_partners takes positional ids in, returns unique ids out.
      part_df <- get_partners(dat, idsIndex,
                              truncate = trace.lookback,
                              only.active.nodes = TRUE)

      if (!is.null(part_df) && nrow(part_df) > 0) {
        # Translate partner unique ids back to positional ids. In a
        # closed population this is a no-op for currently-active nodes,
        # but the conversion is the principled way to do this lookup
        # because get_partners may include unique ids that do not have
        # a current positional id (e.g. departed partners in a vital-
        # dynamics extension).
        partner_pid <- get_posit_ids(dat, part_df$partner)
        partner_pid <- partner_pid[!is.na(partner_pid)]
        partner_pid <- unique(partner_pid)

        # Drop partners who are already diagnosed (their isolation is
        # already in effect via the index pathway) so we do not double-
        # count quarantines.
        partner_pid <- partner_pid[is.na(dx.time[partner_pid])]

        n_traced <- length(partner_pid)

        if (n_traced > 0) {
          # Bernoulli reach: each identified partner is contacted and
          # advised to quarantine with probability trace.reach.prob.
          reached_mask <- rbinom(n_traced, 1, trace.reach.prob) == 1
          reached_pid <- partner_pid[reached_mask]
          n_reached <- length(reached_pid)

          if (n_reached > 0) {
            # Apply quarantine: set quar.until forward `quar.duration`
            # steps. The infection module reads quar.until and shrinks
            # act.rate on edges whose endpoints are still quarantined.
            new_quar <- ifelse(is.na(quar.until[reached_pid]),
                               at + quar.duration,
                               pmax(quar.until[reached_pid],
                                    at + quar.duration))
            quar.until[reached_pid] <- new_quar
            traced.count[reached_pid] <- traced.count[reached_pid] + 1L
            n_quar <- n_reached

            dat <- set_attr(dat, "quar.until", quar.until)
            dat <- set_attr(dat, "traced.count", traced.count)
          }
        }
      }
    }
  }

  # Per-step counters that the analysis pipeline reads later.
  dat <- set_epi(dat, "trace.idx.flow", at, n_traced)
  dat <- set_epi(dat, "trace.reach.flow", at, n_reached)
  dat <- set_epi(dat, "trace.quar.flow", at, n_quar)

  # Current quarantine prevalence is useful for diagnostics.
  is_active <- active == 1
  dat <- set_epi(dat, "quar.num", at,
                 sum(is_active & !is.na(quar.until) & at <= quar.until))

  return(dat)
}
