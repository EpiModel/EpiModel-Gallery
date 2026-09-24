
##
## Partner Notification for an Endemic Bacterial STI
## EpiModel Gallery (https://github.com/EpiModel/EpiModel-Gallery)
##
## Author: Samuel M. Jenness (Emory University)
## Date: September 2026
##


# Attribute initializer ----------------------------------------------------

init_attrs <- function(dat, at) {
  # One-shot setup on the first module pass. EpiModel manages `status`,
  # `active`, `infTime`, `unique_id`, `entrTime`, and `exitTime`, and copies
  # the `sex` and `risk` vertex attributes from the network. Five attributes
  # are added:
  #
  #   symp       -- 1 if the current infection is symptomatic, 0 if not
  #   dx.time    -- step of the most recent diagnosis
  #   tx.time    -- step of the most recent cure of a case diagnosed by
  #                 screening or care-seeking; the clock for repeat
  #                 infection, cleared at the first repeat infection
  #   pn.pending -- 1 while a diagnosed case waits for partner services
  #   pr.visit   -- step on which a partner referred by an index comes in
  #                 for testing and treatment (patient referral only)
  #
  # A run resumed from a burn-in already has these attributes, so the block
  # is skipped.
  if (is.null(get_attr(dat, "dx.time", override.null.error = TRUE))) {
    status <- get_attr(dat, "status")
    sex <- get_attr(dat, "sex")
    n <- length(status)

    symp.prob <- ifelse(sex == "F", get_param(dat, "symp.prob.f"),
                        get_param(dat, "symp.prob.m"))
    symp <- rep(NA_integer_, n)
    seeds <- which(status == "i")
    symp[seeds] <- rbinom(length(seeds), 1, symp.prob[seeds])

    dat <- set_attr(dat, "symp", symp)
    for (nm in c("dx.time", "tx.time", "pr.visit")) {
      dat <- set_attr(dat, nm, rep(NA_integer_, n))
    }
    dat <- set_attr(dat, "pn.pending", rep(0L, n))
  }
  return(dat)
}


# Infection module ---------------------------------------------------------

infect <- function(dat, at) {
  # S -> I across the discordant edges of both partnership layers, with the
  # same weekly transmission probability in every partnership. A new
  # infection in someone cured within the last reinf.window weeks is a
  # repeat infection, and for each one the cumulative edgelist is asked
  # when the transmitting partnership began: on or before the cure means
  # the person was reinfected by a partner they already had when treated.
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  sex <- get_attr(dat, "sex")
  infTime <- get_attr(dat, "infTime")
  symp <- get_attr(dat, "symp")
  tx.time <- get_attr(dat, "tx.time")
  unique_id <- get_attr(dat, "unique_id")

  inf.prob <- get_param(dat, "inf.prob")
  symp.prob.f <- get_param(dat, "symp.prob.f")
  symp.prob.m <- get_param(dat, "symp.prob.m")
  reinf.window <- get_param(dat, "reinf.window")

  # Discordant pairs from each layer, labeled with the layer
  del <- NULL
  for (k in 1:2) {
    del <- rbind(del, discord_edgelist(dat, at, network = k,
                                       include.network = TRUE))
  }

  n_inf <- 0
  n_inf_f <- 0
  n_inf_layer <- c(0, 0)
  n_reinf <- 0
  n_reinf_prior <- 0
  if (!is.null(del)) {
    del <- del[rbinom(nrow(del), 1, inf.prob) == 1, , drop = FALSE]
    # A susceptible exposed by several partners in one week is infected
    # once, by a randomly chosen one of them.
    del <- del[sample.int(nrow(del)), , drop = FALSE]
    del <- del[!duplicated(del$sus), , drop = FALSE]
    n_inf <- nrow(del)

    if (n_inf > 0) {
      new <- del$sus
      n_inf_f <- sum(sex[new] == "F")
      n_inf_layer <- tabulate(del$network, nbins = 2)
      status[new] <- "i"
      infTime[new] <- at
      symp[new] <- rbinom(n_inf, 1, ifelse(sex[new] == "F", symp.prob.f,
                                           symp.prob.m))
      dat <- set_attr(dat, "status", status)
      dat <- set_attr(dat, "infTime", infTime)
      dat <- set_attr(dat, "symp", symp)

      # Repeat infections, and whether the source partnership predates the
      # cure. truncate = 0 returns ongoing partnerships only.
      rep_inf <- !is.na(tx.time[new]) & at - tx.time[new] <= reinf.window
      n_reinf <- sum(rep_inf)
      if (n_reinf > 0) {
        src <- data.frame(index = unique_id[new[rep_inf]],
                          partner = unique_id[del$inf[rep_inf]],
                          network = del$network[rep_inf],
                          tx.time = tx.time[new[rep_inf]])
        part_df <- get_partners(dat, new[rep_inf], truncate = 0)
        src <- merge(src, part_df, by = c("index", "partner", "network"))
        n_reinf_prior <- sum(src$start <= src$tx.time)
        # Each cure is followed until its first repeat infection only
        tx.time[new[rep_inf]] <- NA
        dat <- set_attr(dat, "tx.time", tx.time)
      }
    }
  }

  dat <- set_epi(dat, "si.flow", at, n_inf)
  dat <- set_epi(dat, "si.flow.f", at, n_inf_f)
  dat <- set_epi(dat, "si.flow.main", at, n_inf_layer[1])
  dat <- set_epi(dat, "reinf.flow", at, n_reinf)
  dat <- set_epi(dat, "reinf.prior.flow", at, n_reinf_prior)
  return(dat)
}


# Testing and treatment module ---------------------------------------------

test_treat <- function(dat, at) {
  # Symptomatic infections come to care at symp.test.rate per week, and
  # asymptomatic infections are found by screening at a weekly rate that
  # depends on sex. A positive test is a diagnosis: the case is treated,
  # cured with probability tx.prob, and queued for partner services. Tests
  # of uninfected people change nothing in the model and are not simulated.
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  sex <- get_attr(dat, "sex")
  symp <- get_attr(dat, "symp")
  dx.time <- get_attr(dat, "dx.time")
  tx.time <- get_attr(dat, "tx.time")
  pn.pending <- get_attr(dat, "pn.pending")

  screen.rate.f <- get_param(dat, "screen.rate.f")
  screen.rate.m <- get_param(dat, "screen.rate.m")
  symp.test.rate <- get_param(dat, "symp.test.rate")
  tx.prob <- get_param(dat, "tx.prob")

  ids <- which(active == 1 & status == "i")
  rate <- ifelse(symp[ids] == 1, symp.test.rate,
                 ifelse(sex[ids] == "F", screen.rate.f, screen.rate.m))
  dx <- ids[rbinom(length(ids), 1, rate) == 1]
  cured <- dx[rbinom(length(dx), 1, tx.prob) == 1]

  dx.time[dx] <- at
  pn.pending[dx] <- 1L
  status[cured] <- "s"
  tx.time[cured] <- at

  dat <- set_attr(dat, "status", status)
  dat <- set_attr(dat, "dx.time", dx.time)
  dat <- set_attr(dat, "tx.time", tx.time)
  dat <- set_attr(dat, "pn.pending", pn.pending)

  dat <- set_epi(dat, "dx.flow", at, length(dx))
  dat <- set_epi(dat, "dx.flow.f", at, sum(sex[dx] == "F"))
  dat <- set_epi(dat, "tx.flow", at, length(cured))
  return(dat)
}


# Partner notification module ----------------------------------------------

notify <- function(dat, at) {
  # Partner services for every diagnosed case in the queue. The partners of
  # each index in the last pn.lookback weeks come from the cumulative
  # edgelist of both layers, and each is reached with a probability that
  # depends on whether the partnership is still going and, if it is, on
  # its layer. Treatment is presumptive under both arms. Under expedited
  # partner therapy (EPT) a reached partner takes medication delivered by
  # the index the same week and is never tested. Under patient referral
  # (PR) a reached partner comes in the following week, is tested and
  # treated, and if infected is diagnosed and queued, so that their own
  # partners are notified in turn.
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  sex <- get_attr(dat, "sex")
  dx.time <- get_attr(dat, "dx.time")
  pn.pending <- get_attr(dat, "pn.pending")
  pr.visit <- get_attr(dat, "pr.visit")

  pn.arm <- get_param(dat, "pn.arm")
  pn.lookback <- get_param(dat, "pn.lookback")
  reach <- c(main = get_param(dat, "reach.main"),
             cas = get_param(dat, "reach.cas"),
             ended = get_param(dat, "reach.ended"))
  tx.prob <- get_param(dat, "tx.prob")

  cats <- names(reach)
  n_found <- n_reach <- n_reach_inf <- setNames(c(0, 0, 0), cats)
  n_index <- 0
  n_cured <- 0
  n_visit <- 0
  n_part_dx <- 0
  n_part_dx_f <- 0

  # 1. PR: partners referred last week come in, are tested, and are treated.
  #    This runs under every arm, so partners booked in the last week of the
  #    burn-in still come in during the first week of any scenario.
  visit <- which(active == 1 & pr.visit == at)
  if (length(visit) > 0) {
    n_visit <- length(visit)
    inf <- visit[status[visit] == "i"]
    cured <- inf[rbinom(length(inf), 1, tx.prob) == 1]
    status[cured] <- "s"
    dx.time[inf] <- at
    pn.pending[inf] <- 1L
    n_cured <- length(cured)
    n_part_dx <- length(inf)
    n_part_dx_f <- sum(sex[inf] == "F")
  }

  # 2. Indices: everyone in the queue, including partners diagnosed in step 1
  idsIndex <- which(active == 1 & pn.pending == 1)
  pn.pending[idsIndex] <- 0L

  if (pn.arm != "none" && length(idsIndex) > 0) {
    n_index <- length(idsIndex)

    # 3. Partners in the lookback window, from both layers
    part_df <- get_partners(dat, idsIndex, truncate = pn.lookback,
                            only.active.nodes = TRUE)

    if (nrow(part_df) > 0) {
      # 4. Partner unique ids back to positional ids; one row per partner,
      #    keeping an ongoing partnership over an ended one. Every partner in
      #    the window is counted as found, but partners diagnosed this week
      #    or last are already in care and are not contacted, which also
      #    stops a partner diagnosed through PR from notifying the index who
      #    referred them.
      part_df$pid <- get_posit_ids(dat, part_df$partner)
      part_df <- part_df[order(!is.na(part_df$stop)), ]
      part_df <- part_df[!duplicated(part_df$pid), ]
      type <- ifelse(!is.na(part_df$stop), "ended",
                     ifelse(part_df$network == 1, "main", "cas"))
      n_found[] <- table(factor(type, levels = cats))
      in_care <- !is.na(dx.time[part_df$pid]) & dx.time[part_df$pid] >= at - 1
      part_df <- part_df[!in_care, ]
      type <- type[!in_care]

      # 5. Reach each partner with the probability for its partnership type
      reached <- rbinom(nrow(part_df), 1, reach[type]) == 1
      pid <- part_df$pid[reached]
      infected <- status[pid] == "i"
      n_reach[] <- table(factor(type[reached], levels = cats))
      n_reach_inf[] <- table(factor(type[reached][infected], levels = cats))

      # 6. EPT: treated now. PR: booked to come in next week.
      if (pn.arm == "EPT") {
        inf <- pid[infected]
        cured <- inf[rbinom(length(inf), 1, tx.prob) == 1]
        status[cured] <- "s"
        n_cured <- n_cured + length(cured)
      } else if (pn.arm == "PR") {
        pr.visit[pid] <- at + 1
      }
    }
  }

  dat <- set_attr(dat, "status", status)
  dat <- set_attr(dat, "dx.time", dx.time)
  dat <- set_attr(dat, "pn.pending", pn.pending)
  dat <- set_attr(dat, "pr.visit", pr.visit)

  dat <- set_epi(dat, "pn.index.flow", at, n_index)
  for (k in cats) {
    dat <- set_epi(dat, paste0("pn.found.", k, ".flow"), at, n_found[[k]])
    dat <- set_epi(dat, paste0("pn.reach.", k, ".flow"), at, n_reach[[k]])
    dat <- set_epi(dat, paste0("pn.inf.", k, ".flow"), at, n_reach_inf[[k]])
  }
  dat <- set_epi(dat, "pn.visit.flow", at, n_visit)
  dat <- set_epi(dat, "pn.cured.flow", at, n_cured)
  dat <- set_epi(dat, "pn.dx.flow", at, n_part_dx)
  dat <- set_epi(dat, "pn.dx.flow.f", at, n_part_dx_f)
  return(dat)
}


# Natural clearance module -------------------------------------------------

clear <- function(dat, at) {
  # Untreated infections clear on their own at rec.rate per week.
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  rec.rate <- get_param(dat, "rec.rate")

  ids <- which(active == 1 & status == "i")
  rec <- ids[rbinom(length(ids), 1, rec.rate) == 1]
  status[rec] <- "s"

  dat <- set_attr(dat, "status", status)
  dat <- set_epi(dat, "is.flow", at, length(rec))
  return(dat)
}


# Departure and arrival modules --------------------------------------------

depart <- function(dat, at) {
  # People leave at departure.rate per week, which stands for aging out of
  # the 15 to 24 age band. With tergmLite, departed nodes are deleted at the
  # end of the step, and the positional ids of everyone after them shift.
  active <- get_attr(dat, "active")
  exitTime <- get_attr(dat, "exitTime")
  departure.rate <- get_param(dat, "departure.rate")

  ids <- which(active == 1)
  dep <- ids[rbinom(length(ids), 1, departure.rate) == 1]
  active[dep] <- 0L
  exitTime[dep] <- at

  dat <- set_attr(dat, "active", active)
  dat <- set_attr(dat, "exitTime", exitTime)
  dat <- set_epi(dat, "d.flow", at, length(dep))
  return(dat)
}

arrive <- function(dat, at) {
  # New susceptible people arrive at a rate that balances departures. The
  # custom attributes are extended here; `sex` and `risk` are vertex
  # attributes of the network, so EpiModel assigns them to arrivals itself,
  # following the attr.rules setting in control.net().
  n <- sum(get_attr(dat, "active") == 1)
  arrival.rate <- get_param(dat, "arrival.rate")
  n_arr <- rpois(1, n * arrival.rate)

  if (n_arr > 0) {
    dat <- append_core_attr(dat, at, n_arr)
    dat <- append_attr(dat, "status", "s", n_arr)
    dat <- append_attr(dat, "infTime", NA_integer_, n_arr)
    for (nm in c("symp", "dx.time", "tx.time", "pr.visit")) {
      dat <- append_attr(dat, nm, NA_integer_, n_arr)
    }
    dat <- append_attr(dat, "pn.pending", 0L, n_arr)
  }
  dat <- set_epi(dat, "a.flow", at, n_arr)
  return(dat)
}
