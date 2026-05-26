
##
## Partner Notification for an Endemic STI (SIS)
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##
## Author: Samuel M. Jenness (Emory University)
## Date: May 2026
##


# Screening Module ---------------------------------------------------------

screen <- function(dat, at) {
  # Routine screening for the asymptomatic infected population. Every
  # infected, active node is screened (Bernoulli) at rate screen.rate per
  # timestep. A positive test sets dx.this.step = 1 (the trigger for both
  # treatment and partner notification later in the same timestep) and
  # stamps the dx.time. The diag.status flag is sticky across the run: it
  # remains 1 even after clearance, marking the node as previously
  # diagnosed for the cascade analysis.
  #
  # Initialization happens here on the first call (at == 2). EpiModel
  # reserves at == 1 for init.net setup.

  ## Attributes ##
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")

  if (at == 2) {
    n <- length(active)
    dat <- set_attr(dat, "diag.status",     rep(0, n))
    dat <- set_attr(dat, "dx.time",         rep(NA_integer_, n))
    dat <- set_attr(dat, "dx.this.step",    rep(0, n))
    dat <- set_attr(dat, "tx.this.step",    rep(0, n))
    dat <- set_attr(dat, "pn.notified",     rep(NA_integer_, n))
    dat <- set_attr(dat, "infections",      ifelse(status == "i", 1, 0))
  }

  diag.status <- get_attr(dat, "diag.status")
  dx.time     <- get_attr(dat, "dx.time")

  # Reset the per-step flags before this step's events fire.
  dx.this.step <- rep(0L, length(active))
  tx.this.step <- rep(0L, length(active))

  ## Parameters ##
  screen.rate <- get_param(dat, "screen.rate")
  pn.start    <- get_param(dat, "pn.start")

  ## Screening of infected actives (any prior diagnosis is fine, re-test
  ## is allowed). Susceptibles are not screened in this simple model:
  ## the focus is the partner-notification pathway, not the diagnostic
  ## cascade for negatives.
  nDxNew <- 0
  if (screen.rate > 0) {
    idsElig <- which(active == 1 & status == "i")
    nElig <- length(idsElig)
    if (nElig > 0) {
      vec <- which(rbinom(nElig, 1, screen.rate) == 1)
      if (length(vec) > 0) {
        idsDx <- idsElig[vec]
        nDxNew <- length(idsDx)
        dx.this.step[idsDx] <- 1
        diag.status[idsDx]  <- 1
        dx.time[idsDx]      <- at
      }
    }
  }

  ## Write attributes ##
  dat <- set_attr(dat, "diag.status",  diag.status)
  dat <- set_attr(dat, "dx.time",      dx.time)
  dat <- set_attr(dat, "dx.this.step", dx.this.step)
  dat <- set_attr(dat, "tx.this.step", tx.this.step)

  ## Summary statistics ##
  # Indices are the fresh-positive screening cases this step (the PN
  # triggers). PN is not yet active during burn-in (at < pn.start).
  dat <- set_epi(dat, "n.indices",   at, nDxNew)
  dat <- set_epi(dat, "pn.on",       at, as.numeric(at >= pn.start))
  dat <- set_epi(dat, "n.diag.ever", at,
                 sum(active == 1 & diag.status == 1, na.rm = TRUE))

  return(dat)
}


# Partner Notification Module ----------------------------------------------

partner_services <- function(dat, at) {
  # Cumulative-edgelist partner notification. The recipe is short and
  # reusable for any partner-tracing intervention in EpiModel. Three
  # components do the work:
  #
  #   1. control.net() options:
  #        cumulative.edgelist      = TRUE   (must be on; lets the engine
  #                                           collect past partnerships)
  #        truncate.el.cuml         = pn.lookback  (destructive: edges
  #                                           older than this are
  #                                           permanently dropped from
  #                                           the stored history)
  #        save.cumulative.edgelist = TRUE   (attach the edgelist to the
  #                                           returned sim object so it
  #                                           can be inspected after)
  #
  #   2. Inside this module, call get_partners() on the freshly diagnosed
  #      indices with the same lookback window:
  #
  #        idsIndex <- which(active == 1 & dx.this.step == 1)
  #        part_df  <- get_partners(dat, idsIndex,
  #                                 truncate          = pn.lookback,
  #                                 only.active.nodes = TRUE)
  #
  #   3. ID space conversion. get_partners() returns UNIQUE IDs in the
  #      `partner` column, because a past partner may have already
  #      departed the population and therefore no longer has a positional
  #      ID. The values must be translated back to positional IDs with
  #      get_posit_ids() before indexing into any per-node vector
  #      (status, active, etc.). only.active.nodes = TRUE keeps only
  #      partners still active, which guarantees the round-trip succeeds
  #      (no NAs from get_posit_ids), but the conversion is still required.
  #
  # After partners are converted to positional IDs the module draws a
  # Bernoulli at pn.trace.prob ("successfully notified") and stamps
  # pn.notified = at on the reached partners. The downstream treat module
  # reads pn.notified to determine how to handle each notified partner
  # under the active PN arm.

  ## Attributes ##
  active       <- get_attr(dat, "active")
  status       <- get_attr(dat, "status")
  dx.this.step <- get_attr(dat, "dx.this.step")
  pn.notified  <- get_attr(dat, "pn.notified")

  ## Parameters ##
  pn.lookback   <- get_param(dat, "pn.lookback")
  pn.trace.prob <- get_param(dat, "pn.trace.prob")
  pn.start      <- get_param(dat, "pn.start")

  ## Default counters (in case PN is off this step) ##
  n.partners.elig <- 0
  n.partners.reached <- 0

  if (at >= pn.start && pn.trace.prob > 0) {

    ## Step 1: identify fresh-positive indices this step ##
    idsIndex <- which(active == 1 & dx.this.step == 1)

    if (length(idsIndex) > 0) {

      ## Step 2: pull partners from the cumulative edgelist within the
      ## lookback window, restricted to partners still in the network.
      part_df <- get_partners(
        dat, idsIndex,
        truncate          = pn.lookback,
        only.active.nodes = TRUE
      )

      if (!is.null(part_df) && nrow(part_df) > 0) {

        ## Step 3: ID conversion. get_partners() puts unique IDs in
        ## `partner`. We need positional IDs to index into per-node
        ## vectors. Because we set only.active.nodes = TRUE, every
        ## unique ID maps to a current positional ID; we still drop
        ## any NAs defensively (a partner could have just departed
        ## between the cumulative-edgelist update and this call).
        partner_pid <- get_posit_ids(dat, part_df$partner)
        partner_pid <- partner_pid[!is.na(partner_pid)]

        ## Dedupe in case the same partner appears across multiple
        ## indices or multiple edges within the window.
        partner_pid <- unique(partner_pid)

        ## Drop any partners who happen to also be fresh indices this
        ## step. Their fate is handled by the index-treatment path in
        ## treat(); double-counting them as notified partners would
        ## inflate the cascade stats.
        partner_pid <- setdiff(partner_pid, idsIndex)

        n.partners.elig <- length(partner_pid)

        ## Bernoulli trace ("successfully notified") ##
        if (n.partners.elig > 0) {
          reached <- partner_pid[
            rbinom(n.partners.elig, 1, pn.trace.prob) == 1
          ]
          n.partners.reached <- length(reached)
          if (n.partners.reached > 0) {
            pn.notified[reached] <- at
          }
        }
      }
    }
  }

  dat <- set_attr(dat, "pn.notified", pn.notified)

  dat <- set_epi(dat, "n.partners.elig",    at, n.partners.elig)
  dat <- set_epi(dat, "n.partners.reached", at, n.partners.reached)

  return(dat)
}


# Treatment Module ---------------------------------------------------------

treat <- function(dat, at) {
  # Treatment cascade for both indices and notified partners.
  #
  #   Indices (dx.this.step == 1): always offered treatment. Success
  #     probability tx.efficacy. If successful and the index is still
  #     infected, clear them (status = "s"). Tracked as "index" treatment.
  #
  #   Notified partners (pn.notified == at) under three PN arms:
  #     "none": no partner services. (partner_services emits no
  #             notifications when pn.trace.prob == 0, so this branch
  #             will not fire; the arm is here for symmetry.)
  #     "PR"  : Patient Referral. The partner returns for clinical
  #             evaluation. They are tested first with sensitivity
  #             pn.test.prob; if positive they are then treated at
  #             tx.efficacy. Uninfected partners are not treated.
  #     "EPT" : Expedited Partner Therapy. The partner is dispensed
  #             medication directly without testing. Treat both infected
  #             and uninfected partners at ept.efficacy; uninfected
  #             treatment is recorded as a "wasted dose" in the cascade
  #             stats (no epidemic effect, but real cost).
  #
  # Treatment events flip status from "i" to "s" stochastically. Recovery
  # via natural clearance happens in a separate recovery module.

  ## Attributes ##
  active       <- get_attr(dat, "active")
  status       <- get_attr(dat, "status")
  dx.this.step <- get_attr(dat, "dx.this.step")
  tx.this.step <- get_attr(dat, "tx.this.step")
  pn.notified  <- get_attr(dat, "pn.notified")
  diag.status  <- get_attr(dat, "diag.status")

  ## Parameters ##
  tx.efficacy  <- get_param(dat, "tx.efficacy")
  ept.efficacy <- get_param(dat, "ept.efficacy")
  pn.test.prob <- get_param(dat, "pn.test.prob")
  pn.arm       <- get_param(dat, "pn.arm")

  ## Counters ##
  n.index.tx           <- 0
  n.partner.tx         <- 0   # partners actually treated (any state)
  n.partner.tx.wasted  <- 0   # partners treated who were uninfected
  n.partner.cleared.pn <- 0   # partners moved I -> S via PN this step

  ## --- Index treatment (always on; one path for every fresh dx) --- ##
  idsIndex <- which(active == 1 & dx.this.step == 1)
  if (length(idsIndex) > 0) {
    vec <- which(rbinom(length(idsIndex), 1, tx.efficacy) == 1)
    if (length(vec) > 0) {
      idsClear <- idsIndex[vec]
      idsClear <- idsClear[status[idsClear] == "i"]
      n.index.tx <- length(idsClear)
      if (n.index.tx > 0) {
        status[idsClear]       <- "s"
        tx.this.step[idsClear] <- 1
      }
    }
  }

  ## --- Notified-partner treatment --- ##
  idsPart <- which(active == 1 & pn.notified == at)

  if (length(idsPart) > 0 && pn.arm %in% c("PR", "EPT")) {

    if (pn.arm == "PR") {
      # Patient Referral: test first, then treat infecteds only.
      # Uninfected notified partners are not treated.
      idsPartI <- idsPart[status[idsPart] == "i"]
      nPartI <- length(idsPartI)

      if (nPartI > 0) {
        # Test sensitivity
        testPos <- which(rbinom(nPartI, 1, pn.test.prob) == 1)
        idsTestPos <- idsPartI[testPos]
        # Of those who test positive, treat at tx.efficacy
        if (length(idsTestPos) > 0) {
          rxOk <- which(rbinom(length(idsTestPos), 1, tx.efficacy) == 1)
          idsRx <- idsTestPos[rxOk]
          n.partner.tx <- length(idsRx)
          if (n.partner.tx > 0) {
            status[idsRx]       <- "s"
            diag.status[idsRx]  <- 1   # they were tested and treated
            tx.this.step[idsRx] <- 1
            n.partner.cleared.pn <- n.partner.tx
          }
        }
      }

    } else if (pn.arm == "EPT") {
      # Expedited Partner Therapy: dispense to all notified partners
      # at ept.efficacy. Uninfected partners count as wasted doses.
      rxOk <- which(rbinom(length(idsPart), 1, ept.efficacy) == 1)
      idsRx <- idsPart[rxOk]
      n.partner.tx <- length(idsRx)

      if (n.partner.tx > 0) {
        tx.this.step[idsRx] <- 1
        # Wasted = treated but not infected
        isI <- status[idsRx] == "i"
        n.partner.tx.wasted <- sum(!isI)
        idsClear <- idsRx[isI]
        if (length(idsClear) > 0) {
          status[idsClear]       <- "s"
          n.partner.cleared.pn   <- length(idsClear)
        }
      }
    }
  }

  ## Write attributes ##
  dat <- set_attr(dat, "status",       status)
  dat <- set_attr(dat, "tx.this.step", tx.this.step)
  dat <- set_attr(dat, "diag.status",  diag.status)

  ## Summary statistics ##
  dat <- set_epi(dat, "n.index.tx",            at, n.index.tx)
  dat <- set_epi(dat, "n.partner.tx",          at, n.partner.tx)
  dat <- set_epi(dat, "n.partner.tx.wasted",   at, n.partner.tx.wasted)
  dat <- set_epi(dat, "n.partner.cleared.pn",  at, n.partner.cleared.pn)

  return(dat)
}


# Recovery Module ----------------------------------------------------------

recov <- function(dat, at) {
  # Natural clearance I -> S at the slow rec.rate. Treatment-driven
  # clearance is handled in treat(), so this module is purely the
  # background recovery process.

  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  rec.rate <- get_param(dat, "rec.rate")

  nRec <- 0
  idsElig <- which(active == 1 & status == "i")
  if (length(idsElig) > 0) {
    vec <- which(rbinom(length(idsElig), 1, rec.rate) == 1)
    if (length(vec) > 0) {
      idsRec <- idsElig[vec]
      nRec <- length(idsRec)
      status[idsRec] <- "s"
    }
  }

  dat <- set_attr(dat, "status", status)
  dat <- set_epi(dat, "is.flow.natural", at, nRec)
  return(dat)
}


# Infection Module ---------------------------------------------------------

infect <- function(dat, at) {
  # Standard discordant-edge transmission, with a custom hook so we can
  # increment the per-node `infections` counter on every fresh S -> I.
  # Tracking that lets us compute the reinfection-per-index distribution
  # after the simulation finishes.

  active    <- get_attr(dat, "active")
  status    <- get_attr(dat, "status")
  infTime   <- get_attr(dat, "infTime")
  infections <- get_attr(dat, "infections")

  inf.prob <- get_param(dat, "inf.prob")
  act.rate <- get_param(dat, "act.rate")

  idsInf <- which(active == 1 & status == "i")
  nActive <- sum(active == 1)
  nElig <- length(idsInf)
  nInf <- 0

  if (nElig > 0 && nElig < nActive) {
    del <- discord_edgelist(dat, at)
    if (!is.null(del) && nrow(del) > 0) {
      del$transProb <- inf.prob
      del$actRate   <- act.rate
      del$finalProb <- 1 - (1 - del$transProb)^del$actRate
      transmit <- rbinom(nrow(del), 1, del$finalProb)
      del <- del[which(transmit == 1), ]
      idsNewInf <- unique(del$sus)
      nInf <- length(idsNewInf)
      if (nInf > 0) {
        status[idsNewInf]     <- "i"
        infTime[idsNewInf]    <- at
        infections[idsNewInf] <- infections[idsNewInf] + 1
      }
    }
  }

  if (nInf > 0) {
    dat <- set_attr(dat, "status",     status)
    dat <- set_attr(dat, "infTime",    infTime)
    dat <- set_attr(dat, "infections", infections)
  }

  dat <- set_epi(dat, "si.flow", at, nInf)
  return(dat)
}
