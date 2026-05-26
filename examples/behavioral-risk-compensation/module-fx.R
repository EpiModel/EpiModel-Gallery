
##
## SIR with Behavioral Risk Compensation During Illness
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##
## Author: Samuel M. Jenness (Emory University)
## Date: May 2026
##


# Replacement infection module ---------------------------------------------

infect <- function(dat, at) {
  # Transmission along discordant edges, with the per-edge act count
  # modified by the infected partner's sub-stage. Effective acts =
  # act.rate * mult, where mult depends on whether the infected partner is
  # currently in the early (symptomatic, reduced contacts) or late
  # (recovering, contacts returning to baseline) sub-stage of infection.
  # Setting both multipliers to 1 recovers the standard behavior-naive
  # model.

  ## Idempotent initialization of the inf.stage attribute. ##
  # Initial seeds enter with status == "i" but no inf.stage. Place them
  # all in the early sub-stage so the very first transmissions use the
  # index-period multiplier. Whichever of infect() / progress() runs
  # first on at == 2 creates the attribute.
  if (is.null(get_attr(dat, "inf.stage", override.null.error = TRUE))) {
    status0 <- get_attr(dat, "status")
    dat <- set_attr(dat, "inf.stage",
                    ifelse(status0 == "i", "early", NA))
  }

  ## Attributes ##
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  inf.stage <- get_attr(dat, "inf.stage")
  infTime <- get_attr(dat, "infTime")

  ## Parameters ##
  inf.prob <- get_param(dat, "inf.prob")
  act.rate <- get_param(dat, "act.rate")
  mult.early <- get_param(dat, "mult.early")
  mult.late <- get_param(dat, "mult.late")

  idsInf <- which(active == 1 & status == "i")
  nActive <- sum(active == 1)
  nElig <- length(idsInf)
  nInf <- 0

  if (nElig > 0 && nElig < nActive) {

    del <- discord_edgelist(dat, at)

    if (!is.null(del) && nrow(del) > 0) {
      # Per-edge multiplier indexed by the infected partner's sub-stage.
      # discord_edgelist() places the infected node id in del$inf.
      stg <- inf.stage[del$inf]
      mult <- ifelse(stg == "early", mult.early, mult.late)

      del$transProb <- inf.prob
      del$actRate <- act.rate * mult
      del$finalProb <- 1 - (1 - del$transProb)^del$actRate

      transmit <- rbinom(nrow(del), 1, del$finalProb)
      del <- del[which(transmit == 1), ]
      idsNewInf <- unique(del$sus)
      nInf <- length(idsNewInf)

      if (nInf > 0) {
        status[idsNewInf] <- "i"
        inf.stage[idsNewInf] <- "early"
        infTime[idsNewInf] <- at
        dat <- set_attr(dat, "status", status)
        dat <- set_attr(dat, "inf.stage", inf.stage)
        dat <- set_attr(dat, "infTime", infTime)
      }
    }
  }

  dat <- set_epi(dat, "si.flow", at, nInf)
  return(dat)
}


# Disease progression module -----------------------------------------------

progress <- function(dat, at) {
  # Two stochastic stage transitions per timestep:
  #   inf.stage early -> late at rate el.rate
  #   status   i -> r when in late stage, at rate lr.rate
  # Geometric durations with means 1/el.rate and 1/lr.rate.

  # Idempotent initialization of the inf.stage attribute. Whichever of
  # infect() / progress() runs first on at == 2 creates the attribute;
  # the second sees it already present.
  if (is.null(get_attr(dat, "inf.stage", override.null.error = TRUE))) {
    status0 <- get_attr(dat, "status")
    dat <- set_attr(dat, "inf.stage",
                    ifelse(status0 == "i", "early", NA))
  }

  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  inf.stage <- get_attr(dat, "inf.stage")

  el.rate <- get_param(dat, "el.rate")
  lr.rate <- get_param(dat, "lr.rate")

  ## early -> late ##
  nEL <- 0
  idsEL <- which(active == 1 & status == "i" & inf.stage == "early")
  if (length(idsEL) > 0) {
    vec <- which(rbinom(length(idsEL), 1, el.rate) == 1)
    if (length(vec) > 0) {
      ids <- idsEL[vec]
      nEL <- length(ids)
      inf.stage[ids] <- "late"
    }
  }

  ## late -> recovered ##
  nLR <- 0
  idsLR <- which(active == 1 & status == "i" & inf.stage == "late")
  if (length(idsLR) > 0) {
    vec <- which(rbinom(length(idsLR), 1, lr.rate) == 1)
    if (length(vec) > 0) {
      ids <- idsLR[vec]
      nLR <- length(ids)
      status[ids] <- "r"
      inf.stage[ids] <- NA
    }
  }

  dat <- set_attr(dat, "status", status)
  dat <- set_attr(dat, "inf.stage", inf.stage)

  ## Summary statistics ##
  dat <- set_epi(dat, "el.flow", at, nEL)
  dat <- set_epi(dat, "ir.flow", at, nLR)
  dat <- set_epi(dat, "ie.num", at,
                 sum(active == 1 & status == "i" & inf.stage == "early",
                     na.rm = TRUE))
  dat <- set_epi(dat, "il.num", at,
                 sum(active == 1 & status == "i" & inf.stage == "late",
                     na.rm = TRUE))
  dat <- set_epi(dat, "r.num", at, sum(active == 1 & status == "r"))
  return(dat)
}
