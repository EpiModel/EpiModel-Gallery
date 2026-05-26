
##
## SIR Model with Time-Varying (Phased) Vaccination
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##
## Author: Samuel M. Jenness (Emory University)
## Date: May 2026
##


# Infection Module ---------------------------------------------------------

infect <- function(dat, at) {
  # Standard S -> I transmission along discordant edges. Structurally
  # identical to EpiModel's built-in SI infection module, but written out
  # here so the SIR + V state machine is fully visible in one place.
  # Vaccine-immune individuals (status = "v") are excluded automatically:
  # discord_edgelist() only pairs "s" with "i".
  #
  # Optional exog.inf.prob adds an external (exogenous) force of infection.
  # At each timestep every active susceptible has probability exog.inf.prob
  # of being infected by a source outside the modeled network: spillover
  # from an unmodeled reservoir, environmental exposure, or contacts with
  # individuals not represented in the simulated network. This is *not*
  # importation in the traditional epidemiological sense. Importation
  # describes already-infected individuals arriving in the population,
  # which is an arrivals process and belongs in a vital-dynamics module.
  # Here the population is closed and existing susceptibles transition to
  # infectious at a constant background hazard that is independent of
  # in-population prevalence. Setting exog.inf.prob = 0 disables the
  # exogenous channel and leaves a purely network-driven model.

  ## Attributes ##
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  infTime <- get_attr(dat, "infTime")

  ## Parameters ##
  inf.prob <- get_param(dat, "inf.prob")
  act.rate <- get_param(dat, "act.rate")
  exog.inf.prob <- get_param(dat, "exog.inf.prob")

  ## Find infected nodes ##
  idsInf <- which(active == 1 & status == "i")
  nActive <- sum(active == 1)
  nElig <- length(idsInf)
  nInf <- 0
  nExog <- 0

  ## Endogenous (network) transmission ##
  if (nElig > 0 && nElig < nActive) {
    del <- discord_edgelist(dat, at)
    if (!is.null(del) && nrow(del) > 0) {
      del$transProb <- inf.prob
      del$actRate <- act.rate
      del$finalProb <- 1 - (1 - del$transProb)^del$actRate
      transmit <- rbinom(nrow(del), 1, del$finalProb)
      del <- del[which(transmit == 1), ]
      idsNewInf <- unique(del$sus)
      nInf <- length(idsNewInf)
      if (nInf > 0) {
        status[idsNewInf] <- "i"
        infTime[idsNewInf] <- at
      }
    }
  }

  ## Exogenous force of infection ##
  if (exog.inf.prob > 0) {
    idsSus <- which(active == 1 & status == "s")
    nSus <- length(idsSus)
    if (nSus > 0) {
      vecExog <- which(rbinom(nSus, 1, exog.inf.prob) == 1)
      nExog <- length(vecExog)
      if (nExog > 0) {
        idsExog <- idsSus[vecExog]
        status[idsExog] <- "i"
        infTime[idsExog] <- at
      }
    }
  }

  if (nInf + nExog > 0) {
    dat <- set_attr(dat, "status", status)
    dat <- set_attr(dat, "infTime", infTime)
  }

  dat <- set_epi(dat, "si.flow", at, nInf)
  dat <- set_epi(dat, "exog.flow", at, nExog)

  return(dat)
}


# Recovery Module ----------------------------------------------------------

recov <- function(dat, at) {
  # I -> R at a fixed per-timestep rate (rec.rate).
  #
  # Optional R -> S waning of natural immunity at rate rs.rate. When
  # rs.rate = 0, this gives standard SIR with permanent natural immunity.
  # When rs.rate > 0, the model becomes SIRS, which allows sustained
  # endemic dynamics that the reactive intervention can respond to over
  # many cycles.

  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  rec.rate <- get_param(dat, "rec.rate")
  rs.rate <- get_param(dat, "rs.rate")

  nRec <- 0
  idsElig <- which(active == 1 & status == "i")
  nElig <- length(idsElig)
  if (nElig > 0) {
    vec <- which(rbinom(nElig, 1, rec.rate) == 1)
    if (length(vec) > 0) {
      idsRec <- idsElig[vec]
      nRec <- length(idsRec)
      status[idsRec] <- "r"
    }
  }

  nWane <- 0
  if (rs.rate > 0) {
    idsEligWane <- which(active == 1 & status == "r")
    nEligWane <- length(idsEligWane)
    if (nEligWane > 0) {
      vecWane <- which(rbinom(nEligWane, 1, rs.rate) == 1)
      if (length(vecWane) > 0) {
        idsWane <- idsEligWane[vecWane]
        nWane <- length(idsWane)
        status[idsWane] <- "s"
      }
    }
  }

  dat <- set_attr(dat, "status", status)
  dat <- set_epi(dat, "ir.flow", at, nRec)
  dat <- set_epi(dat, "rs.flow", at, nWane)
  dat <- set_epi(dat, "r.num", at, sum(active == 1 & status == "r"))

  return(dat)
}


# Vaccination Module -------------------------------------------------------

vaccinate <- function(dat, at) {
  # Time-varying vaccination intervention. At each timestep the module
  # decides whether vaccination is currently *active* based on the
  # user-supplied schedule, then applies the per-S vaccination probability
  # to eligible susceptibles.
  #
  # Two scheduling modes are supported:
  #
  #   1. Windowed mode (default). vax.starts and vax.ends are parallel
  #      numeric vectors defining one or more (start, end) activation
  #      windows. Vaccination is active whenever `at` lies inside any
  #      window. Passing vax.starts = -1 (and vax.ends = -1) effectively
  #      disables the program because no timestep ever satisfies
  #      at >= -1 & at <= -1 in a normal simulation.
  #
  #   2. Reactive mode. Set vax.prev.on to a positive prevalence threshold
  #      to enable. Vaccination activates whenever the current prevalence
  #      crosses above vax.prev.on and deactivates only when it falls
  #      below vax.prev.off. The two thresholds create *hysteresis* and
  #      prevent the program from rapidly toggling on and off near the
  #      activation threshold.
  #
  # The vaccine itself is "all-or-nothing": each vaccinated susceptible
  # becomes immune (status "v") with probability vax.efficacy. Vaccinated
  # but unprotected individuals stay susceptible.
  #
  # Optional waning. When vax.wane > 0, vaccine-immune individuals revert
  # to susceptible at that per-timestep rate. The vax.attempted flag is
  # reset so they can be re-vaccinated by future activations. Setting
  # vax.wane = 0 gives the standard SIR + V dynamics with permanent vaccine
  # immunity.

  ## Attributes ##
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  vax.attempted <- get_attr(dat, "vax.attempted",
                            override.null.error = TRUE)

  ## Parameters ##
  vax.rate <- get_param(dat, "vax.rate")
  vax.efficacy <- get_param(dat, "vax.efficacy")
  vax.wane <- get_param(dat, "vax.wane")
  vax.starts <- get_param(dat, "vax.starts")
  vax.ends <- get_param(dat, "vax.ends")
  vax.prev.on <- get_param(dat, "vax.prev.on")
  vax.prev.off <- get_param(dat, "vax.prev.off")

  ## Lazy initialization of per-node vax.attempted attribute ##
  if (is.null(vax.attempted)) {
    vax.attempted <- rep(0, length(active))
    dat <- set_attr(dat, "vax.attempted", vax.attempted)
  }

  ## --- Vaccine waning: V -> S --- ##
  # Runs before the activity check so newly waned individuals are eligible
  # to be re-vaccinated in this same timestep if the program is active.
  nWane <- 0
  if (vax.wane > 0) {
    idsEligWane <- which(active == 1 & status == "v")
    nEligWane <- length(idsEligWane)
    if (nEligWane > 0) {
      vecWane <- which(rbinom(nEligWane, 1, vax.wane) == 1)
      if (length(vecWane) > 0) {
        idsWane <- idsEligWane[vecWane]
        nWane <- length(idsWane)
        status[idsWane] <- "s"
        vax.attempted[idsWane] <- 0
      }
    }
  }

  ## --- Determine whether vaccination is active this timestep --- ##
  if (vax.prev.on > 0) {
    # Reactive mode with hysteresis
    nAct <- sum(active == 1)
    current.prev <- if (nAct > 0) {
      sum(active == 1 & status == "i") / nAct
    } else {
      0
    }
    prev.active <- 0
    if (at > 2) {
      val <- get_epi(dat, "vax.active", at = at - 1,
                     override.null.error = TRUE)
      if (!is.null(val) && !is.na(val)) prev.active <- val
    }
    if (prev.active == 1) {
      vax.active <- as.numeric(current.prev > vax.prev.off) # <1>
    } else {
      vax.active <- as.numeric(current.prev > vax.prev.on)
    }
  } else {
    # Windowed mode: active if `at` lies in any (start, end) window
    vax.active <- as.numeric(any(at >= vax.starts & at <= vax.ends)) # <2>
  }

  ## --- Apply vaccination when active --- ##
  nVax <- 0       # number vaccinated this timestep
  nImmune <- 0    # number who became vaccine-immune

  if (vax.active == 1 && vax.rate > 0) {
    idsElig <- which(active == 1 & status == "s" & vax.attempted == 0)
    nElig <- length(idsElig)
    if (nElig > 0) {
      vecVax <- which(rbinom(nElig, 1, vax.rate) == 1)
      if (length(vecVax) > 0) {
        idsVax <- idsElig[vecVax]
        nVax <- length(idsVax)
        vax.attempted[idsVax] <- 1
        # All-or-nothing protection: stochastic at vax.efficacy
        vecImmune <- which(rbinom(nVax, 1, vax.efficacy) == 1)
        idsImmune <- idsVax[vecImmune]
        nImmune <- length(idsImmune)
        if (nImmune > 0) {
          status[idsImmune] <- "v"
        }
      }
    }
  }

  dat <- set_attr(dat, "status", status)
  dat <- set_attr(dat, "vax.attempted", vax.attempted)

  ## --- Summary statistics --- ##
  dat <- set_epi(dat, "vax.active", at, vax.active)
  dat <- set_epi(dat, "vax.flow", at, nVax)
  dat <- set_epi(dat, "sv.flow", at, nImmune)
  dat <- set_epi(dat, "vs.flow", at, nWane)
  dat <- set_epi(dat, "v.num", at, sum(active == 1 & status == "v"))

  return(dat)
}
