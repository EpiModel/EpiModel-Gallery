
##
## SEIR Model with All-or-Nothing Vaccination and Vital Dynamics
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##
## Authors: Connor M. Van Meter (Emory University)
## Date: October 2018
##


# Vaccine Attribute Initialization Module ----------------------------------

init_vaccine_attrs <- function(dat, at) {
  # Initialize vaccination and protection attributes for all nodes.
  # Runs once at the first module call. Uses override.null.error to detect
  # whether the "vaccination" attribute has been set yet.
  #
  # The starting population may have pre-existing vaccination coverage:
  #   vaccination.rate.initialization: probability each node is vaccinated
  #   protection.rate.initialization: probability vaccine confers protection
  #     (among vaccinated susceptibles only)
  #
  # In the all-or-nothing model, protection is binary: protected individuals
  # are moved to status = "v" (vaccine-immune) in the arrival module and
  # cannot be infected. Vaccinated but unprotected individuals receive
  # protection = "none" and remain fully susceptible.

  if (is.null(get_attr(dat, "vaccination", override.null.error = TRUE))) {

    ## Attributes ##
    active <- get_attr(dat, "active")
    n <- sum(active == 1)
    status <- get_attr(dat, "status")

    ## Parameters ##
    vaccination.rate.init <- get_param(dat, "vaccination.rate.initialization")
    protection.rate.init <- get_param(dat, "protection.rate.initialization")

    ## Initialize vectors (NA = unvaccinated/unprotected) ##
    vaccination <- rep(NA, n)
    protection <- rep(NA, n)

    ## Stochastic vaccination of initial population ##
    idsEligVacInit <- which(active == 1)
    vecVacInit <- rbinom(length(idsEligVacInit), 1, vaccination.rate.init)
    idsVacInit <- idsEligVacInit[which(vecVacInit == 1)]
    vaccination[idsVacInit] <- "initial"

    ## Stochastic protection among vaccinated susceptibles ##
    # Only vaccinated susceptible individuals can receive protection.
    # Vaccinated individuals who are already infected gain no benefit.
    idsEligProtInit <- which(vaccination == "initial" & status == "s")
    vecProtInit <- rbinom(length(idsEligProtInit), 1, protection.rate.init)
    idsProtInit <- idsEligProtInit[which(vecProtInit == 1)]
    idsNoProtInit <- setdiff(idsVacInit, idsProtInit)
    protection[idsProtInit] <- "initial"
    protection[idsNoProtInit] <- "none"

    ## Write out attributes ##
    dat <- set_attr(dat, "vaccination", vaccination)
    dat <- set_attr(dat, "protection", protection)
  }

  return(dat)
}


# Infection Module ---------------------------------------------------------

infect <- function(dat, at) {
  # Simulate S -> E transmission along discordant edges.
  #
  # Protected individuals (status = "v") are excluded from the susceptible
  # pool by discord_edgelist(), which only considers nodes with status "s"
  # as susceptible. This is how all-or-nothing vaccination provides
  # complete immunity -- protected individuals never appear on the
  # discordant edgelist and cannot be infected.

  ## Attributes ##
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  infTime <- get_attr(dat, "infTime")

  ## Parameters ##
  inf.prob <- get_param(dat, "inf.prob")
  act.rate <- get_param(dat, "act.rate")

  ## Find infected nodes ##
  idsInf <- which(active == 1 & status == "i")
  nActive <- sum(active == 1)
  nElig <- length(idsInf)

  ## Initialize incidence counter ##
  nInf <- 0

  ## Transmission process ##
  if (nElig > 0 && nElig < nActive) {

    del <- discord_edgelist(dat, at)

    if (!(is.null(del)) && nrow(del) > 0) {

      # Per-timestep transmission probability:
      # P(transmit) = 1 - (1 - inf.prob)^act.rate
      del$transProb <- inf.prob
      del$actRate <- act.rate
      del$finalProb <- 1 - (1 - del$transProb)^del$actRate

      # Stochastic transmission (Bernoulli trial per discordant edge)
      transmit <- rbinom(nrow(del), 1, del$finalProb)
      del <- del[which(transmit == 1), ]

      idsNewInf <- unique(del$sus)
      nInf <- length(idsNewInf)

      # Newly infected enter the exposed (E) compartment
      if (nInf > 0) {
        status[idsNewInf] <- "e"
        dat <- set_attr(dat, "status", status)
        infTime[idsNewInf] <- at
        dat <- set_attr(dat, "infTime", infTime)
      }
    }
  }

  ## Summary statistics ##
  dat <- set_epi(dat, "se.flow", at, nInf)

  dat <- set_epi(dat, "s.num", at, sum(active == 1 & status == "s"))
  dat <- set_epi(dat, "e.num", at, sum(active == 1 & status == "e"))
  dat <- set_epi(dat, "i.num", at, sum(active == 1 & status == "i"))
  dat <- set_epi(dat, "r.num", at, sum(active == 1 & status == "r"))
  dat <- set_epi(dat, "v.num", at, sum(active == 1 & status == "v"))

  return(dat)
}


# Disease Progression Module -----------------------------------------------

progress <- function(dat, at) {
  # Simulate disease progression: E -> I at ei.rate, I -> R at ir.rate.
  # No R -> S transition (SEIR, not SEIRS) -- recovered individuals retain
  # permanent natural immunity.

  ## Attributes ##
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")

  ## Parameters ##
  ei.rate <- get_param(dat, "ei.rate")
  ir.rate <- get_param(dat, "ir.rate")

  ## E -> I (exposed to infectious) ##
  nProg <- 0
  idsEligProg <- which(active == 1 & status == "e")
  nEligProg <- length(idsEligProg)

  if (nEligProg > 0) {
    vecProg <- which(rbinom(nEligProg, 1, ei.rate) == 1)
    if (length(vecProg) > 0) {
      idsProg <- idsEligProg[vecProg]
      nProg <- length(idsProg)
      status[idsProg] <- "i"
    }
  }

  ## I -> R (infectious to recovered) ##
  nRec <- 0
  idsEligRec <- which(active == 1 & status == "i")
  nEligRec <- length(idsEligRec)

  if (nEligRec > 0) {
    vecRec <- which(rbinom(nEligRec, 1, ir.rate) == 1)
    if (length(vecRec) > 0) {
      idsRec <- idsEligRec[vecRec]
      nRec <- length(idsRec)
      status[idsRec] <- "r"
    }
  }

  ## Write out updated status ##
  dat <- set_attr(dat, "status", status)

  ## Summary statistics ##
  dat <- set_epi(dat, "ei.flow", at, nProg)
  dat <- set_epi(dat, "ir.flow", at, nRec)

  return(dat)
}


# Departure Module ---------------------------------------------------------

dfunc <- function(dat, at) {
  # Simulate departures (mortality) with disease-induced excess mortality.
  # All active nodes face a baseline departure.rate per timestep. Infected
  # individuals (status = "i") have their rate multiplied by
  # departure.disease.mult, representing excess disease-induced mortality.

  ## Attributes ##
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  exitTime <- get_attr(dat, "exitTime")

  ## Parameters ##
  departure.rate <- get_param(dat, "departure.rate")
  dep.rates <- rep(departure.rate, length(active))
  dep.dis.mult <- get_param(dat, "departure.disease.mult")

  ## Determine eligible (alive) nodes ##
  idsElig <- which(active == 1)
  nElig <- length(idsElig)
  nDeaths <- 0

  if (nElig > 0) {

    dep_rates_of_elig <- dep.rates[idsElig]

    # Disease-induced excess mortality for infected individuals
    idsElig.inf <- which(status[idsElig] == "i")
    dep_rates_of_elig[idsElig.inf] <- dep.rates[idsElig.inf] * dep.dis.mult

    ## Stochastic departure process ##
    vecDeaths <- which(rbinom(nElig, 1, dep_rates_of_elig) == 1)
    idsDeaths <- idsElig[vecDeaths]
    nDeaths <- length(idsDeaths)

    if (nDeaths > 0) {
      active[idsDeaths] <- 0
      exitTime[idsDeaths] <- at
    }
  }

  ## Write out updated attributes ##
  dat <- set_attr(dat, "active", active)
  dat <- set_attr(dat, "exitTime", exitTime)

  ## Summary statistics ##
  dat <- set_epi(dat, "d.flow", at, nDeaths)

  return(dat)
}


# Arrival Module -----------------------------------------------------------

afunc <- function(dat, at) {
  # Simulate arrivals (births) with all-or-nothing vaccination.
  #
  # Two vaccination routes operate each timestep:
  #
  # 1. PROGRESSION: Unvaccinated active nodes may receive vaccination at
  #    vaccination.rate.progression. Among newly vaccinated susceptibles,
  #    protection is conferred at protection.rate.progression. This
  #    represents ongoing vaccination campaigns in the existing population.
  #
  # 2. ARRIVALS: New nodes enter as susceptible. Each may be vaccinated at
  #    vaccination.rate.arrivals, with protection at protection.rate.arrivals.
  #    This represents newborn vaccination programs.
  #
  # After all vaccination processing, any susceptible node with vaccine
  # protection is moved to status = "v" (vaccine-immune). This is the
  # all-or-nothing mechanism: protected individuals are completely immune
  # and cannot be infected.
  #
  # Key assumption: individuals may only be vaccinated once. Those who are
  # vaccinated but not protected (protection = "none") will not have
  # another opportunity to become protected.

  ## Attributes ##
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  infTime <- get_attr(dat, "infTime")

  ## Parameters ##
  nActive <- sum(active == 1)
  nTotal <- length(active)
  a.rate <- get_param(dat, "arrival.rate")
  vaccination.rate.arrivals <- get_param(dat, "vaccination.rate.arrivals")
  protection.rate.arrivals <- get_param(dat, "protection.rate.arrivals")
  vaccination.rate.progression <- get_param(dat, "vaccination.rate.progression")
  protection.rate.progression <- get_param(dat, "protection.rate.progression")

  ## Initialize flow counters ##
  nVax.prog <- 0
  nPrt.prog <- 0
  nVax.arrival <- 0
  nPrt.arrival <- 0

  vaccination <- get_attr(dat, "vaccination")
  protection <- get_attr(dat, "protection")

  ## --- VACCINATION PROGRESSION --- ##
  # Unvaccinated active nodes may receive vaccination each timestep
  idsEligVacProg <- which(is.na(vaccination) & active == 1)
  vecVacProg <- rbinom(length(idsEligVacProg), 1, vaccination.rate.progression)
  idsVacProg <- idsEligVacProg[which(vecVacProg == 1)]
  vaccination[idsVacProg] <- "progress"

  # Among newly vaccinated susceptibles, stochastic protection
  idsEligProtProg <- which(vaccination == "progress" & is.na(protection) & status == "s")
  vecProtProg <- rbinom(length(idsEligProtProg), 1, protection.rate.progression)
  idsProtProg <- idsEligProtProg[which(vecProtProg == 1)]
  idsNoProtProg <- setdiff(idsVacProg, idsProtProg)
  protection[idsProtProg] <- "progress"
  protection[idsNoProtProg] <- "none"

  nVax.prog <- length(idsVacProg)
  nPrt.prog <- length(idsProtProg)

  ## --- ARRIVAL PROCESS --- ##
  nArrivalsExp <- nActive * a.rate
  nArrivals <- rpois(1, nArrivalsExp)

  if (nArrivals > 0) {
    dat <- append_core_attr(dat, at, nArrivals)

    newNodes <- (nTotal + 1):(nTotal + nArrivals)
    status <- c(status, rep("s", nArrivals))
    infTime <- c(infTime, rep(NA, nArrivals))
    vaccination <- c(vaccination, rep(NA, nArrivals))
    protection <- c(protection, rep(NA, nArrivals))

    # Vaccinate new arrivals
    vaccinatedNewArrivals <- rbinom(nArrivals, 1, vaccination.rate.arrivals)
    vaccination[newNodes] <- ifelse(vaccinatedNewArrivals == 1, "arrival", NA)
    nVax.arrival <- sum(vaccinatedNewArrivals == 1)

    # Confer protection to vaccinated susceptible arrivals
    idsEligProtArrival <- which(
      vaccination == "arrival" & status == "s" & is.na(protection)
    )
    vecProtArrival <- rbinom(length(idsEligProtArrival), 1, protection.rate.arrivals)
    idsProtArrival <- idsEligProtArrival[which(vecProtArrival == 1)]
    idsNoProtArrival <- idsEligProtArrival[which(vecProtArrival == 0)]
    protection[idsProtArrival] <- "arrival"
    protection[idsNoProtArrival] <- "none"
    nPrt.arrival <- sum(vecProtArrival == 1)
  }

  ## --- UPDATE STATUS: protected susceptibles -> vaccine-immune --- ##
  # This is the all-or-nothing mechanism: any susceptible with vaccine
  # protection is moved to status = "v" (complete immunity to infection).
  active <- get_attr(dat, "active")
  statusV <- which(
    status == "s" & protection %in% c("initial", "progress", "arrival") &
      active == 1
  )
  status[statusV] <- "v"

  ## Write out updated attributes ##
  dat <- set_attr(dat, "status", status)
  dat <- set_attr(dat, "vaccination", vaccination)
  dat <- set_attr(dat, "protection", protection)
  dat <- set_attr(dat, "infTime", infTime)

  ## Summary statistics ##
  dat <- set_epi(dat, "a.flow", at, nArrivals)
  dat <- set_epi(dat, "vac.flow", at, nVax.prog + nVax.arrival)
  dat <- set_epi(dat, "prt.flow", at, nPrt.prog + nPrt.arrival)
  dat <- set_epi(dat, "vac.num", at,
                 sum(active == 1 & vaccination %in%
                       c("initial", "progress", "arrival"), na.rm = TRUE))
  dat <- set_epi(dat, "prt.num", at,
                 sum(active == 1 & protection %in%
                       c("initial", "progress", "arrival"), na.rm = TRUE))

  return(dat)
}
