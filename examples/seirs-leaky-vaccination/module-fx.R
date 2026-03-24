
##
## SEIRS Model with Leaky Vaccination and Vital Dynamics
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##
## Authors: Connor M. Van Meter (Emory University)
## Date: December 2018
##


# Infection Module with Leaky Vaccine Protection ---------------------------

infect <- function(dat, at) {
  # Simulate S -> E transmission along discordant edges, with reduced
  # transmission probability for vaccine-protected susceptibles.
  #
  # In a leaky vaccine model, protected individuals remain susceptible but
  # face a reduced force of infection:
  #   reduced inf.prob = (1 - vaccine.efficacy) * inf.prob
  #
  # This contrasts with the all-or-nothing model where protected individuals
  # are completely immune (status = "v"). Here, protected individuals keep
  # status = "s" and can still be infected -- just at a lower rate.
  #
  # Protection persists through the SEIRS cycle: if a protected individual
  # is infected, progresses through E -> I -> R -> S, they retain their
  # protection attribute and benefit from reduced FOI again.

  ## Attributes ##
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  infTime <- get_attr(dat, "infTime")
  protection <- get_attr(dat, "protection")

  ## Parameters ##
  inf.prob <- get_param(dat, "inf.prob")
  act.rate <- get_param(dat, "act.rate")
  vaccine.efficacy <- get_param(dat, "vaccine.efficacy")

  ## Find infected nodes ##
  idsInf <- which(active == 1 & status == "i")
  nActive <- sum(active == 1)
  nElig <- length(idsInf)

  ## Initialize incidence counter ##
  nInf <- 0

  ## Transmission process ##
  if (nElig > 0 && nElig < nActive) {

    del <- discord_edgelist(dat, at)

    if (!(is.null(del))) {

      # Determine transmission probability for each discordant edge.
      # Vaccine-protected susceptibles have reduced probability.
      inf.prob.reduced <- (1 - vaccine.efficacy) * inf.prob
      idsProtected <- which(
        !is.na(protection) & protection != "none" &
          status == "s" & active == 1
      )

      # Build lookup: protected susceptible IDs -> reduced transProb
      protDF <- data.frame(
        sus = idsProtected,
        transProb = rep(inf.prob.reduced, length(idsProtected))
      )
      del <- merge(del, protDF, by = "sus", all.x = TRUE)

      # Unprotected susceptibles use the full inf.prob
      del$transProb <- ifelse(is.na(del$transProb), inf.prob, del$transProb)

      # Per-timestep transmission probability:
      # P(transmit) = 1 - (1 - transProb)^act.rate
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

  return(dat)
}


# Disease Progression Module -----------------------------------------------

progress <- function(dat, at) {
  # Simulate disease progression: E -> I -> R -> S (SEIRS).
  # The R -> S transition represents waning natural immunity -- recovered
  # individuals eventually become susceptible again. If they have vaccine
  # protection, they retain it and face a reduced force of infection
  # upon re-entering the susceptible pool.

  ## Attributes ##
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")

  ## Parameters ##
  ei.rate <- get_param(dat, "ei.rate")
  ir.rate <- get_param(dat, "ir.rate")
  rs.rate <- get_param(dat, "rs.rate")

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

  ## R -> S (recovered to susceptible, waning immunity) ##
  nSus <- 0
  idsEligSus <- which(active == 1 & status == "r")
  nEligSus <- length(idsEligSus)

  if (nEligSus > 0) {
    vecSus <- which(rbinom(nEligSus, 1, rs.rate) == 1)
    if (length(vecSus) > 0) {
      idsSus <- idsEligSus[vecSus]
      nSus <- length(idsSus)
      status[idsSus] <- "s"
    }
  }

  ## Write out updated status ##
  dat <- set_attr(dat, "status", status)

  ## Summary statistics ##
  dat <- set_epi(dat, "ei.flow", at, nProg)
  dat <- set_epi(dat, "ir.flow", at, nRec)
  dat <- set_epi(dat, "rs.flow", at, nSus)

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
  departure.rates <- rep(departure.rate, length(active))
  departure.dis.mult <- get_param(dat, "departure.disease.mult")

  ## Determine eligible (alive) nodes ##
  idsElig <- which(active == 1)
  nElig <- length(idsElig)
  nDepartures <- 0

  if (nElig > 0) {

    dep_rates_of_elig <- departure.rates[idsElig]

    # Disease-induced excess mortality for infected individuals
    idsElig.inf <- which(status[idsElig] == "i")
    dep_rates_of_elig[idsElig.inf] <- departure.rates[idsElig.inf] *
      departure.dis.mult

    ## Stochastic departure process ##
    vecDeparture <- which(rbinom(nElig, 1, dep_rates_of_elig) == 1)
    idsDeparture <- idsElig[vecDeparture]
    nDepartures <- length(idsDeparture)

    if (nDepartures > 0) {
      active[idsDeparture] <- 0
      exitTime[idsDeparture] <- at
    }
  }

  ## Write out updated attributes ##
  dat <- set_attr(dat, "active", active)
  dat <- set_attr(dat, "exitTime", exitTime)

  ## Summary statistics ##
  dat <- set_epi(dat, "d.flow", at, nDepartures)

  return(dat)
}


# Arrival Module with Leaky Vaccination ------------------------------------

afunc <- function(dat, at) {
  # Simulate arrivals with leaky vaccination. Three vaccination routes:
  #
  # 1. INITIALIZATION (at == 2 only): The starting population may have
  #    pre-existing vaccination coverage. EpiModel convention: modules
  #    first run at timestep 2 (timestep 1 is reserved for init.net).
  #    Custom attributes not set in init.net are initialized here.
  #
  # 2. PROGRESSION: Each timestep, unvaccinated active nodes may receive
  #    vaccination at vaccination.rate.progression. Among newly vaccinated
  #    susceptibles, protection is conferred at protection.rate.progression.
  #
  # 3. ARRIVALS: New nodes enter as susceptible. Each may be vaccinated at
  #    vaccination.rate.arrivals, with protection at protection.rate.arrivals.
  #
  # In the leaky model, vaccine protection REDUCES the transmission
  # probability rather than preventing infection entirely. Protection
  # persists through the SEIRS cycle -- a protected individual who is
  # infected and later returns to S retains their protection.
  #
  # Key assumption: individuals may only be vaccinated once. Those who are
  # vaccinated but not protected (protection = "none") cannot be
  # re-vaccinated for another chance at protection.

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
  vaccination.rate.init <- get_param(dat, "vaccination.rate.initialization")
  protection.rate.init <- get_param(dat, "protection.rate.initialization")
  vaccination.rate.progression <- get_param(dat, "vaccination.rate.progression")
  protection.rate.progression <- get_param(dat, "protection.rate.progression")

  ## Initialize flow counters ##
  nVax.init <- 0
  nPrt.init <- 0
  nVax.prog <- 0
  nPrt.prog <- 0
  nVax.arrival <- 0
  nPrt.arrival <- 0

  ## --- INITIALIZATION (at == 2 only) --- ##
  # EpiModel convention: modules first run at timestep 2. Custom attributes
  # must be initialized here because they cannot be set in init.net().
  if (at == 2) {

    vaccination <- rep(NA, nTotal)
    protection <- rep(NA, nTotal)
    dat <- set_attr(dat, "vaccination", vaccination)
    dat <- set_attr(dat, "protection", protection)

    vaccination <- get_attr(dat, "vaccination")
    protection <- get_attr(dat, "protection")

    # Stochastic vaccination of initial population
    idsEligVacInit <- which(active == 1)
    vecVacInit <- rbinom(length(idsEligVacInit), 1, vaccination.rate.init)
    idsVacInit <- idsEligVacInit[which(vecVacInit == 1)]
    vaccination[idsVacInit] <- "initial"

    # Stochastic protection among vaccinated susceptibles
    idsEligProtInit <- which(vaccination == "initial" & status == "s")
    vecProtInit <- rbinom(length(idsEligProtInit), 1, protection.rate.init)
    idsProtInit <- idsEligProtInit[which(vecProtInit == 1)]
    idsNoProtInit <- setdiff(idsVacInit, idsProtInit)
    protection[idsProtInit] <- "initial"
    protection[idsNoProtInit] <- "none"

    nVax.init <- length(idsVacInit)
    nPrt.init <- length(idsProtInit)

    dat <- set_attr(dat, "vaccination", vaccination)
    dat <- set_attr(dat, "protection", protection)
  }

  vaccination <- get_attr(dat, "vaccination")
  protection <- get_attr(dat, "protection")

  ## --- VACCINATION PROGRESSION --- ##
  # Unvaccinated active nodes may receive vaccination each timestep
  idsEligVacProg <- which(is.na(vaccination) & active == 1)
  vecVacProg <- rbinom(length(idsEligVacProg), 1, vaccination.rate.progression)
  idsVacProg <- idsEligVacProg[which(vecVacProg == 1)]
  vaccination[idsVacProg] <- "progress"

  # Among newly vaccinated susceptibles, stochastic protection
  idsEligProtProg <- which(
    vaccination == "progress" & is.na(protection) & status == "s"
  )
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

  ## Write out updated attributes ##
  dat <- set_attr(dat, "vaccination", vaccination)
  dat <- set_attr(dat, "protection", protection)
  dat <- set_attr(dat, "infTime", infTime)
  dat <- set_attr(dat, "status", status)

  ## Summary statistics ##
  # Re-retrieve active after append_core_attr so lengths match
  active <- get_attr(dat, "active")
  dat <- set_epi(dat, "a.flow", at, nArrivals)
  dat <- set_epi(dat, "vac.flow", at, nVax.init + nVax.prog + nVax.arrival)
  dat <- set_epi(dat, "prt.flow", at, nPrt.init + nPrt.prog + nPrt.arrival)
  dat <- set_epi(dat, "vac.num", at,
                 sum(active == 1 & vaccination %in%
                       c("initial", "progress", "arrival"), na.rm = TRUE))
  dat <- set_epi(dat, "prt.num", at,
                 sum(active == 1 & protection %in%
                       c("initial", "progress", "arrival"), na.rm = TRUE))

  return(dat)
}
