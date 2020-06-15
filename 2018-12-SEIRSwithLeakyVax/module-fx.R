##
## SEIRS Model with Vital Dynamics and a Leaky Vaccine Implementation
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##
## Authors: Connor M. Van Meter
## Date: January 2018
##


# Replacement infection/transmission module -------------------------------

infect <- function(dat, at) {

  ## Uncomment this to function environment interactively
  #browser()

  ## Attributes ##
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  infTime <- get_attr(dat, "infTime")
  protection <- get_attr(dat, "protection")

  ## Parameters ##
  inf.prob <- get_param(dat, "inf.prob")
  act.rate <- get_param(dat,"act.rate")
  vaccine.efficacy <- get_param(dat, "vaccine.efficacy")

  ## Find infected nodes ##
  idsInf <- which(active == 1 & status == "i")
  nActive <- sum(active == 1)
  nElig <- length(idsInf)

  ## Initialize default incidence at 0 ##
  nInf <- 0

  ## If any infected nodes, proceed with transmission ##
  if (nElig > 0 && nElig < nActive) {

    ## Look up discordant edgelist ##
    del <- discord_edgelist(dat, at)

    ## If any discordant pairs, proceed ##
    if (!(is.null(del))) {

      # Set parameters on discordant edgelist data frame
      # If susceptible individuals in discordant edgelist
      # are vaccine protected, then they have reduced transmission probability
      inf.prob.reducedFOI <- (1 - vaccine.efficacy) * inf.prob
      idsReducedFOI <- which(!is.na(protection) & protection != "none" &
                               status == 's' & active == 1)
      reducedFOItransProbDF <-
        data.frame(idsReducedFOI,rep(inf.prob.reducedFOI,length(idsReducedFOI)))
      colnames(reducedFOItransProbDF) <- c("sus","transProb")
      del <- merge(del, reducedFOItransProbDF, by = "sus", all.x = TRUE)
      del$transProb <- ifelse(is.na(del$transProb),inf.prob,del$transProb)
      del$actRate <- act.rate
      del$finalProb <- 1 - (1 - del$transProb)^del$actRate

      # Stochastic transmission process
      transmit <- rbinom(nrow(del), 1, del$finalProb)

      # Keep rows where transmission occurred
      del <- del[which(transmit == 1), ]

      # Look up new ids if any transmissions occurred
      idsNewInf <- unique(del$sus)
      nInf <- length(idsNewInf)

      # Set new attributes for those newly infected
      if (nInf > 0) {
        status[idsNewInf] <- "e"
        dat <- set_attr(dat, "status", status)
        infTime[idsNewInf] <- at
        dat <- set_attr(dat, "infTime", infTime)
      }
    }
  }

  ## Save summary statistic for S->E flow
  dat <- set_epi(dat,"se.flow", at, nInf)

  #Vaccine Protected (Active) Number -
  #equivalent to dat$epi$protection.num.active[at]
  dat <- set_epi(dat, "s.num", at, sum(dat$attr$active == 1 & dat$attr$status == "s"))
  dat <- set_epi(dat, "e.num", at, sum(dat$attr$active == 1 & dat$attr$status == "e"))
  dat <- set_epi(dat, "i.num", at, sum(dat$attr$active == 1 & dat$attr$status == "i"))
  dat <- set_epi(dat, "r.num", at, sum(dat$attr$active == 1 & dat$attr$status == "r"))

  return(dat)
}


# New disease progression module ------------------------------------------
# (Replaces the recovery module)

progress <- function(dat, at) {

  ## Uncomment this to function environment interactively
  #browser()

  ## Attributes ##
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")

  ## Parameters ##
  ei.rate <- get_param(dat, "ei.rate")
  ir.rate <- get_param(dat, "ir.rate")
  rs.rate <- get_param(dat, "rs.rate")

  ## E to I progression process ##
  nInf <- 0
  idsEligInf <- which(active == 1 & status == "e")
  nEligInf <- length(idsEligInf)

  if (nEligInf > 0) {
    vecInf <- which(rbinom(nEligInf, 1, ei.rate) == 1)
    if (length(vecInf) > 0) {
      idsInf <- idsEligInf[vecInf]
      nInf <- length(idsInf)
      status[idsInf] <- "i"
    }
  }

  ## I to R progression process ##
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

  # ## R to S progression process ##
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
  ## Write out updated status attribute ##
  dat <- set_attr(dat, "status", status)

  ## Save summary statistics ##
  dat <- set_epi(dat, "ei.flow", at, nInf)
  dat <- set_epi(dat, "ir.flow", at, nRec)
  dat <- set_epi(dat, "rs.flow", at, nSus)

  return(dat)
}

# Update Departure Module -----------------------------------------------------

dfunc <- function(dat, at) {
  #browser()

  ## Attributes ##
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")

  ## Parameters ##
  departure.rate <- get_param(dat, "departure.rate")
  departure.rates <- rep(departure.rate, network.size(dat$nw[[1]]))
  departure.dis.mult <- get_param(dat, "departure.disease.mult")

  ## Query alive ##
  idsElig <- which(active == 1)
  nElig <- length(idsElig)
  nDepartures <- 0

  if (nElig > 0) {

    departure_rates_of_elig <- departure.rates[idsElig]

    ## Multiply departure rates for diseased persons
    idsElig.inf <- which(status[idsElig] == "i")
    departure_rates_of_elig[idsElig.inf] <- departure.rates[idsElig.inf] *
      departure.dis.mult

    ## Simulate departure process
    vecDeparture <- which(rbinom(nElig, 1, departure_rates_of_elig) == 1)
    idsDeparture <- idsElig[vecDeparture]
    nDepartures <- length(idsDeparture)

    ## Update nodal attributes on attr and networkDynamic object ##
    if (nDepartures > 0) {
      active[idsDeparture] <- 0
      dat$nw[[1]] <- deactivate.vertices(dat$nw[[1]], onset = at, terminus = Inf,
                                         v = idsDeparture, deactivate.edges = TRUE)
    }

    ## Write out updated status attribute ##
    dat <- set_attr(dat, "active", active)
  }

  ## Summary statistics ##
  dat <- set_epi(dat,"d.flow", at, nDepartures)

  return(dat)
}


# Updated Arrival Module ----------------------------------------------------

afunc <- function(dat, at) {

  #Toggle for step-through debugging
  #browser()

  ## Attributes ##
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  if (at == 2) {
    infTime <- rep(NA, length(active))
    infTime[which(status == "i")] <- 1
    dat <- set_attr(dat, "infTime", infTime)
  } else {
    infTime <- get_attr(dat, "infTime")
  }
  entrTime <- get_attr(dat, "entrTime")
  exitTime <- get_attr(dat, "exitTime")

  ## Parameters ##
  n <- network.size(dat$nw[[1]])
  a.rate <- get_param(dat, "arrival.rate")
  vaccination.rate.arrivals <- get_param(dat, "vaccination.rate.arrivals")
  protection.rate.arrivals <- get_param(dat, "protection.rate.arrivals")
  vaccination.rate.initialization <- get_param(dat, "vaccination.rate.initialization")
  protection.rate.initialization <- get_param(dat, "protection.rate.initialization")
  vaccination.rate.progression <- get_param(dat, "vaccination.rate.progression")
  protection.rate.progression <- get_param(dat,"protection.rate.progression")

  ## Initializing Vaccination and Protection Process Flow Count Variables ##
  nVax.init <- 0
  nPrt.init <- 0
  nVax.prog <- 0
  nPrt.prog <- 0
  nVax.arrival <- 0
  nPrt.arrival <- 0

  ## INITIALIZATION OF VACCINATION AND PROTECTION VERTEX (NODE) ATTRIBUTES ##
  if (at == 2) {

    #Initialize vaccination and protection vectors
    vaccination <- rep(NA, n)
    protection <- rep(NA,n)
    dat <- set_attr(dat, "vaccination", rep(NA, n))
    dat <- set_attr(dat, "protection", rep(NA, n))

    vaccination <- get_attr(dat, "vaccination")
    protection <- get_attr(dat, "protection")

    # Determine individuals at time t=2 who are initially vaccinated
    idsEligVacInit <- which(active == 1)
    vecVacInit <- rbinom(length(idsEligVacInit), 1,
                         vaccination.rate.initialization)
    idsVacInit <- idsEligVacInit[which(vecVacInit == 1)]
    vaccination[idsVacInit] <- "initial"

    #Determines if individual is protected based on
    #protection rate and infectious status
    idsEligProtInit <- which(vaccination == "initial" & status == 's')
    vecProtInit <- rbinom(length(idsEligProtInit), 1,
                          protection.rate.initialization)
    idsProtInit <- idsEligProtInit[which(vecProtInit == 1)]
    idsNoProtInit <- setdiff(idsVacInit, idsProtInit)
    protection[idsProtInit] <- "initial"
    protection[idsNoProtInit] <- "none"

    #Captures the number of vaccinated and the number of protected (active)
    #individuals at the time of initialization
    nVax.init <- length(idsVacInit)
    nPrt.init <- length(idsProtInit)

    #Update attribute list
    dat <- set_attr(dat, "vaccination", vaccination)
    dat <- set_attr(dat, "protection", protection)
  }


  vaccination <- get_attr(dat, "vaccination")
  protection <- get_attr(dat, "protection")


  ## VACCINATION PROGRESSION PROCESSES ##

  #Update the vaccination vector through the vaccination progression process
  idsEligVacProg <- which(is.na(vaccination) & active == 1)
  vecVacProg <- rbinom(length(idsEligVacProg), 1, vaccination.rate.progression)
  idsVacProg <- idsEligVacProg[which(vecVacProg == 1)]
  vaccination[idsVacProg] <- "progress"

  #Update the protection vector through the vaccination protection progression
  #process
  idsEligProtProg <- which(vaccination == "progress" &  is.na(protection)
                           & status == 's')
  vecProtProg <- rbinom(length(idsEligProtProg), 1, protection.rate.progression)
  idsProtProg <- idsEligProtProg[which(vecProtProg == 1)]
  idsNoProtProg <- setdiff(idsVacProg, idsProtProg)
  protection[idsProtProg] <- "progress"
  protection[idsNoProtProg] <- "none"

  #Captures the total number of vaccinated and protected individuals
  #after running vaccination and protection processes for current time step
  nVax.prog <- length(idsVacProg)
  nPrt.prog <- length(idsProtProg)


  ## ARRIVALS AND ARRIVAL VACCINATION PROCESSES ##

  #Arrival Process
  nArrivalsExp <- n * a.rate
  nArrivals <- rpois(1, nArrivalsExp)
  nVax.arrival <- 0
  nPrt.arrival <- 0

  if (nArrivals > 0) {
    dat$nw[[1]] <- add.vertices(dat$nw[[1]], nv = nArrivals)
    newNodes <- (n + 1):(n + nArrivals)
    dat$nw[[1]] <- activate.vertices(dat$nw[[1]], onset = at, terminus = Inf,
                                     v = newNodes)
  }

  #Update attributes
  if (nArrivals > 0) {
    active <- c(active, rep(1, nArrivals))
    status <- c(status, rep("s", nArrivals))
    infTime <- c(infTime, rep(NA, nArrivals))
    entrTime <- c(entrTime, rep(at, nArrivals))
    exitTime <- c(exitTime, rep(NA, nArrivals))
    vaccination <- c(vaccination, rep(NA, nArrivals))
    protection <- c(protection, rep(NA, nArrivals))

    # New arrival vaccination process and count
    vaccinatedNewArrivals <- rbinom(nArrivals, 1, vaccination.rate.arrivals)
    vaccination[newNodes] <- ifelse(vaccinatedNewArrivals == 1, "arrival", NA)
    nVax.arrival <- length(which(vaccinatedNewArrivals == 1))

    #Update the protection vector through the vaccination protection at arrival
    #process
    idsEligProtArrival <- which(vaccination == "arrival" & status == 's' &
                                  is.na(protection))
    vecProtArrival <- rbinom(length(idsEligProtArrival), 1,
                             protection.rate.arrivals)
    idsProtArrival <- idsEligProtArrival[which(vecProtArrival == 1)]
    idsNoProtArrival <- idsEligProtArrival[which(vecProtArrival == 0)]
    protection[idsProtArrival] <- "arrival"
    protection[idsNoProtArrival] <- "none"
    nPrt.arrival <- length(which(vecProtArrival == 1))

  }

  ## UPDATE NODE ATTRIBUTES ##
  dat <- set_attr(dat, "active", active, override.length.check = TRUE)
  dat <- set_attr(dat, "vaccination", vaccination)
  dat <- set_attr(dat, "protection", protection)
  dat <- set_attr(dat, "infTime", infTime)
  dat <- set_attr(dat, "entrTime", entrTime)
  dat <- set_attr(dat, "exitTime", exitTime)
  dat <- set_attr(dat, "status", status)

  ## SUMMARY STATISTICS ##

  #Arrivals
  dat <- set_epi(dat, "a.flow", at, nArrivals)

  #Vaccination and Protection
  dat <- set_epi(dat, "vac.flow", at, nVax.init + nVax.prog + nVax.arrival)
  dat <- set_epi(dat, "prt.flow", at, nPrt.init + nPrt.prog + nPrt.arrival)
  dat <- set_epi(dat, "vac.num", at, sum(active == 1 & vaccination %in%
                                           c("initial", "progress", "arrival")))
  dat <- set_epi(dat, "prt.num", at, sum(active == 1 & protection %in%
                                           c("initial", "progress", "arrival")))

  dat <- set_epi(dat, "vac.init.flow", at, nVax.init)
  dat <- set_epi(dat, "prt.init.flow", at, nPrt.init)
  dat <- set_epi(dat, "vac.prog.flow", at, nVax.prog)
  dat <- set_epi(dat, "prt.prog.flow", at, nPrt.prog)
  dat <- set_epi(dat, "vac.arrival.flow", at, nVax.arrival)
  dat <- set_epi(dat, "prt.arrival.flow", at, nPrt.arrival)

  dat <- set_epi(dat, "vac.init.num", at, sum(active == 1 & !is.na(vaccination) &
                                                vaccination == "initial"))
  dat <- set_epi(dat, "prt.init.num", at, sum(active == 1 & !is.na(protection) &
                                                protection == "initial"))
  dat <- set_epi(dat, "vac.prog.num", at, sum(active == 1 & !is.na(vaccination) &
                                                vaccination == "progress"))
  dat <- set_epi(dat, "prt.prog.num", at, sum(active == 1 & !is.na(protection) &
                                                protection == "progress"))
  dat <- set_epi(dat, "vac.arrival.num", at,  sum(active == 1 & !is.na(vaccination) &
                                                    vaccination == "arrival"))
  dat <- set_epi(dat, "prt.arrival.num", at, sum(active == 1 & !is.na(protection) &
                                                   protection == "arrival"))

  return(dat)
}
