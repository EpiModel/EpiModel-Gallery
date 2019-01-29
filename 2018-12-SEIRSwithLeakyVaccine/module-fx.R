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
  active <- dat$attr$active
  status <- dat$attr$status
  protection <- dat$attr$protection

  ## Parameters ##
  inf.prob <- dat$param$inf.prob
  act.rate <- dat$param$act.rate
  vaccine.efficacy <- dat$param$vaccine.efficacy

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
        dat$attr$status[idsNewInf] <- "e"
        dat$attr$infTime[idsNewInf] <- at
      }
    }
  }

  ## Save summary statistic for S->E flow
  dat$epi$se.flow[at] <- nInf

  #Vaccine Protected (Active) Number -
  #equivalent to dat$epi$protection.num.active[at]
  dat$epi$s.num[at] <- sum(dat$attr$active == 1 & dat$attr$status == "s")
  dat$epi$e.num[at] <- sum(dat$attr$active == 1 & dat$attr$status == "e")
  dat$epi$i.num[at] <- sum(dat$attr$active == 1 & dat$attr$status == "i")
  dat$epi$r.num[at] <- sum(dat$attr$active == 1 & dat$attr$status == "r")

  return(dat)
}


# New disease progression module ------------------------------------------
# (Replaces the recovery module)

progress <- function(dat, at) {

  ## Uncomment this to function environment interactively
  #browser()

  ## Attributes ##
  active <- dat$attr$active
  status <- dat$attr$status

  ## Parameters ##
  ei.rate <- dat$param$ei.rate
  ir.rate <- dat$param$ir.rate
  rs.rate <- dat$param$rs.rate

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
  dat$attr$status <- status

  ## Save summary statistics ##
  dat$epi$ei.flow[at] <- nInf
  dat$epi$ir.flow[at] <- nRec
  dat$epi$rs.flow[at] <- nSus

  return(dat)
}

# Update Departure Module -----------------------------------------------------

dfunc <- function(dat, at) {
  #browser()

  ## Attributes ##
  active <- dat$attr$active
  status <- dat$attr$status

  ## Parameters ##
  departure.rates <- rep(dat$param$departure.rate, network.size(dat$nw))
  departure.dis.mult <- dat$param$departure.disease.mult

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
      dat$nw <- deactivate.vertices(dat$nw, onset = at, terminus = Inf,
                                    v = idsDeparture, deactivate.edges = TRUE)
    }

    ## Write out updated status attribute ##
    dat$attr$active <- active
  }

  ## Summary statistics ##
  dat$epi$d.flow[at] <- nDepartures
  if (at == 2) {
    dat$epi$d.num[at] <- sum(active == 0)
  } else {
    dat$epi$d.num[at] <- dat$epi$d.num[at - 1] + sum(active == 0)
  }

  return(dat)
}


# Updated Arrival Module ----------------------------------------------------

afunc <- function(dat, at) {

  #Toggle for step-through debugging
  #browser()

  ## Attributes ##
  active <- dat$attr$active
  status <- dat$attr$status
  infTime <- dat$attr$infTime
  entrTime <- dat$attr$entrTime
  exitTime <- dat$attr$exitTime
  vaccination <- dat$attr$vaccination
  protection <- dat$attr$protection

  ## Parameters ##
  n <- network.size(dat$nw)
  a.rate <- dat$param$arrival.rate
  vaccination.rate.arrivals <- dat$param$vaccination.rate.arrivals
  protection.rate.arrivals <- dat$param$protection.rate.arrivals
  vaccination.rate.initialization <- dat$param$vaccination.rate.initialization
  protection.rate.initialization <- dat$param$protection.rate.initialization
  vaccination.rate.progression <- dat$param$vaccination.rate.progression
  protection.rate.progression <- dat$param$protection.rate.progression

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
    vaccination <- rep(NA,n)
    protection <- rep(NA,n)

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

  }


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
    dat$nw <- add.vertices(dat$nw, nv = nArrivals)
    newNodes <- (n + 1):(n + nArrivals)
    dat$nw <- activate.vertices(dat$nw, onset = at, terminus = Inf,
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
  dat$attr$vaccination <- vaccination
  dat$attr$protection <- protection
  dat$attr$active <- active
  dat$attr$infTime <- infTime
  dat$attr$entrTime <- entrTime
  dat$attr$exitTime <- exitTime
  dat$attr$status <- status


  ## SUMMARY STATISTICS ##

  #Arrivals
  dat$epi$a.flow[at] <- nArrivals
  if (at == 2) {
    dat$epi$a.num[at] <- nArrivals
  } else {
    dat$epi$a.num[at] <- dat$epi$a.num[at - 1] + nArrivals
  }

  #Vaccination and Protection
  dat$epi$vac.flow[at] <- nVax.init + nVax.prog + nVax.arrival
  dat$epi$prt.flow[at] <- nPrt.init + nPrt.prog + nPrt.arrival
  dat$epi$vac.num[at] <- sum(active == 1 & vaccination %in%
                               c("initial", "progress", "arrival"))
  dat$epi$prt.num[at] <- sum(active == 1 & protection %in%
                               c("initial", "progress", "arrival"))

  dat$epi$vac.init.flow[at] <- nVax.init
  dat$epi$prt.init.flow[at] <- nPrt.init
  dat$epi$vac.prog.flow[at] <- nVax.prog
  dat$epi$prt.prog.flow[at] <- nPrt.prog
  dat$epi$vac.arrival.flow[at] <- nVax.arrival
  dat$epi$prt.arrival.flow[at] <- nPrt.arrival

  dat$epi$vac.init.num[at] <- sum(active == 1 & !is.na(vaccination) &
                                    vaccination == "initial")
  dat$epi$prt.init.num[at] <- sum(active == 1 & !is.na(protection) &
                                    protection == "initial")
  dat$epi$vac.prog.num[at] <- sum(active == 1 & !is.na(vaccination) &
                                    vaccination == "progress")
  dat$epi$prt.prog.num[at] <- sum(active == 1 & !is.na(protection) &
                                    protection == "progress")
  dat$epi$vac.arrival.num[at] <- sum(active == 1 & !is.na(vaccination) &
                                       vaccination == "arrival")
  dat$epi$prt.arrival.num[at] <- sum(active == 1 & !is.na(protection) &
                                       protection == "arrival")

  return(dat)
}
