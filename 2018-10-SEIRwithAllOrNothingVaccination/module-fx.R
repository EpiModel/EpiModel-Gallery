##
## SEIR Model with Vital Dynamics and an All or Nothing Vaccine Implementation
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##
## Authors: Connor M. Van Meter
## Date: November 2018
##


# Replacement infection/transmission module -------------------------------

infect <- function(dat, at) {

  ## Uncomment this to function environment interactively
  #browser()

  ## Attributes ##
  active <- dat$attr$active
  status <- dat$attr$status

  ## Parameters ##
  inf.prob <- dat$param$inf.prob
  act.rate <- dat$param$act.rate

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
      del$transProb <- inf.prob
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
  dat$epi$v.num[at] <- sum(dat$attr$active == 1 & dat$attr$status == "v")

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

  ## Write out updated status attribute ##
  dat$attr$status <- status

  ## Save summary statistics ##
  dat$epi$ei.flow[at] <- nInf
  dat$epi$ir.flow[at] <- nRec

  return(dat)
}

# Update Death Module -----------------------------------------------------

dfunc <- function(dat, at) {

  #browser()

  ## Attributes ##
  active <- dat$attr$active
  status <- dat$attr$status

  ## Parameters ##
  mort.rates <- rep(dat$param$mortality.rate, network.size(dat$nw))
  mort.dis.mult <- dat$param$mortality.disease.mult

  ## Query alive ##
  idsElig <- which(active == 1)
  nElig <- length(idsElig)
  nDeaths <- 0

  if (nElig > 0) {

    death_rates_of_elig <- mort.rates[idsElig]

    ## Multiply death rates for diseased persons
    idsElig.inf <- which(status[idsElig] == "i")
    death_rates_of_elig[idsElig.inf] <- mort.rates[idsElig.inf] * mort.dis.mult

    ## Simulate mortality process
    vecDeaths <- which(rbinom(nElig, 1, death_rates_of_elig) == 1)
    idsDeaths <- idsElig[vecDeaths]
    nDeaths <- length(idsDeaths)

    ## Update nodal attributes on attr and networkDynamic object ##
    if (nDeaths > 0) {
      active[idsDeaths] <- 0
      dat$nw <- deactivate.vertices(dat$nw, onset = at, terminus = Inf,
                                    v = idsDeaths, deactivate.edges = TRUE)
    }

    ## Write out updated status attribute ##
    dat$attr$active <- active
  }

  ## Summary statistics ##
  dat$epi$d.flow[at] <- nDeaths
  if (at == 2) {
    dat$epi$d.num[at] <- sum(active == 0)
  } else {
    dat$epi$d.num[at] <- dat$epi$d.num[at - 1] + sum(active == 0)
  }

  return(dat)
}


# Updated Birth Module ----------------------------------------------------

bfunc <- function(dat, at) {

  #Toggle for step-through debugging
  #browser()

  ## Attributes ##
  active <- dat$attr$active
  status <- dat$attr$status
  infTime <- dat$attr$infTime
  entrTime <- dat$attr$entrTime
  exitTime <- dat$attr$exitTime
  vaccination <- dat$attr$vaccination
  vaccination_method <- dat$attr$vaccination_method
  protection <- dat$attr$protection
  protection_method <- dat$attr$protection_method
  infTime <- dat$attr$infTime
  entrTime <- dat$attr$entrTime
  exitTime <- dat$attr$exitTime

  ## Parameters ##
  n <- network.size(dat$nw)
  b.rate <- dat$param$birth.rate
  vaccination.rate.births <- dat$param$vaccination.rate.births
  protection.rate.births <- dat$param$protection.rate.births
  vaccination.rate.initialization <- dat$param$vaccination.rate.initialization
  protection.rate.initialization <- dat$param$protection.rate.initialization
  vaccination.rate.progression <- dat$param$vaccination.rate.progression
  protection.rate.progression <- dat$param$protection.rate.progression

  ## Initializing Vaccination and Protection Process Flow Count Variables ##
  nVax.initialization <- 0
  nPrt.initialization <- 0
  nVax.progression <- 0
  nPrt.progression <- 0
  nVax.birth <- 0
  nPrt.birth <- 0

  ## INITIALIZATION OF VACCINATION AND PROTECTION VERTEX (NODE) ATTRIBUTES ##
  if (at == 2) {

    # Determine individuals at time t=2 who are initially vaccinated
    vaccination <- rbinom(n, 1, vaccination.rate.initialization)
    vaccination_method <- ifelse(vaccination == 0, NA, "initial")

    #Assuption for the All-Or-Nothing Vaccine Model is that Infected/Infectious
    #Individuals Cannot Be Protected at Start of Simulation
    #Create a Vector, infectionStatusAtInitializationVector, in which status of
    #"i" = 0 and not "i" = 1
    infectionStatusAtInitializationVector <- ifelse(status == "i", 0, 1)

    #Determines if individual is protected based on
    #protection rate and infectious status
    protection <- vaccination * rbinom(n, 1, protection.rate.initialization) *
      infectionStatusAtInitializationVector
    protection_method <- ifelse(protection == 0, NA, "initial")

    #Captures the number of vaccinated and the number of protected (active)
    #individuals at the time of initialization
    nVax.initialization <- length(which(vaccination == 1))
    nPrt.initialization <- length(which(protection == 1))

  }

  #Create a vector of 0s and 1s, where 1s are unvaccinated active individuals
  unvaccinated <- ifelse(vaccination == 0 & active == 1, 1, 0)

  #First, create a vaccinationUnvaccinated vector of 0s and 1s,
  #where 1s are newly vaccinated individuals who were active and unvaccinated
  #at start of time period.
  #Next, update the vaccination vector, keeping 1s for those
  #previously vaccinated through the vaccination progression process
  #and then updating remaining 0s based on vaccinationUnvaccinated.
  vaccinationUnvaccinated <- unvaccinated *
    rbinom(n, 1, vaccination.rate.progression)
  vaccination <- ifelse(vaccination == 1, 1, vaccinationUnvaccinated)

  #Set vaccination method to 'progress' for new vaccinations
  noVaccination <- vaccination == 0 & is.na(vaccination_method)
  prevVaccination <- vaccination == 1 & !is.na(vaccination_method)
  vaccination_method <- ifelse(noVaccination | prevVaccination,
                               vaccination_method, "progress")

  #First, create a protectionUnvaccinatedSusceptibles vector of 0s and 1s,
  #where 1s are newly protected individuals who were vaccinated with
  #the vaccinationUnvaccinated process, are susceptible,
  #and confer vaccine immunity based on the user-input
  #vaccination protection rate for the progression process.
  #Next, update the protection progression vector, keeping 1s for those
  #previously protected through the protection progression process
  #and then updating remaining 0s based on protectionUnvaccinatedSusceptibles.
  protectionUnvaccinatedSusceptibles <- vaccinationUnvaccinated *
    rbinom(n, 1, protection.rate.progression) * ifelse(status == "s", 1, 0)
  protection <- ifelse(protection == 1, 1, protectionUnvaccinatedSusceptibles)
  noProtection <- protection == 0 & is.na(protection_method)
  prevProtection <- protection == 1 & !is.na(protection_method)
  protection_method <- ifelse(noProtection | prevProtection,
                              protection_method, "progress")

  #Captures the total number of vaccinated and protected individuals
  #after running vaccination and protection processes for current time step
  nVax.progression <- length(which(vaccinationUnvaccinated == 1))
  nPrt.progression <- length(which(protectionUnvaccinatedSusceptibles == 1))


  ## BIRTH AND BIRTH VACCINATION PROCESSES ##

  #Birth Process
  nBirthsExp <- n * b.rate
  nBirths <- rpois(1, nBirthsExp)
  nvaccination.births <- 0
  nprotection.births <- 0

  if (nBirths > 0) {
    dat$nw <- add.vertices(dat$nw, nv = nBirths)
    newNodes <- (n + 1):(n + nBirths)
    dat$nw <- activate.vertices(dat$nw, onset = at, terminus = Inf, v = newNodes)
  }

  #Update attributes (pre-existing vectors)
  #Continue with vaccination and protection of births
  if (nBirths > 0) {
    active <- c(active, rep(1, nBirths))
    status <- c(status, rep("s", nBirths))
    infTime <- c(infTime, rep(NA, nBirths))
    entrTime <- c(entrTime, rep(at, nBirths))
    exitTime <- c(exitTime, rep(NA, nBirths))

    # New birth vaccination process and count
    vaccinatedNewBirths <- rbinom(nBirths, 1, vaccination.rate.births)
    nVax.birth <- length(which(vaccinatedNewBirths == 1))

    # New birth vaccination protection process and count
    protectedNewBirths <- vaccinatedNewBirths * rbinom(nBirths, 1,
                                                       protection.rate.births)
    nPrt.birth <- length(which(protectedNewBirths == 1))

    # Update vaccination and protection vectors
    vaccination <- c(vaccination, vaccinatedNewBirths)
    protection <- c(protection, protectedNewBirths)
    vaccination_method <- c(vaccination_method, rep(NA, nBirths))
    protection_method <- c(protection_method, rep(NA, nBirths))
    noVaccination <- vaccination == 0 & is.na(vaccination_method)
    prevVaccination <- vaccination == 1 & !is.na(vaccination_method)
    vaccination_method <- ifelse(noVaccination | prevVaccination,
                                 vaccination_method, "birth")
    noProtection <- protection == 0 & is.na(protection_method)
    prevProtection <- protection == 1 & !is.na(protection_method)
    protection_method <- ifelse(noProtection | prevProtection,
                                protection_method, "birth")
  }

  ## UPDATE NODE ATTRIBUTES ##

  dat$attr$status <- ifelse(status == "s" & protection == 1
                            & active == 1, "v", status)
  dat$attr$vaccination <- vaccination
  dat$attr$protection <- protection
  dat$attr$vaccination_method <- vaccination_method
  dat$attr$protection_method <- protection_method
  dat$attr$active <- active
  dat$attr$infTime <- infTime
  dat$attr$entrTime <- entrTime
  dat$attr$exitTime <- exitTime


  ## SUMMARY STATISTICS ##

  #Births
  dat$epi$b.flow[at] <- nBirths
  if (at == 2) {
    dat$epi$b.num[at] <- nBirths
  } else {
    dat$epi$b.num[at] <- dat$epi$b.num[at - 1] + nBirths
  }

  #Vaccination and Protection
  dat$epi$vac.flow[at] <- nVax.initialization + nVax.progression + nVax.birth
  dat$epi$prt.flow[at] <- nPrt.initialization + nPrt.progression + nPrt.birth
  dat$epi$vac.num[at] <- sum(active == 1 & vaccination == 1)
  dat$epi$prt.num[at] <- sum(active == 1 & protection == 1)

  dat$epi$vac.init.flow[at] <- nVax.initialization
  dat$epi$prt.init.flow[at] <- nPrt.initialization
  dat$epi$vac.prog.flow[at] <- nVax.progression
  dat$epi$prt.prog.flow[at] <- nPrt.progression
  dat$epi$vac.birth.flow[at] <- nVax.birth
  dat$epi$prt.birth.flow[at] <- nPrt.birth

  dat$epi$vac.init.num[at] <- sum(active == 1 & !is.na(vaccination_method)
                                  & vaccination_method == "initial")
  dat$epi$prt.init.num[at] <- sum(active == 1 & !is.na(protection_method)
                                  & protection_method == "initial")
  dat$epi$vac.prog.num[at] <- sum(active == 1 & !is.na(vaccination_method)
                                  & vaccination_method == "progress")
  dat$epi$prt.prog.num[at] <- sum(active == 1 & !is.na(protection_method)
                                  & protection_method == "progress")
  dat$epi$vac.birth.num[at] <- sum(active == 1 & !is.na(vaccination_method)
                                   & vaccination_method == "birth")
  dat$epi$prt.birth.num[at] <- sum(active == 1 & !is.na(protection_method)
                                   & protection_method == "birth")

  return(dat)
}
