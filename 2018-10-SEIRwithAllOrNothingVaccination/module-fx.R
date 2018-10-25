##
## SEIR Model with Vital Dynamics and an All or Nothing vaccination Implementation
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##
## Authors: Samuel M. Jenness, Venkata R. Duvvuri, Connor Van Meter
## Date: October 2018
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
  dat$epi$e.num[at] <- sum(dat$attr$active == 1 & dat$attr$status == "e")

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
  dat$epi$s.num[at] <- sum(dat$attr$active == 1 & dat$attr$status == "s")
  dat$epi$i.num[at] <- sum(dat$attr$active == 1 & dat$attr$status == "i")
  dat$epi$r.num[at] <- sum(dat$attr$active == 1 & dat$attr$status == "r")

  return(dat)
}

# Update Death Module -----------------------------------------------------

dfunc <- function(dat, at) {

  #browser()

  ## Attributes ##
  active <- dat$attr$active
  status <- dat$attr$status

  ## Parameters ##
  mort.rates <- rep(dat$param$mortality.rate,network.size(dat$nw))
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
      dat$attr$active[idsDeaths] <- 0
      dat$nw <- deactivate.vertices(dat$nw, onset = at, terminus = Inf,
                                    v = idsDeaths, deactivate.edges = TRUE)
    }
  }

  ## Summary statistics ##
  dat$epi$d.flow[at] <- nDeaths

  return(dat)
}


# Updated Birth Module ----------------------------------------------------

bfunc <- function(dat, at) {

  #Toggle for step-through debugging
  #browser()


  ## Attributes ##
  active <- dat$attr$active
  status <- dat$attr$status


  ## Parameters ##
  n <- network.size(dat$nw)
  b.rate <- dat$param$birth.rate
  nVaccination.initialization = 0
  nProtection.initialization = 0
  nVaccinatedNewBirths = 0
  nProtectedNewBirths = 0

  vaccination.rate <- dat$param$vaccination.rate #Uniform vaccination rate. If vaccination.rate.births, vaccination.rate.initialization, vaccination.rate.progression are specified, they take priority.
  protection.rate <- dat$param$protection.rate #Uniform protection rate. If protection.rate.births, protection.rate.initialization, protection.rate.progression are specified, they take priority.

  vaccination.rate.births <- ifelse(is.na(dat$param$vaccination.rate.births),dat$param$vaccination.rate,dat$param$vaccination.rate.births) #The vaccination rate in new births
  protection.rate.births <- ifelse(is.na(dat$param$protection.rate.births),dat$param$protection.rate,dat$param$protection.rate.births) #The vaccination protection rate in new births

  vaccination.rate.initialization <- ifelse(is.na(dat$param$vaccination.rate.initialization),dat$param$vaccination.rate,dat$param$vaccination.rate.initialization) #The vaccination rate in the initial network population
  protection.rate.initialization <- ifelse(is.na(dat$param$protection.rate.initialization),dat$param$protection.rate,dat$param$protection.rate.initialization) #The vaccination protection rate in the initial network population

  vaccination.rate.progression <- ifelse(is.na(dat$param$vaccination.rate.progression),dat$param$vaccination.rate,dat$param$vaccination.rate.progression) #The vaccination rate in susceptibles over time
  protection.rate.progression <- ifelse(is.na(dat$param$protection.rate.progression),dat$param$protection.rate,dat$param$protection.rate.progression) #The vaccination protection rate in susceptibles over time

  ## Check if any of the vaccination or protection rates are missing for the model and notify user if any rates are missing
  #Check if any vaccination rates are missing
  if (is.na(vaccination.rate.births)|is.na(vaccination.rate.initialization)|is.na(vaccination.rate.progression)) {
    stop("Specify either vaccination.rate OR (vaccination.rate.births AND vaccination.rate.initialization AND vaccination.rate.progression)")
  }
  #Check if any vaccination protection rates are missing
  if (is.na(protection.rate.births)|is.na(protection.rate.initialization)|is.na(protection.rate.progression)) {
    stop("Specify either protection.rate OR (protection.rate.births AND protection.rate.initialization AND protection.rate.progression)")
  }



  ## INITIALIZATION OF VACCINATION AND PROTECTION VERTEX (NODE) ATTRIBUTES ##
  if (at == 2) {

    #Initialize vector of zeros representing no "births" at initialization
    dat$attr$birthed <- rep(0,n)

    # Pull vaccination and protection attributes from the fitted network model
    vaccination.initialization <- rbinom(n, 1, vaccination.rate.initialization)

    #Assuption for the All-Or-Nothing Vaccine Model is that Infected/Infectious Individuals Cannot Be Protected at Start of Simulation
    #Create a Vector, infectionStatusAtInitializationVector, in which status of "i" = 0 and not "i" = 1
    infectionStatusAtInitializationVector <- rep(1,n)
    infInitVec <- which(status == "i")
    if(length(infInitVec) > 1){
      infectionStatusAtInitializationVector[infInitVec] <- 0
    }

    #Determines if individual is protection based on protection rate and infectious status
    protection.initialization <- vaccination.initialization * rbinom(n, 1, protection.rate.initialization) * infectionStatusAtInitializationVector

    #Captures vaccination and protection at initialization status and stores in new attributes
    dat$attr$vaccination.initialization <- vaccination.initialization
    dat$attr$protection.initialization <- protection.initialization
    dat$attr$vaccination.births <- rep(0,n)
    dat$attr$protection.births <- rep(0,n)
    dat$attr$vaccination.progression <- rep(0,n)
    dat$attr$protection.progression <- rep(0,n)
    dat$attr$vaccination <- ifelse(dat$attr$vaccination.initialization==1,1,0)
    dat$attr$protection <- ifelse(dat$attr$protection.initialization==1,1,0)

    #Captures the number of vaccinated and the number of protected (active) individuals at the time of initialization
    nVaccination.initialization <- length(which(vaccination.initialization == 1))
    nProtection.initialization <- length(which(protection.initialization == 1))

  }



  ## VACCINATION AND PROTECTION OF UNVACCINATED (ACTIVE) SUSCEPTIBLES PROCESS ##
  nUnvaccinatedSusceptibles <- length(which(dat$attr$status == "s" & dat$attr$vaccination == 0  & active==1))

  #Create a vector of 0s and 1s, where 1s are unvaccinated and susceptible (active) individuals
  unvaccinatedSusceptibles = ifelse(dat$attr$status == "s" & dat$attr$vaccination == 0 & active==1,1,0)

  #First, create a vaccinationUnvaccinatedSusceptibles vector of 0s and 1s,
  #where 1s are newly vaccinated individuals who were active, susceptible, and unvaccinated at start of time period
  #Next, update the vaccination progression vector, keeping 1s for those
  #previously vaccinated through the vaccination progression process and then updating remaining 0s based on vaccinationUnvaccinatedSusceptibles
  vaccinationUnvaccinatedSusceptibles <- unvaccinatedSusceptibles * rbinom(n, 1, vaccination.rate.progression)
  vaccination.progression <- ifelse(dat$attr$vaccination.progression==1,1,vaccinationUnvaccinatedSusceptibles)
  dat$attr$vaccination.progression <- vaccination.progression

  #First, create a protectionUnvaccinatedSusceptibles vector of 0s and 1s,
  #where 1s are newly protected individuals who were active, susceptible, and unvaccinated at start of time period,
  #but who then became vaccinated with the vaccinationUnvaccinatedSusceptibles step above.
  #Next, update the protection progression vector, keeping 1s for those
  #previously protected through the protection progression process and then updating remaining 0s based on protectionUnvaccinatedSusceptibles
  protectionUnvaccinatedSusceptibles <- vaccinationUnvaccinatedSusceptibles * rbinom(n, 1, protection.rate.progression)
  protection.progression <- ifelse(dat$attr$protection.progression==1,1,protectionUnvaccinatedSusceptibles)
  dat$attr$protection.progression <- protection.progression

  #Captures the total number of vaccinated and protected individuals after running vaccination and protection processes for current time step
  nvaccination.progression <- length(which(vaccination.progression == 1))
  nprotection.progression <- length(which(protection.progression == 1))

  #Captures the number of vaccinated and the number of protected ACTIVE individuals after running vaccination and protection processes for current time step
  nvaccination.progression.active <- length(which(vaccination.progression == 1 & active==1))
  nprotection.progression.active <- length(which(protection.progression == 1 & active==1))

  #Captures the number of newly vaccinated and protected (active) individuals
  #Note: active==1 is not specified in this count as active=1 is implied.
  #active=1 is a requirement for unvaccinatedSusceptibles=1, which is used to determine
  #if vaccinationUnvaccinatedSusceptibles=1 which is used to determine if protectionUnvaccinatedSusceptibles=1
  nvaccination.progression.new <- length(which(vaccinationUnvaccinatedSusceptibles == 1))
  nprotection.progression.new <- length(which(protectionUnvaccinatedSusceptibles == 1))

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

  # Update attributes (pre-existing vectors) and continue with vaccination and protection of births
  if (nBirths > 0) {
    dat$attr$active <- c(active, rep(1, nBirths))
    dat$attr$status <- c(status, rep("s", nBirths))
    dat$attr$infTime <- c(dat$attr$infTime, rep(NA, nBirths))
    dat$attr$entrTime <- c(dat$attr$entrTime, rep(at, nBirths))
    dat$attr$exitTime <- c(dat$attr$exitTime, rep(NA, nBirths))

    dat$attr$birthed <- c(dat$attr$birthed,rep(1,nBirths))

    dat$attr$vaccination.initialization <- c(dat$attr$vaccination.initialization, rep(0,nBirths))
    dat$attr$protection.initialization <- c(dat$attr$protection.initialization, rep(0,nBirths))

    dat$attr$vaccination.progression <- c(dat$attr$vaccination.progression, rep(0,nBirths))
    dat$attr$protection.progression <- c(dat$attr$protection.progression, rep(0,nBirths))

    # New birth vaccination process and count
    vaccinatedNewBirths <- rbinom(nBirths, 1, vaccination.rate.births)
    nVaccinatedNewBirths <- length(which(vaccinatedNewBirths==1))

    # New birth vaccination protection process and count
    protectedNewBirths <- vaccinatedNewBirths * rbinom(nBirths, 1, protection.rate.births)
    nProtectedNewBirths <- length(which(protectedNewBirths==1))

    #Append new birth vaccination and protection attribute vectors to existing birth vaccination and protection attribute vectors
    vaccination.births <- c(dat$attr$vaccination.births,vaccinatedNewBirths)
    protection.births <- c(dat$attr$vaccination.births,protectedNewBirths)
    dat$attr$vaccination.births <- vaccination.births
    dat$attr$protection.births <- protection.births

  }

  #Determining Overall Vaccination and Protection Vectors
  dat$attr$vaccination <- ifelse(dat$attr$vaccination.initialization == 1 | dat$attr$vaccination.progression == 1 | dat$attr$vaccination.births == 1,1,0)
  dat$attr$protection <- ifelse(dat$attr$protection.initialization == 1 | dat$attr$protection.progression == 1 | dat$attr$protection.births == 1,1,0)


  ## S (susceptible) to V (vaccine protected) progression process ##
  status <- ifelse(status=="s" & dat$attr$protection==1 & active==1,"v",status)

  ## Write out updated status attribute ##
  dat$attr$status <- status


  ## SUMMARY STATISTICS ##

  #Births
  dat$epi$b.flow[at] <- nBirths # number of new births
  dat$epi$b.num.active[at] <- sum(dat$attr$birthed == 1 & dat$attr$active == 1) # total number of births that are active
  dat$epi$b.num[at] <- sum(dat$attr$birthed == 1) # total number of births


  #Vaccine Protected (Active) Number - equivalent to dat$epi$protection.num.active[at]
  dat$epi$v.num[at] <- sum(dat$attr$active == 1 & dat$attr$status == "v")


  #Vaccination and Protection - Initialization
  dat$epi$vac.init.flow[at] <- nVaccination.initialization # number of newly vaccinated individuals from vaccination initialization process
  dat$epi$prot.init.flow[at] <- nProtection.initialization # number of newly vaccinated and protected individuals from protection initialization process
  dat$epi$vac.init.num.active[at] <- sum(dat$attr$active == 1 & dat$attr$vaccination.initialization == 1) # total number of active vaccinated individuals from vaccination initialization process
  dat$epi$prot.init.num.active[at] <- sum(dat$attr$active == 1 & dat$attr$protection.initialization == 1) # total number of active vaccinated and protected individuals from protection initialization process
  dat$epi$vac.init.num[at] <- sum(dat$attr$vaccination.initialization == 1) # total number of vaccinated individuals from vaccination initialization process
  dat$epi$prot.init.num[at] <- sum(dat$attr$protection.initialization == 1) # total number of vaccinated and protected individuals from protection initialization process


  #Vaccination and Protection - Progression
  dat$epi$unVacSus.num[at] <- nUnvaccinatedSusceptibles # Number of active, unvaccinated, susceptible individuals
  dat$epi$vac.prog.flow[at] <- nvaccination.progression.new # Number of active, unvaccinated, susceptible individuals who become vaccinated with this timestep
  dat$epi$prot.prog.flow[at] <- nprotection.progression.new # Number of active, unvaccinated, susceptible individuals who become vaccinated and protected with this timestep
  dat$epi$vac.prog.num.active[at] <- nvaccination.progression.active # Total number of active, unvaccinated, susceptible individuals who have become vaccinated through the vaccination progression process who are still active
  dat$epi$prot.prog.num.active[at] <- nprotection.progression.active # Total number of active, unvaccinated, susceptible individuals who have become vaccinated and protected through the protection progression process who are still active
  dat$epi$vac.prog.num[at] <- nvaccination.progression # Total number of active, unvaccinated, susceptible individuals who have become vaccinated through the vaccination progression process
  dat$epi$prot.prog.num[at] <- nprotection.progression # Total number of active, unvaccinated, susceptible individuals who have become vaccinated and protected through the protection progression process

  #Vaccination and Protection - Births
  dat$epi$vac.b.flow[at] <- nVaccinatedNewBirths #New births during this time period that are vaccinated
  dat$epi$prot.b.flow[at] <- nProtectedNewBirths #New births during this time period that are vaccinated and protected
  dat$epi$vac.b.num.active[at] <- sum(dat$attr$active == 1 & dat$attr$vaccination.births == 1) #Total number of active individuals who were "birthed" and are vaccinated
  dat$epi$prot.b.num.active[at] <- sum(dat$attr$active == 1 & dat$attr$protection.births == 1) #Total number of active individuals who were "birthed" and are vaccinated and protected
  dat$epi$vac.b.num[at] <- sum(dat$attr$vaccination.births == 1) #Total number of individuals who were "birthed" and are vaccinated
  dat$epi$prot.b.num[at] <- sum(dat$attr$protection.births == 1) #Total number of individuals who were "birthed" and are vaccinated and protected

  #Vaccination and Protection - Total
  dat$epi$vac.flow[at] <- dat$epi$vaccination.initialization.flow[at] + dat$epi$vaccination.progression.flow[at] + dat$epi$vaccination.births.flow[at]
  dat$epi$prot.flow[at] <- dat$epi$protection.initialization.flow[at] + dat$epi$protection.progression.flow[at] + dat$epi$protection.births.flow[at]
  dat$epi$vac.num.active[at] <- sum(dat$attr$active == 1 & dat$attr$vaccination == 1)
  dat$epi$prot.num.active[at] <- sum(dat$attr$active == 1 & dat$attr$protection == 1)
  dat$epi$vac.num[at] <- sum(dat$attr$vaccination == 1)
  dat$epi$prot.num[at] <- sum(dat$attr$protection == 1)

##Troubleshooting Only - Delete prior to submitting pull request
  # check_for_false <- (sum(active == 1 & status == "s") + sum(active == 1 & status == "e") + sum(active == 1 & status == "i") + sum(active == 1 & status == "r") == network.size(dat$nw))
  # if(isFALSE(check_for_false)){
  #   s.num <- sum(active == 1 & status == "s")
  #   s.num
  #   e.num <- sum(active == 1 & status == "e")
  #   e.num
  #   i.num <- sum(active == 1 & status == "i")
  #   i.num
  #   r.num <- sum(active == 1 & status == "r")
  #   r.num
  #   network.size(dat$nw)
  #   status
  # }
  # dat$attr$entrTime
  # a<-1

  return(dat)
}
