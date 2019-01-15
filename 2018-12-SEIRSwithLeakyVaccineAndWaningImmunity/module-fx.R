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
  browser()

  ## Attributes ##
  active <- dat$attr$active
  status <- dat$attr$status
  protection.trans.prob <- dat$attr$protection.trans.prob

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

      #Create transmission probability dataframe
      transProb <- ifelse(!is.na(protection.trans.prob),protection.trans.prob,inf.prob)
      transProbDF <- data.frame(c(1:length(transProb)),transProb)
      colnames(transProbDF) <- c("sus","transProb")


      # Set parameters on discordant edgelist data frame
      del <- merge(x = del, y = transProbDF, by = "sus", all.x = TRUE)
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
        dat$attr$disease.experience[idsNewInf] <- "experienced"

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
  browser()

  ## Attributes ##
  active <- dat$attr$active
  status <- dat$attr$status
  infTime <- dat$attr$infTime

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
      infTime[idsSus] <- NA
    }
  }

  ## Write out updated status and infTime attributes ##
  dat$attr$status <- status

  ## Save summary statistics ##
  dat$epi$ei.flow[at] <- nInf
  dat$epi$ir.flow[at] <- nRec
  dat$epi$rs.flow[at] <- nSus

  return(dat)
}

# Update Death Module -----------------------------------------------------

dfunc <- function(dat, at) {

  browser()

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
  browser()

  ## Attributes ##
  active <- dat$attr$active
  status <- dat$attr$status
  infTime <- dat$attr$infTime
  entrTime <- dat$attr$entrTime
  exitTime <- dat$attr$exitTime
  vaccination <- dat$attr$vaccination
  protection <- dat$attr$protection
  protection.trans.prob <- dat$attr$protection.trans.prob
  time.of.vaccine.protection <- dat$attr$time.of.vaccine.protection
  disease.experience <- dat$attr$disease.experience

  ## Parameters ##
  n <- network.size(dat$nw)
  b.rate <- dat$param$birth.rate
  vaccination.rate.births <- dat$param$vaccination.rate.births
  protection.rate.births <- dat$param$protection.rate.births
  vaccination.rate.initialization <- dat$param$vaccination.rate.initialization
  protection.rate.initialization <- dat$param$protection.rate.initialization
  vaccination.rate.progression.disease.experienced <- dat$param$vaccination.rate.progression.disease.experienced
  vaccination.rate.progression.disease.naive <- dat$param$vaccination.rate.progression.disease.naive
  protection.rate.progression <- dat$param$protection.rate.progression
  leaky.degree.of.protection.max <- dat$param$leaky.degree.of.protection.max
  leaky.degree.of.protection.min <- dat$param$leaky.degree.of.protection.min
  time.to.vaccine.protection.decay <- dat$param$time.to.vaccine.protection.decay
  decay.error.tolerance <- dat$param$decay.error.tolerance

  ## Initializing Vaccination and Protection Process Flow Count Variables ##
  nVax.init <- 0
  nPrt.init <- 0
  nVax.prog <- 0
  nPrt.prog <- 0
  nVax.birth <- 0
  nPrt.birth <- 0
  lambda = -log(decay.error.tolerance)/time.to.vaccine.protection.decay

  ## INITIALIZATION OF VACCINATION AND PROTECTION VERTEX (NODE) ATTRIBUTES ##
  if (at == 2) {

    #Initialize vaccination and protection vectors
    vaccination <- rep(NA, n)
    protection <- rep(NA, n)
    protection.trans.prob <- rep(NA, n)
    time.of.vaccine.protection <- rep(NA, n)
    disease.experience <- rep(NA, n)

    # Determine individuals at time t=2 who are initially vaccinated
    idsEligVacInit <- which(active == 1)
    vecVacInit <- rbinom(length(idsEligVacInit), 1, vaccination.rate.initialization)
    idsVacInit <- idsEligVacInit[which(vecVacInit == 1)]
    vaccination[idsVacInit] <- "initial"

    #Determines if individual is protected based on
    #protection rate and infectious status
    idsEligProtInit <- which(vaccination == "initial" & status == 's')
    vecProtInit <- rbinom(length(idsEligProtInit), 1, protection.rate.initialization)
    idsProtInit <- idsEligProtInit[which(vecProtInit == 1)]
    idsNoProtInit <- setdiff(idsVacInit, idsProtInit)
    protection[idsProtInit] <- "initial"
    time.of.vaccine.protection[idsProtInit] <- runif(idsProtInit,-time.to.vaccine.protection.decay,at)
    protection.trans.prob[idsProtInit] <- (1 - exp(-lambda*(at - time.of.vaccine.protection[idsProtInit])))*((1 - leaky.degree.of.protection.min) - (1 - leaky.degree.of.protection.max)) + (1 - leaky.degree.of.protection.max)
    protection[idsNoProtInit] <- "none"

    #Determine individuals at time t=2 who are initially disease-experienced
    disease.experience[which(status != "s")] <- "experienced"

    #Captures the number of vaccinated and the number of protected (active)
    #individuals at the time of initialization
    nVax.init <- length(idsVacInit)
    nPrt.init <- length(idsProtInit)

  }


  ## VACCINATION PROGRESSION PROCESSES ##

  #Update the vaccination vector through the vaccination progression process
  idsEligVacProgDiseaseNaive <- which(is.na(vaccination) & active == 1 & is.na(disease.experience))
  idsEligVacProgDiseaseExperienced <- which(is.na(vaccination) & active == 1 & disease.experience == "experienced")
  vecVacProgDiseaseNaive <- rbinom(length(idsEligVacProgDiseaseNaive), 1, vaccination.rate.progression.disease.naive)
  vecVacProgDiseaseExperienced <- rbinom(length(idsEligVacProgDiseaseExperienced), 1, vaccination.rate.progression.disease.experienced)
  idsVacProg <- c(idsEligVacProgDiseaseNaive[which(vecVacProgDiseaseNaive == 1)], idsEligVacProgDiseaseExperienced[which(vecVacProgDiseaseExperienced == 1)])
  vaccination[idsVacProg] <- "progress"

  #Update the protection vector through the vaccination protection progression
  #process
  idsEligProtProg <- which(vaccination == "progress" &  is.na(protection)
                           & status == 's')
  vecProtProg <- rbinom(length(idsEligProtProg), 1, protection.rate.progression)
  idsProtProg <- idsEligProtProg[which(vecProtProg == 1)]
  idsNoProtProg <- setdiff(idsVacProg, idsProtProg)
  protection[idsProtProg] <- "progress"
  time.of.vaccine.protection[idsProtProg] <- at
  protection.trans.prob[idsProtProg] <- (1 - exp(-lambda*(at - time.of.vaccine.protection[idsProtProg])))*((1 - leaky.degree.of.protection.min) - (1 - leaky.degree.of.protection.max)) + (1 - leaky.degree.of.protection.max)
  protection[idsNoProtProg] <- "none"

  #Captures the total number of vaccinated and protected individuals
  #after running vaccination and protection processes for current time step
  nVax.prog <- length(idsVacProg)
  nPrt.prog <- length(idsProtProg)


  ## BIRTH AND BIRTH VACCINATION PROCESSES ##

  #Birth Process
  nBirthsExp <- n * b.rate
  nBirths <- rpois(1, nBirthsExp)
  nVax.birth <- 0
  nPrt.birth <- 0

  if (nBirths > 0) {
    dat$nw <- add.vertices(dat$nw, nv = nBirths)
    newNodes <- (n + 1):(n + nBirths)
    dat$nw <- activate.vertices(dat$nw, onset = at, terminus = Inf, v = newNodes)
  }

  #Update attributes
  if (nBirths > 0) {
    active <- c(active, rep(1, nBirths))
    status <- c(status, rep("s", nBirths))
    infTime <- c(infTime, rep(NA, nBirths))
    entrTime <- c(entrTime, rep(at, nBirths))
    exitTime <- c(exitTime, rep(NA, nBirths))
    vaccination <- c(vaccination, rep(NA, nBirths))
    protection <- c(protection, rep(NA, nBirths))

    # New birth vaccination process and count
    vaccinatedNewBirths <- rbinom(nBirths, 1, vaccination.rate.births)
    vaccination[newNodes] <- ifelse(vaccinatedNewBirths == 1, "birth", NA)
    nVax.birth <- length(which(vaccinatedNewBirths == 1))

    #Update the protection vector through the vaccination protection at birth
    #process
    idsEligProtBirth <- which(vaccination == "birth" & status == "s" & is.na(protection))
    vecProtBirth <- rbinom(length(idsEligProtBirth), 1, protection.rate.births)
    idsProtBirth <- idsEligProtBirth[which(vecProtBirth == 1)]
    idsNoProtBirth <- idsEligProtBirth[which(vecProtBirth == 0)]
    protection[idsProtBirth] <- "birth"
    time.of.vaccine.protection[idsProtBirth] <- at
    protection.trans.prob[idsProtBirth] <- (1 - exp(-lambda*(at - time.of.vaccine.protection[idsProtBirth])))*((1 - leaky.degree.of.protection.min) - (1 - leaky.degree.of.protection.max)) + (1 - leaky.degree.of.protection.max)
    protection[idsNoProtBirth] <- "none"
    nPrt.birth <- length(which(vecProtBirth == 1))

  }

  ## UPDATE NODE ATTRIBUTES ##

  dat$attr$status <- ifelse(status == "s"
                            & protection %in% c("initial", "progress", "birth")
                            & active == 1, "v", status)
  dat$attr$vaccination <- vaccination
  dat$attr$protection <- protection
  dat$attr$time.of.vaccine.protection <- time.of.vaccine.protection
  dat$attr$active <- active
  dat$attr$infTime <- infTime
  dat$attr$entrTime <- entrTime
  dat$attr$exitTime <- exitTime
  dat$attr$protection.trans.prob <- protection.trans.prob
  dat$attr$disease.experience <- disease.experience

  ## SUMMARY STATISTICS ##

  #Births
  dat$epi$b.flow[at] <- nBirths
  if (at == 2) {
    dat$epi$b.num[at] <- nBirths
  } else {
    dat$epi$b.num[at] <- dat$epi$b.num[at - 1] + nBirths
  }

  #Vaccination and Protection
  dat$epi$vac.flow[at] <- nVax.init + nVax.prog + nVax.birth
  dat$epi$prt.flow[at] <- nPrt.init + nPrt.prog + nPrt.birth
  dat$epi$vac.num[at] <- sum(active == 1 & vaccination %in% c("initial", "progress", "birth"))
  dat$epi$prt.num[at] <- sum(active == 1 & protection %in% c("initial", "progress", "birth"))

  dat$epi$vac.init.flow[at] <- nVax.init
  dat$epi$prt.init.flow[at] <- nPrt.init
  dat$epi$vac.prog.flow[at] <- nVax.prog
  dat$epi$prt.prog.flow[at] <- nPrt.prog
  dat$epi$vac.birth.flow[at] <- nVax.birth
  dat$epi$prt.birth.flow[at] <- nPrt.birth

  dat$epi$vac.init.num[at] <- sum(active == 1 & !is.na(vaccination) & vaccination == "initial")
  dat$epi$prt.init.num[at] <- sum(active == 1 & !is.na(protection) & protection == "initial")
  dat$epi$vac.prog.num[at] <- sum(active == 1 & !is.na(vaccination) & vaccination == "progress")
  dat$epi$prt.prog.num[at] <- sum(active == 1 & !is.na(protection) & protection == "progress")
  dat$epi$vac.birth.num[at] <- sum(active == 1 & !is.na(vaccination) & vaccination == "birth")
  dat$epi$prt.birth.num[at] <- sum(active == 1 & !is.na(protection) & protection == "birth")

  return(dat)
}
