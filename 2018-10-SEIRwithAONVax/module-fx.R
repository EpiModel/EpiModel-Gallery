##
## SEIR Model with Vital Dynamics and an All or Nothing Vaccine Implementation
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##
## Authors: Connor M. Van Meter
## Date: November 2018
##


# Attribute initialization function ---------------------------------------

init_new_attrs <- function(dat, at) {

  ## INITIALIZATION OF VACCINATION AND PROTECTION VERTEX (NODE) ATTRIBUTES ##
  if (is.null(get_attr(dat, "vaccination", override.null.error = TRUE))) {

    # Pull attr
    active <- get_attr(dat, "active")
    n <- sum(active == 1)
    status <- get_attr(dat, "status")

    # Pull params
    vaccination.rate.init <- get_param(dat, "vaccination.rate.initialization")
    protection.rate.init <- get_param(dat, "protection.rate.initialization")

    #Initialize vaccination and protection vectors
    vaccination <- rep(NA, n)
    protection <- rep(NA, n)

    # Determine individuals at time t=2 who are initially vaccinated - Sam's method
    idsEligVacInit <- which(active == 1)
    vecVacInit <- rbinom(length(idsEligVacInit), 1, vaccination.rate.init)
    idsVacInit <- idsEligVacInit[which(vecVacInit == 1)]
    vaccination[idsVacInit] <- "initial"

    #Determines if individual is protected based on
    #protection rate and infectious status
    idsEligProtInit <- which(vaccination == "initial" & status == "s")
    vecProtInit <- rbinom(length(idsEligProtInit), 1, protection.rate.init)
    idsProtInit <- idsEligProtInit[which(vecProtInit == 1)]
    idsNoProtInit <- setdiff(idsVacInit, idsProtInit)
    protection[idsProtInit] <- "initial"
    protection[idsNoProtInit] <- "none"

    #Output vaccination/protection attributes
    dat <- set_attr(dat, "vaccination", vaccination)
    dat <- set_attr(dat, "protection", protection)

  }

  return(dat)
}



# Replacement infection/transmission module -------------------------------

infect <- function(dat, at) {

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

  ## Initialize default incidence at 0 ##
  nInf <- 0

  ## If any infected nodes, proceed with transmission ##
  if (nElig > 0 && nElig < nActive) {

    ## Look up discordant edgelist ##
    del <- discord_edgelist(dat, at)

    ## If any discordant pairs, proceed ##
    if (!(is.null(del)) && nrow(del) > 0) {

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
        status[idsNewInf] <- "e"
        dat <- set_attr(dat, "status", status)
        infTime[idsNewInf] <- at
        dat <- set_attr(dat, "infTime", infTime)
      }
    }
  }

  ## Save summary statistics
  dat <- set_epi(dat, "se.flow", at, nInf)

  dat <- set_epi(dat, "s.num", at, sum(active == 1 & status == "s"))
  dat <- set_epi(dat, "e.num", at, sum(active == 1 & status == "e"))
  dat <- set_epi(dat, "i.num", at, sum(active == 1 & status == "i"))
  dat <- set_epi(dat, "r.num", at, sum(active == 1 & status == "r"))
  dat <- set_epi(dat, "v.num", at, sum(active == 1 & status == "v"))

  return(dat)

}


# New disease progression module ------------------------------------------
# (Replaces the recovery module)

progress <- function(dat, at) {

  ## Attributes ##
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")

  ## Parameters ##
  ei.rate <- get_param(dat, "ei.rate")
  ir.rate <- get_param(dat, "ir.rate")

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
  dat <- set_attr(dat, "status", status)

  ## Save summary statistics ##
  dat <- set_epi(dat, "ei.flow", at, nInf)
  dat <- set_epi(dat, "ir.flow", at, nRec)

  return(dat)
}

# Update Death Module -----------------------------------------------------

dfunc <- function(dat, at) {

  ## Attributes ##
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  exitTime <- get_attr(dat, "exitTime")

  ## Parameters ##
  mortality.rate <- get_param(dat, "mortality.rate")
  mort.rates <- rep(mortality.rate, length(active))
  mort.dis.mult <- get_param(dat, "mortality.disease.mult")

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

    ## Update nodal attributes
    if (nDeaths > 0) {
      active[idsDeaths] <- 0
      exitTime[idsDeaths] <- at
    }

  }

  ## Write out updated attributes
  dat <- set_attr(dat, "active", active)
  dat <- set_attr(dat, "exitTime", exitTime)

  ## Summary statistics
  dat <- set_epi(dat, "d.flow", at, nDeaths)
  dat <- set_epi(dat, "d.num", at, sum(active == 0))

  return(dat)
}


# Updated Birth Module ----------------------------------------------------

bfunc <- function(dat, at) {

  ## Attributes ##
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  infTime <- get_attr(dat, "infTime")

  ## Parameters ##
  n <- length(active)
  b.rate <- get_param(dat, "birth.rate")
  vaccination.rate.births <- get_param(dat, "vaccination.rate.births")
  protection.rate.births <- get_param(dat, "protection.rate.births")
  vaccination.rate.progression <- get_param(dat, "vaccination.rate.progression")
  protection.rate.progression <- get_param(dat, "protection.rate.progression")

  ## Initializing Vaccination and Protection Process Flow Count Variables ##
  nVax.prog <- 0
  nPrt.prog <- 0
  nVax.birth <- 0
  nPrt.birth <- 0

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
  idsEligProtProg <- which(vaccination == "progress" &  is.na(protection) & status == "s")
  vecProtProg <- rbinom(length(idsEligProtProg), 1, protection.rate.progression)
  idsProtProg <- idsEligProtProg[which(vecProtProg == 1)]
  idsNoProtProg <- setdiff(idsVacProg, idsProtProg)
  protection[idsProtProg] <- "progress"
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

  #Update attributes
  if (nBirths > 0) {
    dat <- append_core_attr(dat, at, nBirths)

    newNodes <- (n + 1):(n + nBirths)
    status <- c(status, rep("s", nBirths))
    infTime <- c(infTime, rep(NA, nBirths))
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
    protection[idsNoProtBirth] <- "none"
    nPrt.birth <- length(which(vecProtBirth == 1))

  }

  ## UPDATE NODE ATTRIBUTES ##
  active <- get_attr(dat, "active")
  statusV <- which(status == "s" & protection %in% c("initial", "progress", "birth") &
                   active == 1)
  status[statusV] <- "v"

  dat <- set_attr(dat, "status", status)
  dat <- set_attr(dat, "vaccination", vaccination)
  dat <- set_attr(dat, "protection", protection)
  dat <- set_attr(dat, "infTime", infTime)


  ## SUMMARY STATISTICS ##

  # Births
  dat <- set_epi(dat, "a.flow", at, nBirths)
  cumm.births <- get_epi(dat, "a.flow")
  dat <- set_epi(dat, "b.num", at, sum(cumm.births, na.rm = TRUE))

  # Vaccination and Protection
  dat <- set_epi(dat, "vac.flow", at, nVax.prog + nVax.birth)
  dat <- set_epi(dat, "prt.flow", at, nPrt.prog + nPrt.birth)
  dat <- set_epi(dat, "vac.num", at,
                 sum(active == 1 & vaccination %in% c("initial", "progress", "birth")))
  dat <- set_epi(dat, "prt.num", at,
                 sum(active == 1 & protection %in% c("initial", "progress", "birth")))

  dat <- set_epi(dat, "vac.prog.flow", at, nVax.prog)
  dat <- set_epi(dat, "prt.prog.flow", at, nPrt.prog)
  dat <- set_epi(dat, "vac.birth.flow", at, nVax.birth)
  dat <- set_epi(dat, "prt.birth.flow", at, nPrt.birth)

  dat <- set_epi(dat, "vac.init.num", at,
                 sum(active == 1 & !is.na(vaccination) & vaccination == "initial"))
  dat <- set_epi(dat, "prt.init.num", at,
                 sum(active == 1 & !is.na(protection) & protection == "initial"))
  dat <- set_epi(dat, "vac.prog.num", at,
                 sum(active == 1 & !is.na(vaccination) & vaccination == "progress"))
  dat <- set_epi(dat, "prt.prog.num", at,
                 sum(active == 1 & !is.na(protection) & protection == "progress"))
  dat <- set_epi(dat, "vac.birth.num", at,
                 sum(active == 1 & !is.na(vaccination) & vaccination == "birth"))
  dat <- set_epi(dat, "prt.birth.num", at,
                 sum(active == 1 & !is.na(protection) & protection == "birth"))

  return(dat)
}
