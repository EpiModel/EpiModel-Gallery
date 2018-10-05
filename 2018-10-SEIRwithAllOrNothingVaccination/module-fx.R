##
## SEIR Model with Vital Dynamics and an All or Nothing Vaccine Implementation
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
  idsEligInf <- which(active == 1 & status == "e" & dat$attr$protected == "vulnerable")
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

  # browser()

#Toggle for step-through debugging

  ## Attributes ##
  active <- dat$attr$active
  status <- dat$attr$status

  ## Parameters ##
  n <- network.size(dat$nw)
  b.rate <- dat$param$birth.rate
  vaccine.rate <- dat$param$vaccine.rate
  vaccine.efficacy <- dat$param$vaccine.efficacy

  ## Initialization of Vaccinated and Protected Vertex (Node) Attributes
  if (at == 2) {

    # Pull vaccination and protection attributes from the fitted network model
    vaccination.process <- rbinom(n, 1, vaccine.rate)

    #Vector in which status of "i"=0 and not"i" = 1
    numericInfectionStatusVector <- rep(1,n)
    infInitVec <- which(status == "i")
    if(length(infInitVec) > 1){
      numericInfectionStatusVector[infInitVec] <- 0
    }

    #Determines if individual is protected based on vaccine efficacy and status
    protection.process <- vaccination.process * rbinom(n, 1, vaccine.efficacy) * numericInfectionStatusVector

    dat$attr$vaccinated <- c(replace(replace(vaccination.process, vaccination.process==0,"unvaccinated"),vaccination.process==1,"vaccinated"))
    dat$attr$protected <- c(replace(replace(protection.process, protection.process==0,"vulnerable"),protection.process==1,"protected"))

    nVaccinated <- length(which(vaccination.process == 1))
    nProtected <- length(which(protection.process == 1))

  }

  ## Process ##
  nBirthsExp <- n * b.rate
  nBirths <- rpois(1, nBirthsExp)
  nVaccinated <- 0
  nProtected <- 0

  if (nBirths > 0) {
    dat$nw <- add.vertices(dat$nw, nv = nBirths)
    newNodes <- (n + 1):(n + nBirths)
    dat$nw <- activate.vertices(dat$nw, onset = at, terminus = Inf, v = newNodes)
  }

  # Update attributes
  if (nBirths > 0) {
    dat$attr$active <- c(active, rep(1, nBirths))
    dat$attr$status <- c(status, rep("s", nBirths))
    dat$attr$infTime <- c(dat$attr$infTime, rep(NA, nBirths))
    dat$attr$entrTime <- c(dat$attr$entrTime, rep(at, nBirths))
    dat$attr$exitTime <- c(dat$attr$exitTime, rep(NA, nBirths))

    vaccination.process <- rbinom(nBirths, 1, vaccine.rate)
    dat$attr$vaccinated <- c(dat$attr$vaccinated, replace(replace(vaccination.process, vaccination.process==0,"unvaccinated"),vaccination.process==1,"vaccinated"))

    protection.process <- vaccination.process * rbinom(nBirths, 1, vaccine.efficacy)
    dat$attr$protected <- c(dat$attr$protected, replace(replace(protection.process, protection.process==0,"vulnerable"),protection.process==1,"protected"))

    nVaccinated <- length(which(vaccination.process == 1))
    nProtected <- length(which(protection.process == 1))

  }

  ## Summary statistics ##
  dat$epi$b.flow[at] <- nBirths
  dat$epi$vaccinated.flow[at] <- nVaccinated
  dat$epi$protected.flow[at] <- nProtected
  dat$epi$vaccinated.num[at] <- sum(dat$attr$active == 1 & dat$attr$vaccinated == "vaccinated")
  dat$epi$protected.num[at] <- sum(dat$attr$active == 1 & dat$attr$protected == "protected")


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
