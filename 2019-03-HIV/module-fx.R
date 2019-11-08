##
## HIV Transmission Model
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##
## Authors: Connor Van Meter, Samuel M. Jenness, Yuan Zhao
## Date: March 2019
##


# Replacement infection/transmission module -------------------------------

infect <- function(dat, at) {

  ## Uncomment this to function environment interactively
  #browser()

  
  ## Attributes ##
  active <- dat$attr$active
  status <- dat$attr$status
  ART.status <- dat$attr$ART.status
  stage <- dat$attr$stage
  stage.time <- dat$attr$stage.time
  ART.time <- dat$attr$ART.time


  ## Parameters ##
  inf.prob.chronic <- dat$param$inf.prob.chronic
  relative.inf.prob.acute <- dat$param$relative.inf.prob.acute
  relative.inf.prob.AIDS <- dat$param$relative.inf.prob.AIDS
  relative.inf.prob.ART <- dat$param$relative.inf.prob.ART
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
      # dependent on HIV status and ART status

        #Create Empty transProb Column
        del$transProb <- rep(NA,length(del$inf))

        #Trans Prob - Acute, No ART
        idsAcuteNoART <- which(stage[del$inf] == "acute" &
                                 ART.status[del$inf] == 0)
        del$transProb[idsAcuteNoART] <- inf.prob.chronic *
          relative.inf.prob.acute

        #Trans Prob - Chronic, No ART
        idsChronicNoART <- which(stage[del$inf] %in%
                                   c("chronic1", "chronic2") &
                                   ART.status[del$inf] == 0)
        del$transProb[idsChronicNoART] <- inf.prob.chronic

        #Trans Prob - AIDS, No ART
        idsAIDSNoART <- which(stage[del$inf] == "AIDS" &
                                ART.status[del$inf] == 0)
        del$transProb[idsAIDSNoART] <- inf.prob.chronic * relative.inf.prob.AIDS

        #Trans Prob - Acute, ART
        idsAcuteART <- which(stage[del$inf] == "acute" &
                               ART.status[del$inf] == 1)
        del$transProb[idsAcuteART] <- inf.prob.chronic *
                                      relative.inf.prob.acute *
                                      relative.inf.prob.ART

        #Trans Prob - Chronic, ART
        idsChronicART <- which(stage[del$inf] %in% c("chronic1","chronic2") &
                                 ART.status[del$inf] == 1)
        del$transProb[idsChronicART] <- inf.prob.chronic * relative.inf.prob.ART

        #Trans Prob - AIDS, ART
        idsAIDSART <- which(stage[del$inf] == "AIDS" & ART.status[del$inf] == 1)
        del$transProb[idsAIDSART] <- inf.prob.chronic *
                                     relative.inf.prob.AIDS *
                                     relative.inf.prob.ART

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
        stage[idsNewInf] <- "acute"
        status[idsNewInf] <- "i"
        stage.time[idsNewInf] <- 0
        ART.status[idsNewInf] <- 0
        ART.time[idsNewInf] <- 0

      }
    }
  }

  
  #Increment ART time by 1
  ART.time = ifelse(!is.na(ART.time), ART.time + 1, ART.time)

  #Increment status time by 1
  stage.time = ifelse(!is.na(stage.time), stage.time + 1, stage.time)

  
  ## Update attributes
  dat$attr$stage <- stage
  dat$attr$status <- status
  dat$attr$stage.time <- stage.time
  dat$attr$ART.status <- ART.status
  dat$attr$ART.time <- ART.time

  
  ## Save summary statistics
  dat$epi$acute.flow[at] <- nInf
  dat$epi$s.num[at] <- sum(dat$attr$active == 1 &
                             dat$attr$status == "s")
  dat$epi$acute.ART.num[at] <- sum(dat$attr$active == 1 &
                                     dat$attr$stage == "acute" &
                                     dat$attr$ART.status == 1, na.rm = TRUE)
  dat$epi$acute.NoART.num[at] <- sum(dat$attr$active == 1 &
                                       dat$attr$stage == "acute" &
                                       dat$attr$ART.status == 0, na.rm = TRUE)
  dat$epi$chronic1.ART.num[at] <- sum(dat$attr$active == 1 &
                                        dat$attr$stage == "chronic1" &
                                        dat$attr$ART.status == 1, na.rm = TRUE)
  dat$epi$chronic1.NoART.num[at] <- sum(dat$attr$active == 1 &
                                          dat$attr$stage == "chronic1" &
                                          dat$attr$ART.status == 0,
                                          na.rm = TRUE)
  dat$epi$chronic2.ART.num[at] <- sum(dat$attr$active == 1 &
                                        dat$attr$stage == "chronic2" &
                                        dat$attr$ART.status == 1,
                                        na.rm = TRUE)
  dat$epi$chronic2.NoART.num[at] <- sum(dat$attr$active == 1 &
                                          dat$attr$stage == "chronic2" &
                                          dat$attr$ART.status == 0,
                                          na.rm = TRUE)
  dat$epi$AIDS.ART.num[at] <- sum(dat$attr$active == 1 &
                                     dat$attr$stage == "AIDS" &
                                     dat$attr$ART.status == 1, na.rm = TRUE)
  dat$epi$AIDS.NoART.num[at] <- sum(dat$attr$active == 1 &
                                       dat$attr$stage == "AIDS" &
                                       dat$attr$ART.status == 0, na.rm = TRUE)

  return(dat)

}


# HIV progression module ------------------------------------------
# (Replaces the recovery module)

progress <- function(dat, at) {

  ## Uncomment this to function environment interactively
  #browser()
  

  ## Attributes ##
  active <- dat$attr$active
  status <- dat$attr$status


  ## Initiating indicators to track stage of HIV and ART status ##
  if (at == 2) {

    dat$attr$stage <- ifelse(dat$attr$status == "i","acute",NA)
    dat$attr$ART.status <- ifelse(dat$attr$status == "i",0,NA)
    dat$attr$ART.time <- ifelse(!is.na(dat$attr$ART.status) &
                                                 dat$attr$ART.status == 0,
                                               1, NA)
    dat$attr$stage.time <- ifelse(dat$attr$status == "i",1,NA)
  }


  ## Attributes ##
  ART.status <- dat$attr$ART.status
  stage <- dat$attr$stage
  ART.time <- dat$attr$ART.time
  stage.time <- dat$attr$stage.time


  ## Progression Parameters ##
  AcuteToChronic1.Rate <- dat$param$AcuteToChronic1.Rate
  Chronic1ToChronic2.Rate <- dat$param$Chronic1ToChronic2.Rate
  Chronic2ToAIDS.Rate <- dat$param$Chronic2ToAIDS.Rate
  ART.Treatment.Rate <- dat$param$ART.Treatment.Rate
  ART.Discontinuance.Rate <- dat$param$ART.Discontinuance.Rate
  ART.Progression.Reduction.Rate <- dat$param$ART.Progression.Reduction.Rate

  
  ## Initialize vectors ##
  vecART <- vector()
  vecARTDisc <- vector()
  idsART <- vector()
  idsARTDisc <- vector()
  vecChronic1ART <- vector()
  idsChronic1ART <- vector()
  vecChronic1NoART <- vector()
  idsChronic1NoART <- vector()
  vecChronic2ART <- vector()
  idsChronic2ART <- vector()
  vecChronic2NoART <- vector()
  idsChronic2NoART <- vector()
  vecAIDSART <- vector()
  idsAIDSART <- vector()
  vecAIDSNoART <- vector()
  idsAIDSNoART <- vector()

  
  ## ART Treatment ##
  idsEligART <- which(active == 1 & status == "i" &
                        ART.time != 0 &
                        ART.status == 0 & !is.na(ART.status))

  if (length(idsEligART) > 0) {
    vecART <- which(rbinom(length(idsEligART), 1, ART.Treatment.Rate) == 1)
    if (length(vecART) > 0) {
      idsART <- idsEligART[vecART]
      ART.status[idsART] <- 1
      ART.time[idsART] <- 0
    }
  }

  
  ## ART Discontinuance ##
  idsEligARTDisc <- which(active == 1 & status == "i" &
                            ART.time != 0 &
                            ART.status == 1 & !is.na(ART.status))

  if (length(idsEligARTDisc) > 0) {
    vecARTDisc <- which(rbinom(length(idsEligARTDisc), 1,
                               ART.Discontinuance.Rate) == 1)
    if (length(vecARTDisc) > 0) {
      idsARTDisc <- idsEligARTDisc[vecARTDisc]
      ART.status[idsARTDisc] <- 0
      ART.time[idsARTDisc] <- 0
    }
  }

  
  ##ART Treatment/Discontinuance Flows within HIV Status Subcompartment
  acute.ART.treatment.flow <- length(intersect(idsART, which(stage == "acute")))
  acute.ART.discontinuance.flow <- length(intersect(idsARTDisc, 
                                                    which(stage == "acute")))
  chronic1.ART.treatment.flow <- length(intersect(idsART, 
                                                  which(stage == "chronic1")))
  chronic1.ART.discontinuance.flow <- length(intersect(idsARTDisc, 
                                                    which(stage == "chronic1")))
  chronic2.ART.treatment.flow <- length(intersect(idsART, 
                                                  which(stage == "chronic2")))
  chronic2.ART.discontinuance.flow <- length(intersect(idsARTDisc, 
                                                    which(stage == "chronic2")))
  AIDS.ART.treatment.flow <- length(intersect(idsART, which(stage == "AIDS")))
  AIDS.ART.discontinuance.flow <- length(intersect(idsARTDisc, 
                                                   which(stage == "AIDS")))


  ## Acute to chronic 1 stage progression process ##
  idsEligChronic1ART <- which(active == 1 & stage == "acute" &
                                stage.time != 0 & ART.status == 1 &
                                !is.na(ART.status))
  idsEligChronic1NoART <- which(active == 1 & stage == "acute" &
                                  stage.time != 0 & ART.status == 0 &
                                  !is.na(ART.status))

  if (length(idsEligChronic1ART) > 0) {
    vecChronic1ART <- which(rbinom(length(idsEligChronic1ART), 1, 
                                   AcuteToChronic1.Rate *
                                     ART.Progression.Reduction.Rate) == 1)

    if (length(vecChronic1ART) > 0) {
      idsChronic1ART <- idsEligChronic1ART[vecChronic1ART]
      stage[idsChronic1ART] <- "chronic1"
      stage.time[idsChronic1ART] <- 0
    }
  }

  if (length(idsEligChronic1NoART) > 0) {
    vecChronic1NoART <- which(rbinom(length(idsEligChronic1NoART), 1,
                                     AcuteToChronic1.Rate) == 1)

    if (length(vecChronic1NoART) > 0) {
      idsChronic1NoART <- idsEligChronic1NoART[vecChronic1NoART]
      stage[idsChronic1NoART] <- "chronic1"
      stage.time[idsChronic1NoART] <- 0
    }
  }

  
  ## Chronic 1 to chronic 2 stage progression process ##
  idsEligChronic2ART <- which(active == 1 & stage == "chronic1" &
                                stage.time != 0 & ART.status == 1 &
                                !is.na(ART.status))
  idsEligChronic2NoART <- which(active == 1 & stage == "chronic1" &
                                  stage.time != 0 & ART.status == 0 &
                                  !is.na(ART.status))

  if (length(idsEligChronic2ART) > 0) {
    vecChronic2ART <- which(rbinom(length(idsEligChronic2ART), 1,
                                   Chronic1ToChronic2.Rate *
                                     ART.Progression.Reduction.Rate) == 1)

    if (length(vecChronic2ART) > 0) {
      idsChronic2ART <- idsEligChronic2ART[vecChronic2ART]
      stage[idsChronic2ART] <- "chronic2"
      stage.time[idsChronic2ART] <- 0
    }
  }

  if (length(idsEligChronic2NoART) > 0) {
    vecChronic2NoART <- which(rbinom(length(idsEligChronic2NoART), 1,
                                     Chronic1ToChronic2.Rate) == 1)

    if (length(vecChronic2NoART) > 0) {
      idsChronic2NoART <- idsEligChronic2NoART[vecChronic2NoART]
      stage[idsChronic2NoART] <- "chronic2"
      stage.time[idsChronic2NoART] <- 0
    }
  }

  
  ## Chronic 2 to AIDS stage progression process ##
  idsEligAIDSART <- which(active == 1 & stage == "chronic2" &
                            stage.time != 0 & ART.status == 1 &
                             !is.na(ART.status))
  idsEligAIDSNoART <- which(active == 1 & stage == "chronic2" &
                              stage.time != 0 & ART.status == 0 &
                               !is.na(ART.status))

  if (length(idsEligAIDSART) > 0) {
    vecAIDSART <- which(rbinom(length(idsEligAIDSART), 1, Chronic2ToAIDS.Rate *
                                  ART.Progression.Reduction.Rate) == 1)

    if (length(vecAIDSART) > 0) {
      idsAIDSART <- idsEligAIDSART[vecAIDSART]
      stage[idsAIDSART] <- "AIDS"
      stage.time[idsAIDSART] <- 0
    }
  }

  if (length(idsEligAIDSNoART) > 0) {
    vecAIDSNoART <- which(rbinom(length(idsEligAIDSNoART), 1,
                                  Chronic2ToAIDS.Rate) == 1)

    if (length(vecAIDSNoART) > 0) {
      idsAIDSNoART <- idsEligAIDSNoART[vecAIDSNoART]
      stage[idsAIDSNoART] <- "AIDS"
      stage.time[idsAIDSNoART] <- 0
    }
  }


  #Save attributes
  dat$attr$stage <- stage
  dat$attr$stage.time <- stage.time
  dat$attr$ART.status <- ART.status
  dat$attr$ART.time <- ART.time

  #Update summary statistics
  dat$epi$acute.ART.treatment.flow[at] <- acute.ART.treatment.flow
  dat$epi$acute.ART.discontinuance.flow[at] <- acute.ART.discontinuance.flow
  dat$epi$chronic1.ART.treatment.flow[at] <- chronic1.ART.treatment.flow
  dat$epi$chronic1.ART.discontinuance.flow[at] <-
    chronic1.ART.discontinuance.flow
  dat$epi$chronic2.ART.treatment.flow[at] <- chronic2.ART.treatment.flow
  dat$epi$chronic2.ART.discontinuance.flow[at] <-
    chronic2.ART.discontinuance.flow
  dat$epi$AIDS.ART.treatment.flow[at] <- AIDS.ART.treatment.flow
  dat$epi$AIDS.ART.discontinuance.flow[at] <- AIDS.ART.discontinuance.flow
  dat$epi$chronic1.ART.flow[at] <- ifelse(length(idsEligChronic1ART) > 0 &
                                          length(vecChronic1ART) > 0,
                                          length(idsChronic1ART), 0)
  dat$epi$chronic1.NoART.flow[at] <- ifelse(length(idsEligChronic1NoART) > 0 &
                                            length(vecChronic1NoART) > 0,
                                            length(idsChronic1NoART), 0)
  dat$epi$chronic2.ART.flow[at] <- ifelse(length(idsEligChronic2ART) > 0 &
                                          length(vecChronic2ART) > 0,
                                          length(idsChronic2ART), 0)
  dat$epi$chronic2.NoART.flow[at] <- ifelse(length(idsEligChronic2NoART) > 0 &
                                            length(vecChronic2NoART) > 0,
                                            length(idsChronic2NoART), 0)
  dat$epi$AIDS.ART.flow[at] <- ifelse(length(idsEligAIDSART) > 0 &
                                      length(vecAIDSART) > 0,
                                      length(idsAIDSART), 0)
  dat$epi$AIDS.NoART.flow[at] <- ifelse(length(idsEligAIDSNoART) > 0 &
                                        length(vecAIDSNoART) > 0,
                                        length(idsAIDSNoART), 0)

  return(dat)
}


dfunc <- function(dat, at) {

  #browser()

  
  ## Attributes ##
  active <- dat$attr$active
  status <- dat$attr$status
  ART.status <- dat$attr$ART.status
  stage <- dat$attr$stage
  stage.time <- dat$attr$stage.time

  
  ## Parameters ##
  Departure.rates <- rep(dat$param$departure.rate, network.size(dat$nw))
  ART.Progression.Reduction.Rate <- dat$param$ART.Progression.Reduction.Rate
  AIDSToDepart.Rate <- dat$param$AIDSToDepart.Rate

  
  ## Initialize vectors ##
  vecDeparture <- vector()
  idsDeparture <- vector()
  vecDepartAIDSART <- vector()
  idsDepartAIDSART <- vector()
  vecDepartAIDSNoART <- vector()
  idsDepartAIDSNoART <- vector()


  ## Query active individuals not in the AIDS stage of HIV ##
  idsEligDepartStandard <- which(active == 1 & ((stage != "AIDS" &
                                                  !is.na(stage)) |
                                                  status == "s"))

  
  ## Departure process for individuals not in the AIDS stage of HIV
  if (length(idsEligDepartStandard) > 0) {
    Departure_rates_of_elig <- Departure.rates[idsEligDepartStandard]
    vecDeparture <- which(rbinom(idsEligDepartStandard, 1,
                                 Departure_rates_of_elig) == 1)

    if (length(vecDeparture) > 0) {
      idsDeparture <- idsEligDepartStandard[vecDeparture]
    }
  }

  
  ## Departure process for individuals in the AIDS stage of HIV on ART
  idsEligDepartAIDSART <- which(active == 1 & stage == "AIDS" &
                                  stage.time != 0 & ART.status == 1 &
                                   !is.na(ART.status))

  if (length(idsEligDepartAIDSART) > 0) {
    vecDepartAIDSART <- which(rbinom(length(idsEligDepartAIDSART), 1,
                                      AIDSToDepart.Rate *
                                        ART.Progression.Reduction.Rate) == 1)

    if (length(vecDepartAIDSART) > 0) {
      idsDepartAIDSART <- idsEligDepartAIDSART[vecDepartAIDSART]
    }
  }

  
  ## Departure process for individuals in the AIDS stage of HIV not on ART
  idsEligDepartAIDSNoART <- which(active == 1 & stage == "AIDS" &
                                    stage.time != 0 & ART.status == 0 &
                                    !is.na(ART.status))

  if (length(idsEligDepartAIDSNoART) > 0) {
    vecDepartAIDSNoART <- which(rbinom(length(idsEligDepartAIDSNoART), 1,
                                        AIDSToDepart.Rate) == 1)

    if (length(vecDepartAIDSNoART) > 0) {
      idsDepartAIDSNoART <- idsEligDepartAIDSNoART[vecDepartAIDSNoART]
  }
}

  
##Save departure summary statistics
  dat$epi$depart.standard.ART.flow[at] <- 
    ifelse(length(idsEligDepartStandard) > 0 & length(vecDeparture) > 0,
           length(which(ART.status[idsDeparture] == 1)), 0)
  dat$epi$depart.standard.NoART.flow[at] <- 
    ifelse(length(idsEligDepartStandard) > 0 & length(vecDeparture) > 0, 
           length(which(ART.status[idsDeparture] == 0 |
                          is.na(ART.status[idsDeparture]))), 0)
  dat$epi$depart.AIDS.ART.flow[at] <- 
    ifelse(length(idsEligDepartAIDSART) > 0 & length(vecDepartAIDSART) > 0,
           length(idsDepartAIDSART), 0)
  dat$epi$depart.AIDS.NoART.flow[at] <- 
    ifelse(length(idsEligDepartAIDSNoART) > 0 & length(vecDepartAIDSART) > 0,
           length(idsDepartAIDSNoART), 0)

  
  ## Update nodal attributes on attr and networkDynamic object ##
  if (length(idsDeparture) > 0 | length(idsDepartAIDSART) > 0 |
      length(idsDepartAIDSNoART) > 0) {
    idsDeparted <- c(idsDeparture, idsDepartAIDSART, idsDepartAIDSNoART)
    active[idsDeparted] <- 0
    dat$nw <- deactivate.vertices(dat$nw, onset = at, terminus = Inf,
                                  v = idsDeparted, deactivate.edges = TRUE)
  }

  
  ## Save updated status attribute ##
  dat$attr$active <- active

  return(dat)
}


afunc <- function(dat, at) {

  #Toggle for step-through debugging
  #browser()

  
  ## Attributes ##
  active <- dat$attr$active
  status <- dat$attr$status
  ART.status <- dat$attr$ART.status
  stage <- dat$attr$stage
  stage.time <- dat$attr$stage.time
  ART.time <- dat$attr$ART.time


  ## Parameters ##
  n <- network.size(dat$nw)
  a.rate <- dat$param$arrival.rate

  
  #Arrival Process
  nArrivalsExp <- n * a.rate
  nArrivals <- rpois(1, nArrivalsExp)

  if (nArrivals > 0) {
    dat$nw <- add.vertices(dat$nw, nv = nArrivals)
    newNodes <- (n + 1):(n + nArrivals)
    dat$nw <- activate.vertices(dat$nw, onset = at, terminus = Inf,
                                v = newNodes)
  }

  
  # Add attributes for new arrivals
  if (nArrivals > 0) {
    active <- c(active, rep(1, nArrivals))
    status <- c(status, rep("s", nArrivals))
    ART.status <- c(ART.status, rep(0, nArrivals))
    stage <- c(stage, rep(NA, nArrivals))
    stage.time <- c(stage.time, rep(NA, nArrivals))
    ART.time <- c(ART.time, rep(NA, nArrivals))
  }

  
  ## UPDATE NODE ATTRIBUTES ##
  dat$attr$active <- active
  dat$attr$status <- status
  dat$attr$ART.status <- ART.status
  dat$attr$stage <- stage
  dat$attr$stage.time <- stage.time
  dat$attr$ART.time <- ART.time


  ## SUMMARY STATISTICS ##
  #Arrivals
  dat$epi$a.flow[at] <- nArrivals

  return(dat)
}
