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
  Acute.Time <- dat$attr$Acute.Time
  ART.status <- dat$attr$ART.status
  ART.Discontinuance.Time <- dat$attr$ART.Discontinuance.Time


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
        idsAcuteNoART <- which(stage[del$inf] == "acute" & ART.status[del$inf] == 0)
        del$transProb[idsAcuteNoART] <- inf.prob.chronic * relative.inf.prob.acute

        #Trans Prob - Chronic, No ART
        idsChronicNoART <- which(stage[del$inf] %in% c("chronic1", "chronic2") &
                                   ART.status[del$inf] == 0)
        del$transProb[idsChronicNoART] <- inf.prob.chronic

        #Trans Prob - AIDS, No ART
        idsAIDSNoART <- which(stage[del$inf] == "AIDS" & ART.status[del$inf] == 0)
        del$transProb[idsAIDSNoART] <- inf.prob.chronic * relative.inf.prob.AIDS

        #Trans Prob - Acute, ART
        idsAcuteART <- which(stage[del$inf] == "acute" & ART.status[del$inf] == 1)
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
        Acute.Time[idsNewInf] <- at
        ART.status[idsNewInf] <- 0
        ART.Discontinuance.Time[idsNewInf] <- 0

      }
    }
  }

  ## Update attributes
  dat$attr$stage <- stage
  dat$attr$status <- status
  dat$attr$Acute.Time <- Acute.Time
  dat$attr$ART.status <- ART.status
  dat$attr$ART.Discontinuance.Time <- ART.Discontinuance.Time

  ## Save summary statistics
  dat$epi$si.flow[at] <- nInf
  dat$epi$acute.flow[at] <- nInf
  dat$epi$acute.num[at] <- sum(dat$attr$active == 1 &
                                 dat$attr$stage == "acute", na.rm = TRUE)
  dat$epi$acute.ART.num[at] <- sum(dat$attr$active == 1 &
                                     dat$attr$stage == "acute" &
                                     dat$attr$ART.status == 1 &
                                     !is.na(dat$attr$stage))
  dat$epi$acute.NoART.num[at] <- sum(dat$attr$active == 1 &
                                       dat$attr$stage == "acute" &
                                       dat$attr$ART.status == 0 &
                                       !is.na(dat$attr$stage))
  dat$epi$chronic1.num[at] <- sum(dat$attr$active == 1 &
                                    dat$attr$stage == "chronic1" &
                                    !is.na(dat$attr$stage))
  dat$epi$chronic1.ART.num[at] <- sum(dat$attr$active == 1 &
                                        dat$attr$stage == "chronic1" &
                                        dat$attr$ART.status == 1 &
                                        !is.na(dat$attr$stage))
  dat$epi$chronic1.NoART.num[at] <- sum(dat$attr$active == 1 &
                                          dat$attr$stage == "chronic1" &
                                          dat$attr$ART.status == 0 &
                                          !is.na(dat$attr$stage))
  dat$epi$chronic2.num[at] <- sum(dat$attr$active == 1
                                  & dat$attr$stage == "chronic2"
                                  & !is.na(dat$attr$stage))
  dat$epi$chronic2.ART.num[at] <- sum(dat$attr$active == 1 &
                                        dat$attr$stage == "chronic2" &
                                        dat$attr$ART.status == 1 &
                                        !is.na(dat$attr$stage))
  dat$epi$chronic2.NoART.num[at] <- sum(dat$attr$active == 1 &
                                          dat$attr$stage == "chronic2" &
                                          dat$attr$ART.status == 0 &
                                          !is.na(dat$attr$stage))
  dat$epi$AIDS.num[at] <- sum(dat$attr$active == 1 &
                                 dat$attr$stage == "AIDS" &
                                 !is.na(stage))
  dat$epi$AIDS.ART.num[at] <- sum(dat$attr$active == 1 &
                                     dat$attr$stage == "AIDS" &
                                     dat$attr$ART.status == 1 &
                                     !is.na(dat$attr$stage))
  dat$epi$AIDS.NoART.num[at] <- sum(dat$attr$active == 1 &
                                       dat$attr$stage == "AIDS" &
                                       dat$attr$ART.status == 0 &
                                       !is.na(dat$attr$stage))
  dat$epi$s.num[at] <- sum(dat$attr$active == 1 &
                             dat$attr$status == "s")
  dat$epi$i.num[at] <- sum(dat$attr$active == 1 &
                             dat$attr$status == "i")
  dat$epi$ART.num[at] <- sum(dat$attr$active == 1 &
                               dat$attr$ART.status == 1 &
                               !is.na(dat$attr$ART.status))
  dat$epi$NoART.num[at] <- sum(dat$attr$active == 1 &
                                 dat$attr$ART.status == 0 &
                                 !is.na(dat$attr$ART.status))

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
    dat$attr$ART.Discontinuance.Time <- ifelse(!is.na(dat$attr$ART.status) &
                                                 dat$attr$ART.status  == 0,
                                               0, NA)
    dat$attr$ART.Treatment.Time <- rep(NA, length(active))
    dat$attr$Acute.Time <- ifelse(dat$attr$status == "i",1,dat$attr$infTime)
    dat$attr$Chronic1.Time <- rep(NA, length(active))
    dat$attr$Chronic2.Time <- rep(NA, length(active))
    dat$attr$AIDS.Time <- rep(NA, length(active))

  }


  ## Attributes ##
  ART.status <- dat$attr$ART.status
  stage <- dat$attr$stage
  ART.Treatment.Time <- dat$attr$ART.Treatment.Time
  ART.Discontinuance.Time <- dat$attr$ART.Discontinuance.Time
  Acute.Time <- dat$attr$Acute.Time
  Chronic1.Time <- dat$attr$Chronic1.Time
  Chronic2.Time <- dat$attr$Chronic2.Time
  AIDS.Time <- dat$attr$AIDS.Time


  ## Progression Parameters ##
  AcuteToChronic1.Rate <- dat$param$AcuteToChronic1.Rate
  Chronic1ToChronic2.Rate <- dat$param$Chronic1ToChronic2.Rate
  Chronic2ToAIDS.Rate <- dat$param$Chronic2ToAIDS.Rate
  ART.Treatment.Rate <- dat$param$ART.Treatment.Rate
  ART.Discontinuance.Rate <- dat$param$ART.Discontinuance.Rate
  ART.Progression.Reduction.Rate <- dat$param$ART.Progression.Reduction.Rate


  #Trackers to calculate ART Treatment Flows
  acute.ART.count <- sum(active == 1 & ART.status == 1 &
                           stage == "acute" & !is.na(ART.status))
  chronic1.ART.count <- sum(active == 1 & ART.status == 1 &
                              stage == "chronic1" & !is.na(ART.status))
  chronic2.ART.count <- sum(active == 1 & ART.status == 1 &
                              stage == "chronic2" & !is.na(ART.status))
  AIDS.ART.count <- sum(active == 1 & ART.status == 1 &
                           stage == "AIDS" & !is.na(ART.status))
  ART.count <- sum(active == 1 & ART.status == 1 & !is.na(ART.status))


  ## ART Treatment ##
  nART = 0
  idsEligART <- which(active == 1 & status == "i" &
                        ART.Discontinuance.Time < at &
                        ART.status == 0 & !is.na(ART.status))
  nEligART <- length(idsEligART)

  if (nEligART > 0) {
    vecART <- which(rbinom(nEligART, 1, ART.Treatment.Rate) == 1)
    if (length(vecART) > 0) {
      idsART <- idsEligART[vecART]
      nART  <- length(idsART)
      ART.status[idsART] <- 1
      ART.Treatment.Time[idsART] <- at
    }
  }
  dat$attr$ART.status <- ART.status
  dat$attr$ART.Treatment.Time <- ART.Treatment.Time
  dat$epi$ART.flow[at] <- nART


  ## ART Discontinuance ##
  nARTDisc = 0
  idsEligARTDisc <- which(active == 1 & status == "i" &
                            ART.Treatment.Time < at &
                            ART.status == 1 & !is.na(ART.status))
  nEligARTDisc <- length(idsEligARTDisc)

  if (nEligARTDisc > 0) {
    vecARTDisc <- which(rbinom(nEligARTDisc, 1, ART.Discontinuance.Rate) == 1)
    if (length(vecARTDisc) > 0) {
      idsARTDisc <- idsEligARTDisc[vecARTDisc]
      nARTDisc  <- length(idsARTDisc)
      ART.status[idsARTDisc] <- 0
      ART.Discontinuance.Time[idsARTDisc] <- at
    }
  }
  dat$attr$ART.status <- ART.status
  dat$attr$ART.Discontinuance.Time <- ART.Discontinuance.Time
  dat$epi$NoART.flow[at] <- nARTDisc

  ##Net ART Treatment/Discontinuance Flows within HIV Status Subcompartment
  #acute.ART.flow <- intersect(idsART, which(stage == "acute"))


  dat$epi$acute.ART.net.flow[at] <- sum(active == 1 & ART.status == 1 &
                                          stage == "acute" &
                                          !is.na(ART.status)) -
                                          acute.ART.count
  dat$epi$chronic1.ART.net.flow[at] <- sum(active == 1 & ART.status == 1 &
                                             stage == "chronic1" &
                                             !is.na(ART.status)) -
                                             chronic1.ART.count
  dat$epi$chronic2.ART.net.flow[at] <- sum(active == 1 & ART.status == 1 &
                                             stage == "chronic2" &
                                             !is.na(ART.status)) -
                                             chronic2.ART.count
  dat$epi$AIDS.ART.net.flow[at] <- sum(active == 1 & ART.status == 1 &
                                          stage == "AIDS" &
                                          !is.na(ART.status)) -
                                          AIDS.ART.count
  dat$epi$ART.net.flow[at] <- sum(active == 1 & ART.status == 1 &
                                    !is.na(ART.status)) - ART.count


  ## Acute to chronic 1 stage progression process ##
  nEligChronic1ART <- 0
  nEligChronic1NoART <- 0
  nChronic1ART <- 0
  nChronic1NoART <- 0
  idsEligChronic1ART <- which(active == 1 & stage == "acute" &
                                Acute.Time < at & ART.status == 1 &
                                !is.na(ART.status))
  idsEligChronic1NoART <- which(active == 1 & stage == "acute" &
                                  Acute.Time < at & ART.status == 0 &
                                  !is.na(ART.status))
  nEligChronic1ART <- length(idsEligChronic1ART)
  nEligChronic1NoART <- length(idsEligChronic1NoART)

  if (nEligChronic1ART > 0) {
    vecChronic1ART <- which(rbinom(nEligChronic1ART, 1, AcuteToChronic1.Rate *
                                     ART.Progression.Reduction.Rate) == 1)

    if (length(vecChronic1ART) > 0) {
      idsChronic1ART <- idsEligChronic1ART[vecChronic1ART]
      nChronic1ART  <- length(idsChronic1ART)
      stage[idsChronic1ART] <- "chronic1"
      Chronic1.Time[idsChronic1ART] <- at
    }
  }

  if (nEligChronic1NoART > 0) {
    vecChronic1NoART <- which(rbinom(nEligChronic1NoART, 1,
                                     AcuteToChronic1.Rate) == 1)

    if (length(vecChronic1NoART) > 0) {
      idsChronic1NoART <- idsEligChronic1NoART[vecChronic1NoART]
      nChronic1NoART  <- length(idsChronic1NoART)
      stage[idsChronic1NoART] <- "chronic1"
      Chronic1.Time[idsChronic1NoART] <- at
    }
  }
  dat$attr$stage <- stage
  dat$attr$Chronic1.Time <- Chronic1.Time
  dat$epi$chronic1.flow[at] <- nChronic1ART + nChronic1NoART
  dat$epi$chronic1.ART.flow[at] <- nChronic1ART
  dat$epi$chronic1.NoART.flow[at] <- nChronic1NoART


  ## Chronic 1 to chronic 2 stage progression process ##
  nEligChronic2ART <- 0
  nEligChronic2NoART <- 0
  nChronic2ART <- 0
  nChronic2NoART <- 0
  idsEligChronic2ART <- which(active == 1 & stage == "chronic1" &
                                Chronic1.Time < at & ART.status == 1 &
                                !is.na(ART.status))
  idsEligChronic2NoART <- which(active == 1 & stage == "chronic1" &
                                  Chronic1.Time < at & ART.status == 0 &
                                  !is.na(ART.status))
  nEligChronic2ART <- length(idsEligChronic2ART)
  nEligChronic2NoART <- length(idsEligChronic2NoART)

  if (nEligChronic2ART > 0) {
    vecChronic2ART <- which(rbinom(nEligChronic2ART, 1,
                                   Chronic1ToChronic2.Rate *
                                     ART.Progression.Reduction.Rate) == 1)

    if (length(vecChronic2ART) > 0) {
      idsChronic2ART <- idsEligChronic2ART[vecChronic2ART]
      nChronic2ART  <- length(idsChronic2ART)
      stage[idsChronic2ART] <- "chronic2"
      Chronic2.Time[idsChronic2ART] <- at
    }
  }

  if (nEligChronic2NoART > 0) {
    vecChronic2NoART <- which(rbinom(nEligChronic2NoART, 1,
                                     Chronic1ToChronic2.Rate) == 1)

    if (length(vecChronic2NoART) > 0) {
      idsChronic2NoART <- idsEligChronic2NoART[vecChronic2NoART]
      nChronic2NoART  <- length(idsChronic2NoART)
      stage[idsChronic2NoART] <- "chronic2"
      Chronic2.Time[idsChronic2NoART] <- at
    }
  }
  dat$attr$stage <- stage
  dat$attr$Chronic2.Time <- Chronic2.Time
  dat$epi$chronic2.flow[at] <- nChronic2ART + nChronic2NoART
  dat$epi$chronic2.ART.flow[at] <- nChronic2ART
  dat$epi$chronic2.NoART.flow[at] <- nChronic2NoART


  ## Chronic 2 to AIDS stage progression process ##
  nEligAIDSNoART <- 0
  nAIDSART <- 0
  nAIDSNoART <- 0
  idsEligAIDSART <- which(active == 1 & stage == "chronic2" &
                             Chronic2.Time < at & ART.status == 1 &
                             !is.na(ART.status))
  idsEligAIDSNoART <- which(active == 1 & stage == "chronic2" &
                               Chronic2.Time < at & ART.status == 0 &
                               !is.na(ART.status))
  nEligAIDSART <- length(idsEligAIDSART)
  nEligAIDSNoART <- length(idsEligAIDSNoART)

  if (length(idsEligAIDSART) > 0) {
    vecAIDSART <- which(rbinom(nEligAIDSART, 1, Chronic2ToAIDS.Rate *
                                  ART.Progression.Reduction.Rate) == 1)

    if (length(vecAIDSART) > 0) {
      idsAIDSART <- idsEligAIDSART[vecAIDSART]
      nAIDSART  <- length(idsAIDSART)
      stage[idsAIDSART] <- "AIDS"
      AIDS.Time[idsAIDSART] <- at
    }
  }

  if (nEligAIDSNoART > 0) {
    vecAIDSNoART <- which(rbinom(nEligAIDSNoART, 1,
                                  Chronic2ToAIDS.Rate) == 1)

    if (length(vecAIDSNoART) > 0) {
      idsAIDSNoART <- idsEligAIDSNoART[vecAIDSNoART]
      nAIDSNoART  <- length(idsAIDSNoART)
      stage[idsAIDSNoART] <- "AIDS"
      AIDS.Time[idsAIDSNoART] <- at
    }
  }
  dat$attr$stage <- stage
  dat$attr$AIDS.Time <- AIDS.Time
  dat$epi$AIDS.flow[at] <- nAIDSART + nAIDSNoART
  dat$epi$AIDS.ART.flow[at] <- nAIDSART
  dat$epi$AIDS.NoART.flow[at] <- nAIDSNoART

  return(dat)
}


dfunc <- function(dat, at) {

  ## Attributes ##
  active <- dat$attr$active
  status <- dat$attr$status
  ART.status <- dat$attr$ART.status
  stage <- dat$attr$stage
  AIDS.Time <- dat$attr$AIDS.Time

  ## Parameters ##
  Departure.rates <- rep(dat$param$departure.rate, network.size(dat$nw))
  Departure.dis.mult <- dat$param$departure.disease.mult
  ART.Progression.Reduction.Rate <- dat$param$ART.Progression.Reduction.Rate
  AIDSToDepart.Rate <- dat$param$AIDSToDepart.Rate

  ## Query active individuals not in the AIDS stage of HIV ##
  idsEligDepartStandard <- which(active == 1 & (stage != "AIDS" &
                                                  !is.na(stage) |
                                                  status == "s"))
  nEligDepartStandard <- length(idsEligDepartStandard)
  nDeparturesStandard <- 0
  nDeparturesStandardART <- 0
  nDeparturesStandardNoART <- 0

  ## Establish empty vectors for departure IDs
  idsDeparture <- vector()
  idsDepartAIDSART <- vector()
  idsDepartAIDSNoART <- vector()


  ## Departure process for individuals not in the AIDS stage of HIV
  if (nEligDepartStandard > 0) {

    Departure_rates_of_elig <- Departure.rates[idsEligDepartStandard]

    ## Multiply Departure rates for diseased persons
    idsEligDepartStandard.inf <- which(status[idsEligDepartStandard] == "i")
    Departure_rates_of_elig[idsEligDepartStandard.inf] <-
      Departure.rates[idsEligDepartStandard.inf] * Departure.dis.mult

    ## Simulate Departure process
    vecDeparture <- which(rbinom(idsEligDepartStandard, 1,
                                 Departure_rates_of_elig) == 1)
    idsDeparture <- idsEligDepartStandard[vecDeparture]
    idsDepartureART <- which(ART.status[idsDeparture] == 1)
    idsDepartureNoART <- which(ART.status[idsDeparture] == 0
                               | is.na(ART.status[idsDeparture]))
    nDeparturesStandard <- length(idsDeparture)
    nDeparturesStandardART <- length(idsDepartureART)
    nDeparturesStandardNoART <- length(idsDepartureNoART)

  }

  ## Query active individuals in the AIDS stage of HIV ##
  idsEligDepartAIDS <- which(active == 1 & stage == "AIDS" &
                                AIDS.Time < at & !is.na(stage))
  nEligDepartAIDS <- length(idsEligDepartAIDS)
  nDeparturesAIDSART <- 0
  nDeparturesAIDSNoART <- 0

  ## Departure process for individuals in the AIDS stage of HIV
  if (nEligDepartAIDS > 0) {

    idsEligDepartAIDSART <- which(active == 1 & stage == "AIDS" &
                                     AIDS.Time < at & ART.status == 1 &
                                     !is.na(ART.status))
    idsEligDepartAIDSNoART <- which(active == 1 & stage == "AIDS" &
                                       AIDS.Time < at & ART.status == 0 &
                                       !is.na(ART.status))
    nEligDepartAIDSART <- length(idsEligDepartAIDSART)
    nEligDepartAIDSNoART <- length(idsEligDepartAIDSNoART)

    if (nEligDepartAIDSART > 0) {
      vecDepartAIDSART <- which(rbinom(nEligDepartAIDSART, 1,
                                        AIDSToDepart.Rate *
                                          ART.Progression.Reduction.Rate) == 1)

      if (length(vecDepartAIDSART) > 0) {
        idsDepartAIDSART <- idsEligDepartAIDSART[vecDepartAIDSART]
        nDeparturesAIDSART  <- length(idsDepartAIDSART)
        active[idsDepartAIDSART] <- 0
      }
    }

    if (nEligDepartAIDSNoART > 0) {
      vecDepartAIDSNoART <- which(rbinom(nEligDepartAIDSNoART, 1,
                                          AIDSToDepart.Rate) == 1)

      if (length(vecDepartAIDSNoART) > 0) {
        idsDepartAIDSNoART <- idsEligDepartAIDSNoART[vecDepartAIDSNoART]
        nDeparturesAIDSNoART  <- length(idsDepartAIDSNoART)
        active[idsDepartAIDSNoART] <- 0
      }
    }
  }

  ##Save departure summary statistics
  dat$epi$d.flow[at] <- nDeparturesStandard + nDeparturesAIDSART +
                        nDeparturesAIDSNoART
  dat$epi$depart.standard.flow[at] <- nDeparturesStandard
  dat$epi$depart.standard.ART.flow[at] <- nDeparturesStandardART
  dat$epi$depart.standard.NoART.flow[at] <- nDeparturesStandardNoART
  dat$epi$depart.AIDS.flow[at] <- nDeparturesAIDSART + nDeparturesAIDSNoART
  dat$epi$depart.AIDS.ART.flow[at] <- nDeparturesAIDSART
  dat$epi$depart.AIDS.NoART.flow[at] <- nDeparturesAIDSNoART

  ## Update nodal attributes on attr and networkDynamic object ##
  if (nDeparturesStandard > 0 | nDeparturesAIDSART > 0 |
      nDeparturesAIDSNoART > 0) {
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
  infTime <- dat$attr$infTime
  entrTime <- dat$attr$entrTime
  exitTime <- dat$attr$exitTime
  ART.status <- dat$attr$ART.status
  stage <- dat$attr$stage
  Acute.Time <- dat$attr$Acute.Time
  Chronic1.Time <- dat$attr$Chronic1.Time
  Chronic2.Time <- dat$attr$Chronic2.Time
  AIDS.Time <- dat$attr$AIDS.Time
  ART.Treatment.Time <- dat$attr$ART.Treatment.Time
  ART.Discontinuance.Time <- dat$attr$ART.Discontinuance.Time


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
    infTime <- c(infTime, rep(NA, nArrivals))
    entrTime <- c(entrTime, rep(at, nArrivals))
    exitTime <- c(exitTime, rep(NA, nArrivals))
    ART.status <- c(ART.status, rep(0, nArrivals))
    stage <- c(stage, rep(NA, nArrivals))
    Acute.Time <- c(Acute.Time, rep(NA, nArrivals))
    Chronic1.Time <- c(Chronic1.Time, rep(NA, nArrivals))
    Chronic2.Time <- c(Chronic2.Time, rep(NA, nArrivals))
    AIDS.Time <- c(AIDS.Time, rep(NA, nArrivals))
    ART.Treatment.Time <- c(ART.Treatment.Time, rep(NA, nArrivals))
    ART.Discontinuance.Time <- c(ART.Discontinuance.Time, rep(NA, nArrivals))
  }

  ## UPDATE NODE ATTRIBUTES ##
  dat$attr$active <- active
  dat$attr$infTime <- infTime
  dat$attr$entrTime <- entrTime
  dat$attr$exitTime <- exitTime
  dat$attr$status <- status
  dat$attr$ART.status <- ART.status
  dat$attr$stage <- stage
  dat$attr$Acute.Time <- Acute.Time
  dat$attr$Chronic1.Time <- Chronic1.Time
  dat$attr$Chronic2.Time <- Chronic2.Time
  dat$attr$AIDS.Time <- AIDS.Time
  dat$attr$ART.Treatment.Time <- ART.Treatment.Time
  dat$attr$ART.Discontinuance.Time <- ART.Discontinuance.Time


  ## SUMMARY STATISTICS ##

  #Arrivals
  dat$epi$a.flow[at] <- nArrivals

  return(dat)
}
