##
## HIV Progression Module (One-Mode)
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##
## Authors: Samuel M. Jenness, Yuan Zhao, Connor Van Meter
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
  HIV.status <- dat$attr$HIV.status
  Acute.Time <- dat$attr$Acute.Time
  ART.status <- dat$attr$ART.status
  ART.Discontinuance.Time <- dat$attr$ART.Discontinuance.Time


  ## Parameters ##
  inf.prob.chronic <- dat$param$inf.prob.chronic
  relative.inf.prob.acute <- dat$param$relative.inf.prob.acute
  relative.inf.prob.final <- dat$param$relative.inf.prob.final
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
        idsAcuteNoART <- which(HIV.status[del$inf] == "acute"
                               & ART.status[del$inf] == 0)
        del$transProb[idsAcuteNoART] <- inf.prob.chronic * relative.inf.prob.acute

        #Trans Prob - Chronic, No ART
        idsChronicNoART <- which(HIV.status[del$inf] %in% c("chronic1","chronic2")
                                 & ART.status[del$inf] == 0)
        del$transProb[idsChronicNoART] <- inf.prob.chronic

        #Trans Prob - Final, No ART
        idsFinalNoART <- which(HIV.status[del$inf] == "final"
                               & ART.status[del$inf] == 0)
        del$transProb[idsFinalNoART] <- inf.prob.chronic * relative.inf.prob.final

        #Trans Prob - Acute, ART
        idsAcuteART <- which(HIV.status[del$inf] == "acute"
                             & ART.status[del$inf] == 1)
        del$transProb[idsAcuteART] <- inf.prob.chronic * relative.inf.prob.acute * relative.inf.prob.ART

        #Trans Prob - Chronic, ART
        idsChronicART <- which(HIV.status[del$inf] %in% c("chronic1","chronic2")
                               & ART.status[del$inf] == 1)
        del$transProb[idsChronicART] <- inf.prob.chronic * relative.inf.prob.ART

        #Trans Prob - Final, ART
        idsFinalART <- which(HIV.status[del$inf] == "final"
                             & ART.status[del$inf] == 1)
        del$transProb[idsFinalART] <- inf.prob.chronic * relative.inf.prob.final * relative.inf.prob.ART

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
        HIV.status[idsNewInf] <- "acute"
        status[idsNewInf] <- "i"
        Acute.Time[idsNewInf] <- at
        ART.status[idsNewInf] <- 0
        ART.Discontinuance.Time[idsNewInf] <- 0

      }

    }

  }

  ## Update attributes
  dat$epi$Acute.flow[at] <- nInf
  dat$attr$HIV.status <- HIV.status
  dat$attr$status <- status
  dat$attr$Acute.Time <- Acute.Time
  dat$attr$ART.status <- ART.status
  dat$attr$ART.Discontinuance.Time <- ART.Discontinuance.Time

  ## Save summary statistics
  dat$epi$acute.num[at] <- sum(dat$attr$active == 1 & dat$attr$HIV.status == "acute")

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

    dat$attr$HIV.status <- ifelse(dat$attr$status == "i","acute",NA)
    dat$attr$ART.status <- ifelse(dat$attr$status == "i",0,NA)

    dat$attr$ART.Treatment.Time <- rep(NA, length(active))
    dat$attr$ART.Discontinuance.Time <- rep(NA, length(active))
    dat$attr$Acute.Time <- ifelse(dat$attr$status == "i",1,dat$attr$infTime)
    dat$attr$Chronic1.Time <- rep(NA, length(active))
    dat$attr$Chronic2.Time <- rep(NA, length(active))
    dat$attr$Final.Time <- rep(NA, length(active))

  }


  ## Attributes ##
  ART.status <- dat$attr$ART.status
  HIV.status <- dat$attr$HIV.status
  ART.Treatment.Time <- dat$attr$ART.Treatment.Time
  ART.Discontinuance.Time <- dat$attr$ART.Discontinuance.Time
  Acute.Time <- dat$attr$Acute.Time
  Chronic1.Time <- dat$attr$Chronic1.Time
  Chronic2.Time <- dat$attr$Chronic2.Time
  Final.Time <- dat$attr$Final.Time


  ## Progression Parameters ##
  AcuteToChronic1.Rate <- dat$param$AcuteToChronic1.Rate
  Chronic1ToChronic2.Rate <- dat$param$Chronic1ToChronic2.Rate
  Chronic2ToFinal.Rate <- dat$param$Chronic2ToFinal.Rate
  ART.Treatment.Rate <- dat$param$ART.Treatment.Rate
  ART.Discontinuance.Rate <- dat$param$ART.Discontinuance.Rate
  ART.Progression.Reduction.Rate <- dat$param$ART.Progression.Reduction.Rate


  ## ART Treatment ##
  nART = 0
  idsEligART <- which(active == 1 & status == "i" & ART.Discontinuance.Time < at & ART.status == 0)
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


  ## ART Discontinuance ##
  nARTDisc = 0
  idsEligARTDisc <- which(active == 1 & status == "i" & ART.Treatment.Time < at & ART.status == 1)
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


  ## Acute to chronic 1 stage progression process ##
  nEligChronic1ART <- 0
  nEligChronic1NoART <- 0
  nChronic1ART <- 0
  nChronic1NoART <- 0
  idsEligChronic1ART <- which(active == 1 & HIV.status == "acute" & Acute.Time < at & ART.status == 1)
  idsEligChronic1NoART <- which(active == 1 & HIV.status == "acute" & Acute.Time < at & ART.status == 0)
  nEligChronic1ART <- length(idsEligChronic1ART)
  nEligChronic1NoART <- length(idsEligChronic1NoART)

  if (nEligChronic1ART > 0) {
    vecChronic1ART <- which(rbinom(nEligChronic1ART, 1, AcuteToChronic1.Rate * ART.Progression.Reduction.Rate) == 1)

    if (length(vecChronic1ART) > 0) {
      idsChronic1ART <- idsEligChronic1ART[vecChronic1ART]
      nChronic1ART  <- length(idsChronic1ART)
      HIV.status[idsChronic1ART] <- "chronic1"
      Chronic1.Time[idsChronic1ART] <- at
    }
  }

  if (nEligChronic1NoART > 0) {
    vecChronic1NoART <- which(rbinom(nEligChronic1NoART, 1, AcuteToChronic1.Rate) == 1)

    if (length(vecChronic1NoART) > 0) {
      idsChronic1NoART <- idsEligChronic1NoART[vecChronic1NoART]
      nChronic1NoART  <- length(idsChronic1NoART)
      HIV.status[idsChronic1NoART] <- "chronic1"
      Chronic1.Time[idsChronic1NoART] <- at
    }
  }
  dat$attr$HIV.status <- HIV.status
  dat$attr$Chronic1.Time <- Chronic1.Time
  dat$epi$Chronic1.flow <- nChronic1ART + nChronic1NoART
  dat$epi$Chronic1.ART.flow <- nChronic1ART
  dat$epi$Chronic1.NoART.flow <- nChronic1NoART


  ## Chronic 1 to chronic 2 stage progression process ##
  nEligChronic2ART <- 0
  nEligChronic2NoART <- 0
  nChronic2ART <- 0
  nChronic2NoART <- 0
  idsEligChronic2ART <- which(active == 1 & HIV.status == "chronic1" & Chronic1.Time < at & ART.status == 1)
  idsEligChronic2NoART <- which(active == 1 & HIV.status == "chronic1" & Chronic1.Time < at & ART.status == 0)
  nEligChronic2ART <- length(idsEligChronic2ART)
  nEligChronic2NoART <- length(idsEligChronic2NoART)

  if (nEligChronic2ART > 0) {
    vecChronic2ART <- which(rbinom(nEligChronic2ART, 1, Chronic1ToChronic2.Rate * ART.Progression.Reduction.Rate) == 1)

    if (length(vecChronic2ART) > 0) {
      idsChronic2ART <- idsEligChronic2ART[vecChronic2ART]
      nChronic2ART  <- length(idsChronic2ART)
      HIV.status[idsChronic2ART] <- "chronic2"
      Chronic2.Time[idsChronic2ART] <- at
    }
  }

  if (nEligChronic2NoART > 0) {
    vecChronic2NoART <- which(rbinom(nEligChronic2NoART, 1, Chronic1ToChronic2.Rate) == 1)

    if (length(vecChronic2NoART) > 0) {
      idsChronic2NoART <- idsEligChronic2NoART[vecChronic2NoART]
      nChronic2NoART  <- length(idsChronic2NoART)
      HIV.status[idsChronic2NoART] <- "chronic2"
      Chronic2.Time[idsChronic2NoART] <- at
    }
  }
  dat$attr$HIV.status <- HIV.status
  dat$attr$Chronic2.Time <- Chronic2.Time
  dat$epi$Chronic2.flow <- nChronic2ART + nChronic2NoART
  dat$epi$Chronic2.ART.flow <- nChronic2ART
  dat$epi$Chronic2.NoART.flow <- nChronic2NoART


  ## Chronic 2 to Final stage progression process ##
  nEligFinalART <- 0
  nEligFinalNoART <- 0
  nFinalART <- 0
  nFinalNoART <- 0
  idsEligFinalART <- which(active == 1 & HIV.status == "chronic2" & Chronic2.Time < at & ART.status == 1)
  idsEligFinalNoART <- which(active == 1 & HIV.status == "chronic2" & Chronic2.Time < at & ART.status == 0)
  nEligFinalART <- length(idsEligFinalART)
  nEligFinalNoART <- length(idsEligFinalNoART)

  if (nEligFinalART > 0) {
    vecFinalART <- which(rbinom(nEligFinalART, 1, Chronic2ToFinal.Rate * ART.Progression.Reduction.Rate) == 1)

    if (length(vecFinalART) > 0) {
      idsFinalART <- idsEligFinalART[vecFinalART]
      nFinalART  <- length(idsFinalART)
      HIV.status[idsFinalART] <- "final"
      Final.Time[idsFinalART] <- at
    }
  }

  if (nEligFinalNoART > 0) {
    vecFinalNoART <- which(rbinom(nEligFinalNoART, 1, Chronic2ToFinal.Rate) == 1)

    if (length(vecFinalNoART) > 0) {
      idsFinalNoART <- idsEligFinalNoART[vecFinalNoART]
      nFinalNoART  <- length(idsFinalNoART)
      HIV.status[idsFinalNoART] <- "final"
      Final.Time[idsFinalNoART] <- at
    }
  }
  dat$attr$HIV.status <- HIV.status
  dat$attr$Final.Time <- Final.Time
  dat$epi$Final.flow <- nFinalART + nFinalNoART
  dat$epi$Final.ART.flow <- nFinalART
  dat$epi$Final.NoART.flow <- nFinalNoART


  ## Save summary statistics ##
  dat$epi$chronic1.num[at] <- sum(dat$attr$active == 1 & dat$attr$HIV.status == "chronic1")
  dat$epi$chronic2.num[at] <- sum(dat$attr$active == 1 & dat$attr$HIV.status == "chronic2")


  return(dat)
}


dfunc <- function(dat, at) {

  ## Attributes ##
  active <- dat$attr$active
  status <- dat$attr$status
  ART.status <- dat$attr$ART.status
  HIV.status <- dat$attr$HIV.status
  Final.Time <- dat$attr$Final.Time

  ## Parameters ##
  Departure.rates <- rep(dat$param$departure.rate, network.size(dat$nw))
  Departure.dis.mult <- dat$param$departure.disease.mult
  ART.Progression.Reduction.Rate <- dat$param$ART.Progression.Reduction.Rate
  FinalToEnd.Rate <- dat$param$FinalToEnd.Rate

  ## Query active individuals not in the final stage of HIV ##
  idsEligEndStandard <- which(active == 1 & HIV.status != "final")
  nEligStandard <- length(idsEligEndStandard)
  nDeparturesStandard <- 0
  nEndART <- 0
  nEndNoART <- 0

  if (nEligStandard > 0) {

    Departure_rates_of_elig <- Departure.rates[idsEligEndStandard]

    ## Multiply Departure rates for diseased persons
    idsEligEndStandard.inf <- which(status[idsEligEndStandard] == "i")
    Departure_rates_of_elig[idsEligEndStandard.inf] <- Departure.rates[idsEligEndStandard.inf] *
      Departure.dis.mult

    ## Simulate Departure process
    vecDeparture <- which(rbinom(idsEligEndStandard, 1, Departure_rates_of_elig) == 1)
    idsDeparture <- idsEligEndStandard[vecDeparture]
    nDeparturesStandard <- length(idsDeparture)


    ## Final to End stage progression process ##
    nEligEndART <- 0
    nEligEndNoART <- 0
    idsEligEndART <- which(active == 1 & HIV.status == "final" & Final.Time < at & ART.status == 1)
    idsEligEndNoART <- which(active == 1 & HIV.status == "final" & Final.Time < at & ART.status == 0)
    nEligEndART <- length(idsEligEndART)
    nEligEndNoART <- length(idsEligEndNoART)

    if (nEligEndART > 0) {
      vecEndART <- which(rbinom(nEligEndART, 1, FinalToEnd.Rate * ART.Progression.Reduction.Rate) == 1)

      if (length(vecEndART) > 0) {
        idsEndART <- idsEligEndART[vecEndART]
        nEndART  <- length(idsEndART)
        active[idsEndART] <- 0
      }
    }

    if (nEligEndNoART > 0) {
      vecEndNoART <- which(rbinom(nEligEndNoART, 1, FinalToEnd.Rate) == 1)

      if (length(vecEndNoART) > 0) {
        idsEndNoART <- idsEligEndNoART[vecEndNoART]
        nEndNoART  <- length(idsEndNoART)
        active[idsEndNoART] <- 0
      }
    }
    dat$epi$End.flow <- nEndART + nEndNoART
    dat$epi$End.ART.flow <- nEndART
    dat$epi$End.NoART.flow <- nEndNoART

    ## Update nodal attributes on attr and networkDynamic object ##

    #Vertices to remove
    idsDeparted <- c(idsDeparture, idsEligEndART, idsEligEndNoART)

    if (nDeparturesStandard > 0 | nEndART > 0 | nEndNoART > 0) {
      active[idsDeparted] <- 0
      dat$nw <- deactivate.vertices(dat$nw, onset = at, terminus = Inf,
                                    v = idsDeparted, deactivate.edges = TRUE)
    }

    ## Write out updated status attribute ##
    dat$attr$active <- active
  }

  ## Summary statistics ##
  dat$epi$d.flow[at] <- nDeparturesStandard
  dat$epi$final.num[at] <- sum(dat$attr$active == 1 & dat$attr$HIV.status == "final")
  dat$epi$end.num[at] <- sum(dat$attr$active == 0 & dat$attr$HIV.status == "final")

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
  HIV.status <- dat$attr$HIV.status
  Acute.Time <- dat$attr$Acute.Time
  Chronic1.Time <- dat$attr$Chronic1.Time
  Chronic2.Time <- dat$attr$Chronic2.Time
  Final.Time <- dat$attr$Final.Time
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
    HIV.status <- c(HIV.status, rep(NA, nArrivals))
    Acute.Time <- c(Acute.Time, rep(NA, nArrivals))
    Chronic1.Time <- c(Chronic1.Time, rep(NA, nArrivals))
    Chronic2.Time <- c(Chronic2.Time, rep(NA, nArrivals))
    Final.Time <- c(Final.Time, rep(NA, nArrivals))
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
  dat$attr$HIV.status <- HIV.status
  dat$attr$Acute.Time <- Acute.Time
  dat$attr$Chronic1.Time <- Chronic1.Time
  dat$attr$Chronic2.Time <- Chronic2.Time
  dat$attr$Final.Time <- Final.Time
  dat$attr$ART.Treatment.Time <- ART.Treatment.Time
  dat$attr$ART.Discontinuance.Time <- ART.Discontinuance.Time


  ## SUMMARY STATISTICS ##

  #Arrivals
  dat$epi$a.flow[at] <- nArrivals

  return(dat)
}
