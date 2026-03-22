
##
## HIV Transmission Model
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##
## Authors: Connor Van Meter, Samuel M. Jenness, Yuan Zhao, Emeli Anderson
## Date: March 2019
##


# Replacement infection/transmission module -------------------------------

infect <- function(dat, at) {

  ## Attributes ##
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  ART.status <- get_attr(dat, "ART.status")
  stage <- get_attr(dat, "stage")
  stage.time <- get_attr(dat, "stage.time")
  ART.time <- get_attr(dat, "ART.time")

  ## Parameters ##
  inf.prob.chronic <- get_param(dat, "inf.prob.chronic")
  relative.inf.prob.acute <- get_param(dat, "relative.inf.prob.acute")
  relative.inf.prob.AIDS <- get_param(dat, "relative.inf.prob.AIDS")
  relative.inf.prob.ART <- get_param(dat, "relative.inf.prob.ART")
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
    if (!is.null(del)) {

      # Transmission probability follows a multiplicative structure:
      #   base chronic rate * stage multiplier * ART multiplier
      # Acute and AIDS stages scale infectiousness relative to chronic;
      # ART independently reduces infectiousness across all stages.
      stage_mult <- ifelse(stage[del$inf] == "acute", relative.inf.prob.acute,
                   ifelse(stage[del$inf] == "AIDS", relative.inf.prob.AIDS, 1))
      art_mult <- ifelse(ART.status[del$inf] == 1, relative.inf.prob.ART, 1)
      del$transProb <- inf.prob.chronic * stage_mult * art_mult

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
        status[idsNewInf] <- "i"
        stage[idsNewInf] <- "acute"
        ART.status[idsNewInf] <- 0
        stage.time[idsNewInf] <- 0
        ART.time[idsNewInf] <- 0
      }
    }
  }

  ## Update attributes ##
  dat <- set_attr(dat, "status", status)
  dat <- set_attr(dat, "stage", stage)
  dat <- set_attr(dat, "ART.status", ART.status)
  dat <- set_attr(dat, "stage.time", stage.time)
  dat <- set_attr(dat, "ART.time", ART.time)

  ## Save summary statistics ##
  dat <- set_epi(dat, "acute.flow", at, nInf)
  dat <- set_epi(dat, "s.num", at, sum(active == 1 & status == "s"))
  dat <- set_epi(dat, "acute.ART.num", at,
                 sum(active == 1 & stage == "acute" &
                       ART.status == 1, na.rm = TRUE))
  dat <- set_epi(dat, "acute.NoART.num", at,
                 sum(active == 1 & stage == "acute" &
                       ART.status == 0, na.rm = TRUE))
  dat <- set_epi(dat, "chronic1.ART.num", at,
                 sum(active == 1 & stage == "chronic1" &
                       ART.status == 1, na.rm = TRUE))
  dat <- set_epi(dat, "chronic1.NoART.num", at,
                 sum(active == 1 & stage == "chronic1" &
                       ART.status == 0, na.rm = TRUE))
  dat <- set_epi(dat, "chronic2.ART.num", at,
                 sum(active == 1 & stage == "chronic2" &
                       ART.status == 1, na.rm = TRUE))
  dat <- set_epi(dat, "chronic2.NoART.num", at,
                 sum(active == 1 & stage == "chronic2" &
                       ART.status == 0, na.rm = TRUE))
  dat <- set_epi(dat, "AIDS.ART.num", at,
                 sum(active == 1 & stage == "AIDS" &
                       ART.status == 1, na.rm = TRUE))
  dat <- set_epi(dat, "AIDS.NoART.num", at,
                 sum(active == 1 & stage == "AIDS" &
                       ART.status == 0, na.rm = TRUE))

  return(dat)
}


# HIV progression module --------------------------------------------------
# (Replaces the recovery module)

progress <- function(dat, at) {

  ## Attributes ##
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")

  ## Initialize stage/ART tracking attributes at simulation start ##
  if (at == 2) {
    dat <- set_attr(dat, "stage",
                    ifelse(status == "i", "acute", NA))
    dat <- set_attr(dat, "ART.status",
                    ifelse(status == "i", 0, NA))
    dat <- set_attr(dat, "stage.time",
                    ifelse(status == "i", 0, NA))
    dat <- set_attr(dat, "ART.time",
                    ifelse(status == "i", 0, NA))
  }

  ## Attributes (read after initialization) ##
  ART.status <- get_attr(dat, "ART.status")
  stage <- get_attr(dat, "stage")
  ART.time <- get_attr(dat, "ART.time")
  stage.time <- get_attr(dat, "stage.time")

  ## Parameters ##
  AcuteToChronic1.Rate <- get_param(dat, "AcuteToChronic1.Rate")
  Chronic1ToChronic2.Rate <- get_param(dat, "Chronic1ToChronic2.Rate")
  Chronic2ToAIDS.Rate <- get_param(dat, "Chronic2ToAIDS.Rate")
  ART.Treatment.Rate <- get_param(dat, "ART.Treatment.Rate")
  ART.Discontinuance.Rate <- get_param(dat, "ART.Discontinuance.Rate")
  ART.Progression.Reduction.Rate <- get_param(dat, "ART.Progression.Reduction.Rate")


  ## Increment time-in-stage and time-on/off-ART counters ##
  ## The != 0 guard in eligibility checks below prevents individuals from
  ## transitioning in the same timestep they entered their current state.
  ART.time <- ifelse(!is.na(ART.time), ART.time + 1, ART.time)
  stage.time <- ifelse(!is.na(stage.time), stage.time + 1, stage.time)


  ## ---- ART Treatment ---- ##
  ## Eligible: infected, not on ART, not in first timestep of current state
  idsEligART <- which(active == 1 & status == "i" &
                        ART.status == 0 & ART.time != 0 & !is.na(ART.status))
  idsART <- idsEligART[rbinom(length(idsEligART), 1, ART.Treatment.Rate) == 1]

  if (length(idsART) > 0) {
    ART.status[idsART] <- 1
    ART.time[idsART] <- 0
  }


  ## ---- ART Discontinuance ---- ##
  ## Eligible: infected, on ART, not in first timestep of current state
  idsEligARTDisc <- which(active == 1 & status == "i" &
                            ART.status == 1 & ART.time != 0 & !is.na(ART.status))
  idsARTDisc <- idsEligARTDisc[rbinom(length(idsEligARTDisc), 1,
                                       ART.Discontinuance.Rate) == 1]

  if (length(idsARTDisc) > 0) {
    ART.status[idsARTDisc] <- 0
    ART.time[idsARTDisc] <- 0
  }


  ## Save ART flow statistics (before progression modifies stage) ##
  dat <- set_epi(dat, "acute.ART.treatment.flow", at,
                 sum(stage[idsART] == "acute"))
  dat <- set_epi(dat, "acute.ART.discontinuance.flow", at,
                 sum(stage[idsARTDisc] == "acute"))
  dat <- set_epi(dat, "chronic1.ART.treatment.flow", at,
                 sum(stage[idsART] == "chronic1"))
  dat <- set_epi(dat, "chronic1.ART.discont.flow", at,
                 sum(stage[idsARTDisc] == "chronic1"))
  dat <- set_epi(dat, "chronic2.ART.treatment.flow", at,
                 sum(stage[idsART] == "chronic2"))
  dat <- set_epi(dat, "chronic2.ART.discont.flow", at,
                 sum(stage[idsARTDisc] == "chronic2"))
  dat <- set_epi(dat, "AIDS.ART.treatment.flow", at,
                 sum(stage[idsART] == "AIDS"))
  dat <- set_epi(dat, "AIDS.ART.discontinuance.flow", at,
                 sum(stage[idsARTDisc] == "AIDS"))


  ## ---- Acute -> Chronic 1 Progression ---- ##
  ## ART reduces progression rate by ART.Progression.Reduction.Rate
  idsEligC1 <- which(active == 1 & stage == "acute" &
                       stage.time != 0 & !is.na(ART.status))
  prog.rate.C1 <- ifelse(ART.status[idsEligC1] == 1,
                          AcuteToChronic1.Rate * ART.Progression.Reduction.Rate,
                          AcuteToChronic1.Rate)
  idsChronic1 <- idsEligC1[rbinom(length(idsEligC1), 1, prog.rate.C1) == 1]

  if (length(idsChronic1) > 0) {
    stage[idsChronic1] <- "chronic1"
    stage.time[idsChronic1] <- 0
  }


  ## ---- Chronic 1 -> Chronic 2 Progression ---- ##
  idsEligC2 <- which(active == 1 & stage == "chronic1" &
                       stage.time != 0 & !is.na(ART.status))
  prog.rate.C2 <- ifelse(ART.status[idsEligC2] == 1,
                          Chronic1ToChronic2.Rate * ART.Progression.Reduction.Rate,
                          Chronic1ToChronic2.Rate)
  idsChronic2 <- idsEligC2[rbinom(length(idsEligC2), 1, prog.rate.C2) == 1]

  if (length(idsChronic2) > 0) {
    stage[idsChronic2] <- "chronic2"
    stage.time[idsChronic2] <- 0
  }


  ## ---- Chronic 2 -> AIDS Progression ---- ##
  idsEligAIDS <- which(active == 1 & stage == "chronic2" &
                          stage.time != 0 & !is.na(ART.status))
  prog.rate.AIDS <- ifelse(ART.status[idsEligAIDS] == 1,
                            Chronic2ToAIDS.Rate * ART.Progression.Reduction.Rate,
                            Chronic2ToAIDS.Rate)
  idsAIDS <- idsEligAIDS[rbinom(length(idsEligAIDS), 1, prog.rate.AIDS) == 1]

  if (length(idsAIDS) > 0) {
    stage[idsAIDS] <- "AIDS"
    stage.time[idsAIDS] <- 0
  }


  ## Save attributes ##
  dat <- set_attr(dat, "stage", stage)
  dat <- set_attr(dat, "stage.time", stage.time)
  dat <- set_attr(dat, "ART.status", ART.status)
  dat <- set_attr(dat, "ART.time", ART.time)

  ## Save progression flow statistics ##
  dat <- set_epi(dat, "chronic1.ART.flow", at,
                 sum(ART.status[idsChronic1] == 1))
  dat <- set_epi(dat, "chronic1.NoART.flow", at,
                 sum(ART.status[idsChronic1] == 0))
  dat <- set_epi(dat, "chronic2.ART.flow", at,
                 sum(ART.status[idsChronic2] == 1))
  dat <- set_epi(dat, "chronic2.NoART.flow", at,
                 sum(ART.status[idsChronic2] == 0))
  dat <- set_epi(dat, "AIDS.ART.flow", at,
                 sum(ART.status[idsAIDS] == 1))
  dat <- set_epi(dat, "AIDS.NoART.flow", at,
                 sum(ART.status[idsAIDS] == 0))

  return(dat)
}


# Departure module ---------------------------------------------------------

dfunc <- function(dat, at) {

  ## Attributes ##
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  ART.status <- get_attr(dat, "ART.status")
  stage <- get_attr(dat, "stage")
  stage.time <- get_attr(dat, "stage.time")
  exitTime <- get_attr(dat, "exitTime")

  ## Parameters ##
  departure.rate <- get_param(dat, "departure.rate")
  ART.Progression.Reduction.Rate <- get_param(dat, "ART.Progression.Reduction.Rate")
  AIDSToDepart.Rate <- get_param(dat, "AIDSToDepart.Rate")


  ## ---- Standard departure (susceptible + infected non-AIDS) ---- ##
  ## All individuals not in the AIDS stage depart at the background rate.
  idsEligStandard <- which(active == 1 & (is.na(stage) | stage != "AIDS"))
  idsDepart <- idsEligStandard[rbinom(length(idsEligStandard), 1,
                                       departure.rate) == 1]


  ## ---- AIDS departure on ART (reduced rate) ---- ##
  idsEligAIDSART <- which(active == 1 & stage == "AIDS" &
                            stage.time != 0 & ART.status == 1 & !is.na(ART.status))
  idsDepartAIDSART <- idsEligAIDSART[rbinom(length(idsEligAIDSART), 1,
                                              AIDSToDepart.Rate *
                                                ART.Progression.Reduction.Rate) == 1]


  ## ---- AIDS departure not on ART ---- ##
  idsEligAIDSNoART <- which(active == 1 & stage == "AIDS" &
                              stage.time != 0 & ART.status == 0 &
                              !is.na(ART.status))
  idsDepartAIDSNoART <- idsEligAIDSNoART[rbinom(length(idsEligAIDSNoART), 1,
                                                   AIDSToDepart.Rate) == 1]


  ## Combine all departures and deactivate ##
  idsDeparted <- c(idsDepart, idsDepartAIDSART, idsDepartAIDSNoART)
  if (length(idsDeparted) > 0) {
    active[idsDeparted] <- 0
    exitTime[idsDeparted] <- at
  }

  ## Save attributes ##
  dat <- set_attr(dat, "active", active)
  dat <- set_attr(dat, "exitTime", exitTime)

  ## Save summary statistics ##
  dat <- set_epi(dat, "depart.standard.ART.flow", at,
                 sum(ART.status[idsDepart] == 1, na.rm = TRUE))
  dat <- set_epi(dat, "depart.standard.NoART.flow", at,
                 sum(ART.status[idsDepart] == 0 |
                       is.na(ART.status[idsDepart])))
  dat <- set_epi(dat, "depart.AIDS.ART.flow", at,
                 length(idsDepartAIDSART))
  dat <- set_epi(dat, "depart.AIDS.NoART.flow", at,
                 length(idsDepartAIDSNoART))

  return(dat)
}


# Arrival module -----------------------------------------------------------

afunc <- function(dat, at) {

  ## Parameters ##
  n <- network.size(dat$run$nw[[1]])
  a.rate <- get_param(dat, "arrival.rate")

  ## Arrival process ##
  nArrivals <- rpois(1, n * a.rate)

  if (nArrivals > 0) {
    # Extend attribute vectors (must read before append_core_attr modifies dat)
    status <- c(get_attr(dat, "status"), rep("s", nArrivals))
    stage <- c(get_attr(dat, "stage"), rep(NA, nArrivals))
    stage.time <- c(get_attr(dat, "stage.time"), rep(NA, nArrivals))
    ART.status <- c(get_attr(dat, "ART.status"), rep(NA, nArrivals))
    ART.time <- c(get_attr(dat, "ART.time"), rep(NA, nArrivals))

    # Update core attributes (active, entrTime, exitTime)
    dat <- append_core_attr(dat, at, nArrivals)

    # Set extended attribute vectors
    dat <- set_attr(dat, "status", status)
    dat <- set_attr(dat, "stage", stage)
    dat <- set_attr(dat, "stage.time", stage.time)
    dat <- set_attr(dat, "ART.status", ART.status)
    dat <- set_attr(dat, "ART.time", ART.time)
  }

  ## Summary statistics ##
  dat <- set_epi(dat, "a.flow", at, nArrivals)

  return(dat)
}
