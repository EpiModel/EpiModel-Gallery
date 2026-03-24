##
## SEIR Model extension: Syphilis Progression Model
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##
## Authors: Samuel M. Jenness, Yuan Zhao
## Date: November 2018
##


# Replacement infection/transmission module -------------------------------

infect <- function(dat, at) {

  ## Attributes ##
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")


  ## Initialize all custom attributes at the first timestep ##
  ## EpiModel convention: modules first run at timestep 2 (timestep 1 is
  ## reserved for network initialization via init.net). So at == 2 is where
  ## we set up any custom attributes not handled by EpiModel's built-in init.
  if (at == 2) {
    # Syphilis stage: NA for susceptible, "incubating" for initially infected
    dat <- set_attr(dat, "syph.stage",
                    ifelse(status == "i", "incubating", NA))
    # Symptomatic indicator: 0 = asymptomatic (default for all)
    dat <- set_attr(dat, "syph.symp", rep(0, length(active)))
    # Infection time
    dat <- set_attr(dat, "infTime", ifelse(status == "i", 1, NA))
    # Stage transition timestamps
    dat <- set_attr(dat, "priTime", rep(NA, length(active)))
    dat <- set_attr(dat, "secTime", rep(NA, length(active)))
    dat <- set_attr(dat, "elTime", rep(NA, length(active)))
    dat <- set_attr(dat, "llTime", rep(NA, length(active)))
    dat <- set_attr(dat, "terTime", rep(NA, length(active)))
    # Duration tracking per stage
    dat <- set_attr(dat, "syph.dur", rep(NA, length(active)))
    dat <- set_attr(dat, "syph2.dur", rep(NA, length(active)))
    dat <- set_attr(dat, "syph3.dur", rep(NA, length(active)))
    dat <- set_attr(dat, "syph4.dur", rep(NA, length(active)))
    dat <- set_attr(dat, "syph5.dur", rep(NA, length(active)))
    dat <- set_attr(dat, "syph6.dur", rep(NA, length(active)))
    # Treatment and screening indicators
    dat <- set_attr(dat, "syph.trt", rep(NA, length(active)))
    dat <- set_attr(dat, "syph.scr", rep(0, length(active)))
    dat <- set_attr(dat, "trtTime", rep(NA, length(active)))
    dat <- set_attr(dat, "scrTime", rep(NA, length(active)))
  }

  syph.stage <- get_attr(dat, "syph.stage")
  infTime <- get_attr(dat, "infTime")

  ## Parameters ##
  inf.prob.early <- get_param(dat, "inf.prob.early")
  inf.prob.latent <- get_param(dat, "inf.prob.latent")
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
    if (!(is.null(del))) {

      # Transmission probability depends on the infector's syphilis stage:
      #   Incubating, primary, secondary: higher probability (inf.prob.early)
      #   Early latent: lower probability (inf.prob.latent)
      #   Late latent, tertiary: not infectious (probability = 0)
      del$transProb <- ifelse(
        syph.stage[del$inf] %in% c("incubating", "primary", "secondary"),
        inf.prob.early, inf.prob.latent)
      del$transProb <- ifelse(
        syph.stage[del$inf] %in% c("late_latent", "tertiary"),
        0, del$transProb)

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
        syph.stage[idsNewInf] <- "incubating"
        status[idsNewInf] <- "i"
        infTime[idsNewInf] <- at
      }
    }
  }


  ## Save updated attributes ##
  dat <- set_attr(dat, "status", status)
  dat <- set_attr(dat, "syph.stage", syph.stage)
  dat <- set_attr(dat, "infTime", infTime)

  ## Save summary statistics ##
  dat <- set_epi(dat, "si.flow", at, nInf)
  syph.dur <- get_attr(dat, "syph.dur")
  syph.dur <- ifelse(syph.stage == "incubating", (at - infTime), syph.dur)
  dat <- set_attr(dat, "syph.dur", syph.dur)
  dat <- set_epi(dat, "syph.dur", at, mean(syph.dur, na.rm = TRUE))

  return(dat)
}


# Syphilis progression module ------------------------------------------
# (Replaces the recovery module)

progress <- function(dat, at) {

  ## Attributes ##
  active <- get_attr(dat, "active")
  syph.stage <- get_attr(dat, "syph.stage")
  syph.symp <- get_attr(dat, "syph.symp")
  infTime <- get_attr(dat, "infTime")
  priTime <- get_attr(dat, "priTime")
  secTime <- get_attr(dat, "secTime")
  elTime <- get_attr(dat, "elTime")
  llTime <- get_attr(dat, "llTime")
  terTime <- get_attr(dat, "terTime")

  ## Parameters of stage transition ##
  ipr.rate <- get_param(dat, "ipr.rate")
  prse.rate <- get_param(dat, "prse.rate")
  seel.rate <- get_param(dat, "seel.rate")
  elll.rate <- get_param(dat, "elll.rate")
  llter.rate <- get_param(dat, "llter.rate")

  ## Parameters of symptomatic progression ##
  pri.sym <- get_param(dat, "pri.sym")
  sec.sym <- get_param(dat, "sec.sym")


  ## Incubation to primary stage progression process ##
  nPrim <- 0
  idsEligPri <- which(active == 1 & syph.stage == "incubating" & infTime < at)
  nEligPri <- length(idsEligPri)

  if (nEligPri > 0) {
    vecPri <- which(rbinom(nEligPri, 1, ipr.rate) == 1)
    if (length(vecPri) > 0) {
      idsPri <- idsEligPri[vecPri]
      nPrim  <- length(idsPri)
      syph.stage[idsPri] <- "primary"
      priTime[idsPri] <- at
      ## Primary stage symptomatic progression ##
      syph.symp[idsPri] <- sample(0:1, size = length(vecPri),
                                  prob = c(1 - pri.sym, pri.sym), replace = TRUE)
    }
  }

  syph2.dur <- get_attr(dat, "syph2.dur")
  dat <- set_attr(dat, "syph2.dur",
                  ifelse(syph.stage == "primary", (at - priTime), syph2.dur))


  ## Primary to secondary stage progression ##
  nSec <- 0
  idsEligSec <- which(active == 1 & syph.stage == "primary" & priTime < at)
  nEligSec <- length(idsEligSec)

  if (nEligSec > 0) {
    vecSec <- which(rbinom(nEligSec, 1, prse.rate) == 1)
    if (length(vecSec) > 0) {
      idsSec <- idsEligSec[vecSec]
      nSec <- length(idsSec)
      syph.stage[idsSec] <- "secondary"
      syph.symp[idsSec] <- 0
      secTime[idsSec] <- at
      ## Secondary stage symptomatic progression ##
      syph.symp[idsSec] <- sample(0:1, size = length(vecSec),
                                  prob = c(1 - sec.sym, sec.sym), replace = TRUE)
    }
  }

  syph3.dur <- get_attr(dat, "syph3.dur")
  dat <- set_attr(dat, "syph3.dur",
                  ifelse(syph.stage == "secondary", (at - secTime), syph3.dur))


  ## Secondary to early latent progression ##
  nEarL <- 0
  idsEligLat <- which(active == 1 & syph.stage == "secondary" & secTime < at)
  nEligLat <- length(idsEligLat)

  if (nEligLat > 0) {
    vecLat <- which(rbinom(nEligLat, 1, seel.rate) == 1)
    if (length(vecLat) > 0) {
      idsLat <- idsEligLat[vecLat]
      nEarL <- length(idsLat)
      syph.stage[idsLat] <- "early_latent"
      syph.symp[idsLat] <- 0
      elTime[idsLat] <- at
    }
  }

  syph4.dur <- get_attr(dat, "syph4.dur")
  dat <- set_attr(dat, "syph4.dur",
                  ifelse(syph.stage == "early_latent", (at - elTime), syph4.dur))


  ## Early latent to late latent progression ##
  nLaL <- 0
  idsEligLaL <- which(active == 1 & syph.stage == "early_latent" & elTime < at)
  nEligLaL <- length(idsEligLaL)

  if (nEligLaL > 0) {
    vecLal <- which(rbinom(nEligLaL, 1, elll.rate) == 1)
    if (length(vecLal) > 0) {
      idsLal <- idsEligLaL[vecLal]
      nLaL <- length(idsLal)
      syph.stage[idsLal] <- "late_latent"
      syph.symp[idsLal] <- 0
      llTime[idsLal] <- at
    }
  }

  syph5.dur <- get_attr(dat, "syph5.dur")
  llTime <- get_attr(dat, "llTime")
  dat <- set_attr(dat, "syph5.dur",
                  ifelse(syph.stage == "late_latent", (at - llTime), syph5.dur))


  ## Late latent to tertiary progression ##
  nTer <- 0
  idsEligTer <- which(active == 1 & syph.stage == "late_latent" & llTime < at)
  nEligTer <- length(idsEligTer)

  if (nEligTer > 0) {
    vecTer <- which(rbinom(nEligTer, 1, llter.rate) == 1)
    if (length(vecTer) > 0) {
      idsTer <- idsEligTer[vecTer]
      nTer <- length(idsTer)
      syph.stage[idsTer] <- "tertiary"
      syph.symp[idsTer] <- 1
      terTime[idsTer] <- at
    }
  }

  syph6.dur <- get_attr(dat, "syph6.dur")
  dat <- set_attr(dat, "syph6.dur",
                  ifelse(syph.stage == "tertiary", (at - terTime), syph6.dur))


  ## Save updated attributes ##
  dat <- set_attr(dat, "syph.stage", syph.stage)
  dat <- set_attr(dat, "syph.symp", syph.symp)
  dat <- set_attr(dat, "priTime", priTime)
  dat <- set_attr(dat, "secTime", secTime)
  dat <- set_attr(dat, "elTime", elTime)
  dat <- set_attr(dat, "llTime", llTime)
  dat <- set_attr(dat, "terTime", terTime)


  ## Save summary statistics ##
  dat <- set_epi(dat, "ipr.flow", at, nPrim)
  dat <- set_epi(dat, "prse.flow", at, nSec)
  dat <- set_epi(dat, "seel.flow", at, nEarL)
  dat <- set_epi(dat, "elll.flow", at, nLaL)
  dat <- set_epi(dat, "llter.flow", at, nTer)

  dat <- set_epi(dat, "inc.num", at,
                 sum(active == 1 & syph.stage == "incubating", na.rm = TRUE))
  dat <- set_epi(dat, "pr.num", at,
                 sum(active == 1 & syph.stage == "primary", na.rm = TRUE))
  dat <- set_epi(dat, "se.num", at,
                 sum(active == 1 & syph.stage == "secondary", na.rm = TRUE))
  dat <- set_epi(dat, "el.num", at,
                 sum(active == 1 & syph.stage == "early_latent", na.rm = TRUE))
  dat <- set_epi(dat, "ll.num", at,
                 sum(active == 1 & syph.stage == "late_latent", na.rm = TRUE))
  dat <- set_epi(dat, "ter.num", at,
                 sum(active == 1 & syph.stage == "tertiary", na.rm = TRUE))
  dat <- set_epi(dat, "sym.num", at,
                 sum(active == 1 & syph.symp == 1, na.rm = TRUE))

  dat <- set_epi(dat, "syph2.dur", at, mean(syph2.dur, na.rm = TRUE))
  dat <- set_epi(dat, "syph3.dur", at, mean(syph3.dur, na.rm = TRUE))
  dat <- set_epi(dat, "syph4.dur", at, mean(syph4.dur, na.rm = TRUE))
  dat <- set_epi(dat, "syph5.dur", at, mean(syph5.dur, na.rm = TRUE))
  dat <- set_epi(dat, "syph6.dur", at, mean(syph6.dur, na.rm = TRUE))

  return(dat)
}


# Treatment and screening module ------------------------------------------

tnt <- function(dat, at) {

  ## Attributes ##
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  syph.stage <- get_attr(dat, "syph.stage")
  syph.symp <- get_attr(dat, "syph.symp")
  syph.trt <- get_attr(dat, "syph.trt")
  syph.scr <- get_attr(dat, "syph.scr")
  trtTime <- get_attr(dat, "trtTime")
  scrTime <- get_attr(dat, "scrTime")

  ## Parameters ##
  early.trt <- get_param(dat, "early.trt")
  late.trt <- get_param(dat, "late.trt")
  scr.rate <- get_param(dat, "scr.rate")


  # --- Symptomatic treatment initiation ---

  ## Primary symptomatic patients receiving treatment ##
  priTime <- get_attr(dat, "priTime")
  idsPriTrt <- which(active == 1 & syph.stage == "primary" & syph.symp == 1 &
                       is.na(syph.trt) & priTime < at)
  if (length(idsPriTrt) > 0) {
    vecPriTrt <- which(rbinom(length(idsPriTrt), 1, early.trt) == 1)
    if (length(vecPriTrt) > 0) {
      idsPriTrt <- idsPriTrt[vecPriTrt]
      syph.trt[idsPriTrt] <- 1
      trtTime[idsPriTrt] <- at
    }
  }

  ## Secondary symptomatic patients receiving treatment ##
  secTime <- get_attr(dat, "secTime")
  idsSecTrt <- which(active == 1 & syph.stage == "secondary" & syph.symp == 1 &
                       is.na(syph.trt) & secTime < at)
  if (length(idsSecTrt) > 0) {
    vecSecTrt <- which(rbinom(length(idsSecTrt), 1, early.trt) == 1)
    if (length(vecSecTrt) > 0) {
      idsSecTrt <- idsSecTrt[vecSecTrt]
      syph.trt[idsSecTrt] <- 1
      trtTime[idsSecTrt] <- at
    }
  }

  ## Tertiary symptomatic patients receiving treatment ##
  terTime <- get_attr(dat, "terTime")
  idsTerTrt <- which(active == 1 & syph.stage == "tertiary" & syph.symp == 1 &
                       is.na(syph.trt) & terTime < at)
  if (length(idsTerTrt) > 0) {
    vecTerTrt <- which(rbinom(length(idsTerTrt), 1, late.trt) == 1)
    if (length(vecTerTrt) > 0) {
      idsTerTrt <- idsTerTrt[vecTerTrt]
      syph.trt[idsTerTrt] <- 1
      trtTime[idsTerTrt] <- at
    }
  }


  # --- Symptomatic treatment recovery ---

  ## Recover 1 week after treatment (primary) ##
  idsPriRec <- which(active == 1 & syph.stage == "primary" & syph.trt == 1 &
                       trtTime <= at - 1)
  if (length(idsPriRec) > 0) {
    status[idsPriRec] <- "s"
    syph.stage[idsPriRec] <- NA
    syph.trt[idsPriRec] <- NA
    syph.symp[idsPriRec] <- 0
  }

  ## Recover 1 week after treatment (secondary) ##
  idsSecRec <- which(active == 1 & syph.stage == "secondary" & syph.trt == 1 &
                       trtTime <= at - 1)
  if (length(idsSecRec) > 0) {
    status[idsSecRec] <- "s"
    syph.stage[idsSecRec] <- NA
    syph.trt[idsSecRec] <- NA
    syph.symp[idsSecRec] <- 0
  }

  ## Recover 3 weeks after treatment (tertiary) ##
  idsTerRec <- which(active == 1 & syph.stage == "tertiary" & syph.trt == 1 &
                       trtTime <= at - 3)
  if (length(idsTerRec) > 0) {
    status[idsTerRec] <- "s"
    syph.stage[idsTerRec] <- NA
    syph.trt[idsTerRec] <- NA
    syph.symp[idsTerRec] <- 0
  }


  # --- Screening of asymptomatic infected population ---

  nScr <- 0
  idsEligScr <- which(active == 1 & status == "i" & syph.symp == 0 &
                        is.na(syph.trt))
  nEligScr <- length(idsEligScr)

  if (nEligScr > 0) {
    vecScr <- which(rbinom(nEligScr, 1, scr.rate) == 1)
    if (length(vecScr) > 0) {
      idsScr <- idsEligScr[vecScr]
      nScr <- length(idsScr)
      syph.scr[idsScr] <- 1
      syph.trt[idsScr] <- 1
      trtTime[idsScr] <- at
      scrTime[idsScr] <- at
    }
  }


  # --- Screening-detected recovery ---

  ## Recovery for screening-detected patients in non-tertiary stages.
  ## Primary/secondary symptomatic cases are typically caught by the blocks
  ## above; this primarily handles incubating, early latent, and late latent
  ## patients identified through screening. Tertiary is handled above.
  idsRecScr <- which(active == 1 &
                       syph.stage %in% c("incubating", "primary", "secondary",
                                         "early_latent", "late_latent") &
                       syph.trt == 1 & trtTime < at - 1)
  if (length(idsRecScr) > 0) {
    status[idsRecScr] <- "s"
    syph.stage[idsRecScr] <- NA
    syph.trt[idsRecScr] <- NA
    syph.symp[idsRecScr] <- 0
  }


  ## Count total recoveries this timestep ##
  nRec <- length(idsPriRec) + length(idsSecRec) + length(idsTerRec) +
    length(idsRecScr)


  ## Save all modified attributes ##
  dat <- set_attr(dat, "status", status)
  dat <- set_attr(dat, "syph.stage", syph.stage)
  dat <- set_attr(dat, "syph.symp", syph.symp)
  dat <- set_attr(dat, "syph.trt", syph.trt)
  dat <- set_attr(dat, "syph.scr", syph.scr)
  dat <- set_attr(dat, "trtTime", trtTime)
  dat <- set_attr(dat, "scrTime", scrTime)

  ## Save summary statistics ##
  dat <- set_epi(dat, "rec.flow", at, nRec)
  dat <- set_epi(dat, "scr.flow", at, nScr)
  dat <- set_epi(dat, "scr.num", at,
                 sum(active == 1 & syph.scr == 1, na.rm = TRUE))
  dat <- set_epi(dat, "trt.num", at,
                 sum(active == 1 & syph.trt == 1, na.rm = TRUE))

  return(dat)
}
