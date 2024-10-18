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


  ## Initiating an indicator of syphilis stage##
  if (at == 2) {
    status <- get_attr(dat, "status")
    dat  <- set_attr(dat, "syph.stage",
                     ifelse(status == "i", 1, 0))
    dat <- set_attr(dat, "syph.symp",
                    rep(0, length(active)))

    status <- get_attr(dat, "status")
    dat <- set_attr(dat, "infTime", ifelse(status == "i", 1, NA))
    dat <- set_attr(dat, "priTime", rep(NA, length(active)))
    dat <- set_attr(dat, "secTime", rep(NA, length(active)))
    dat <- set_attr(dat, "elTime", rep(NA, length(active)))
    dat <- set_attr(dat, "llTime", rep(NA, length(active)))
    dat <- set_attr(dat, "terTime", rep(NA, length(active)))
    dat <- set_attr(dat, "syph.dur", rep(NA, length(active)))
    dat <- set_attr(dat, "syph2.dur", rep(NA, length(active)))
    dat <- set_attr(dat, "syph3.dur", rep(NA, length(active)))
    dat <- set_attr(dat, "syph4.dur", rep(NA, length(active)))
    dat <- set_attr(dat, "syph5.dur", rep(NA, length(active)))
    dat <- set_attr(dat, "syph6.dur", rep(NA, length(active)))
  }

  syph.stage <- get_attr(dat, "syph.stage")
  infTime <- get_attr(dat, "infTime")

  ## Parameters ##
  inf.prob1 <- get_param(dat, "inf.prob1")
  inf.prob2 <- get_param(dat, "inf.prob2")
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

      # Set parameters on discordant edgelist data frame
      del$transProb <- ifelse(syph.stage[del$inf] < 4, inf.prob1, inf.prob2)
      del$transProb <- ifelse(syph.stage[del$inf] > 4, 0, del$transProb)

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
        syph.stage[idsNewInf] <- 1
        status[idsNewInf] <- "i"
        dat <- set_attr(dat, "status", status)
        infTime[idsNewInf] <- at
        dat <- set_attr(dat, "infTime", infTime)
      }
    }
  }


  ## Save summary statistic for s->i flow
  dat <- set_epi(dat, "si.flow", at, nInf)
  dat <- set_attr(dat, "syph.stage", syph.stage)
  infTime <- get_attr(dat, "infTime")
  syph.dur <- get_attr(dat, "syph.dur")
  syph.dur <- ifelse(syph.stage == 1, (at - infTime), syph.dur)
  dat <- set_attr(dat, "syph.dur", syph.dur)
  dat <- set_epi(dat, "syph.dur", at, mean(syph.dur, na.rm = TRUE))

  return(dat)

}


# Syphilis progression module ------------------------------------------
# (Replaces the recovery module)

progress <- function(dat, at) {

  ## Attributes: syphilis stage, symptomatic or not, treatment and
  ## screening inidicators ##
  active <- get_attr(dat, "active")
  syph.stage <- get_attr(dat, "syph.stage")
  syph.symp <- get_attr(dat, "syph.symp")
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

  ## Parameters of symptomatic ##
  pri.sym <- get_param(dat, "pri.sym")
  sec.sym <- get_param(dat, "sec.sym")

  ## Incubation to primary stage progression process ##
  nPrim <- 0
  idsEligPri <- which(active == 1 & syph.stage == 1)
  nEligPri <- length(idsEligPri)

  if (nEligPri > 0) {
    vecPri <- which(rbinom(nEligPri, 1, ipr.rate) == 1)
    if (length(vecPri) > 0) {
      idsPri <- idsEligPri[vecPri]
      nPrim  <- length(idsPri)
      syph.stage[idsPri] <- 2
      dat <- set_attr(dat, "priTime", at, posit_ids = idsPri)
      ## Primary Stage Symptomatic Progression ##
      syph.symp[idsPri] <- sample(0:1, size = length(vecPri),
                                  prob = c(1 - pri.sym, pri.sym), replace = TRUE)
    }
  }

  syph2.dur <- get_attr(dat, "syph2.dur")
  dat <- set_attr(dat, "syph2.dur",
                  ifelse(syph.stage == 2, (at - priTime), syph2.dur))

  ## Primary to secondary stage progression ##
  nSec <- 0
  idsEligSec <- which(active == 1 & syph.stage == 2 & priTime < at)
  nEligSec <- length(idsEligSec)

  if (nEligSec > 0) {
    vecSec <- which(rbinom(nEligSec, 1, prse.rate) == 1)
    if (length(vecSec) > 0) {
      idsSec <- idsEligSec[vecSec]
      nSec <- length(idsSec)
      syph.stage[idsSec] <- 3
      syph.symp[idsSec] <- 0
      secTime[idsSec] <- at
      dat <- set_attr(dat, "secTime", secTime)
      ## Secondary Stage Symptomatic Progression ##
      syph.symp[idsSec] <- sample(0:1, size = length(vecSec),
                                  prob = c(1 - sec.sym, sec.sym), replace = TRUE)
    }
  }

  syph3.dur <- get_attr(dat, "syph3.dur")
  dat <- set_attr(dat, "syph3.dur",
                  ifelse(syph.stage == 3, (at - secTime), syph3.dur))

  ## Secondary to early latent progression ##
  nEarL <- 0
  idsEligLat <- which(active == 1 & syph.stage == 3 & secTime < at)
  nEligLat <- length(idsEligLat)

  if (nEligLat > 0) {
    vecLat <- which(rbinom(nEligLat, 1, seel.rate) == 1)
    if (length(vecLat) > 0) {
      idsLat <- idsEligLat[vecLat]
      nEarL <- length(idsLat)
      syph.stage[idsLat] <- 4
      syph.symp[idsLat] <- 0
      elTime[idsLat] <- at
      dat <- set_attr(dat, "elTime", elTime)
    }
  }

  syph4.dur <- get_attr(dat, "syph4.dur")
  dat <- set_attr(dat, "syph4.dur",
                  ifelse(syph.stage == 4, (at - elTime), syph4.dur))

  ## Early latent to late latent progression ##
  nLaL <- 0
  idsEligLaL <- which(active == 1 & syph.stage == 4 & elTime < at)
  nEligLaL <- length(idsEligLaL)
  if (nEligLaL > 0) {
    vecLal <- which(rbinom(nEligLaL, 1, elll.rate) == 1)
    if (length(vecLal) > 0) {
      idsLal <- idsEligLaL[vecLal]
      nLaL <- length(idsLal)
      syph.stage[idsLal] <- 5
      syph.symp[idsLal] <- 0
      llTime[idsLal] <- at
      dat <- set_attr(dat, "llTime", llTime)
    }
  }

  syph5.dur <- get_attr(dat, "syph5.dur")
  llTime <- get_attr(dat, "llTime")
  dat <- set_attr(dat, "syph5.dur",
                  ifelse(syph.stage == 5, (at - llTime), syph5.dur))

  ## Late latent to Tertiary progression ##
  nTer <- 0
  idsEligTer <- which(active == 1 & syph.stage == 5 & llTime < at)
  nEligTer <- length(idsEligTer)

  if (nEligTer > 0) {
    vecTer <- which(rbinom(nEligTer, 1, llter.rate) == 1)
    if (length(vecTer) > 0) {
      idsTer <- idsEligTer[vecTer]
      nTer <- length(idsTer)
      syph.stage[idsTer] <- 6
      syph.symp[idsTer] <- 1
      terTime[idsTer] <- at
      dat <- set_attr(dat, "terTime", terTime)
    }
  }

  syph6.dur <- get_attr(dat, "syph6.dur")
  dat <- set_attr(dat, "syph6.dur",
                  ifelse(syph.stage == 6, (at - terTime), syph6.dur))

  dat <- set_attr(dat, "syph.stage", syph.stage)
  dat <- set_attr(dat, "syph.symp", syph.symp)


  ## Save summary statistics ##
  dat <- set_epi(dat, "ipr.flow", at, nPrim)
  dat <- set_epi(dat, "prse.flow", at, nSec)
  dat <- set_epi(dat, "seel.flow", at, nEarL)
  dat <- set_epi(dat, "elll.flow", at, nLaL)
  dat <- set_epi(dat, "llter.flow", at, nTer)

  dat <- set_epi(dat, "inc.num", at, sum(active == 1 & syph.stage == 1))
  dat <- set_epi(dat, "pr.num", at, sum(active == 1 & syph.stage == 2))
  dat <- set_epi(dat, "se.num", at, sum(active == 1 & syph.stage == 3))
  dat <- set_epi(dat, "el.num", at, sum(active == 1 & syph.stage == 4))
  dat <- set_epi(dat, "ll.num", at, sum(active == 1 & syph.stage == 5))
  dat <- set_epi(dat, "ter.num", at, sum(active == 1 & syph.stage == 6))
  dat <- set_epi(dat, "sym.num", at, sum(active == 1 & syph.symp == 1))

  dat <- set_epi(dat, "syph2.dur", at, mean(syph2.dur, na.rm = TRUE))
  dat <- set_epi(dat, "syph3.dur", at, mean(syph3.dur, na.rm = TRUE))
  dat <- set_epi(dat, "syph4.dur", at, mean(syph4.dur, na.rm = TRUE))
  dat <- set_epi(dat, "syph5.dur", at, mean(syph5.dur, na.rm = TRUE))
  dat <- set_epi(dat, "syph6.dur", at, mean(syph6.dur, na.rm = TRUE))

  return(dat)
}

# Treatment and testing module ------------------------------------------
tnt <- function(dat, at) {

  ## Attributes: syphilis stage, symptomatic or not, treatment and
  ## screening indicators ##
  active <- get_attr(dat, "active")
  syph.stage <- get_attr(dat, "syph.stage")
  syph.symp <- get_attr(dat, "syph.symp")
  status <- get_attr(dat, "status")

  if (at == 2) {
    dat <- set_attr(dat, "syph.trt", rep(NA, length(active)))
    dat <- set_attr(dat, "syph.scr", rep(0, length(active)))
    dat <- set_attr(dat, "trtTime", rep(NA, length(active)))
    dat <- set_attr(dat, "scrTime", rep(NA, length(active)))
  }

  syph.trt <- get_attr(dat, "syph.trt")
  syph.scr <- get_attr(dat, "syph.scr")
  priTime <- get_attr(dat, "priTime")
  trtTime <- get_attr(dat, "trtTime")
  scrTime <- get_attr(dat, "scrTime")

  ## Parameters of treatment and screening ##
  early.trt <- get_param(dat, "early.trt")
  late.trt <- get_param(dat, "late.trt")
  scr.rate <- get_param(dat, "scr.rate")


  ## Primary symptomatic patients receiving treatment
  idsPriTrt <- which(active == 1 & syph.stage == 2 & syph.symp == 1 &
                       is.na(syph.trt) & priTime < at)
  nPriTrt <- length(idsPriTrt)

  if (nPriTrt > 0) {
    vecPriTrt <- which(rbinom(nPriTrt, 1, early.trt) == 1)
    if (length(vecPriTrt) > 0) {
      idsPriTrt <- idsPriTrt[vecPriTrt]
      syph.trt[idsPriTrt] <- 1
      trtTime[idsPriTrt] <- at
      dat <- set_attr(dat, "trtTime", trtTime)
    }
  }

  ## Recover in 1 week after treatment ##
  trtTime <- get_attr(dat, "trtTime")
  idsRec <- which(active == 1 & syph.stage == 2 & syph.trt == 1 &
                    trtTime <= at - 1)
  status[idsRec] <- "s"
  dat <- set_attr(dat, "status", status)
  syph.stage[idsRec] <- 0
  syph.trt[idsRec] <- NA
  syph.symp[idsRec] <- NA


  ## Secondary symptomatic patients receiving treatment ##
  secTime <- get_attr(dat, "secTime")
  idsSecTrt <- which(active == 1 & syph.stage == 3 & syph.symp == 1 &
                       is.na(syph.trt) & secTime < at)
  nSecTrt <- length(idsSecTrt)

  if (nSecTrt > 0) {
    vecSecTrt <- which(rbinom(nSecTrt, 1, early.trt) == 1)
    if (length(vecSecTrt) > 0) {
      idsSecTrt <- idsSecTrt[vecSecTrt]
      syph.trt[idsSecTrt] <- 1
      trtTime[idsSecTrt] <- at
      dat <- set_attr(dat, "trtTime", trtTime)
    }
  }

  ## Recover in 1 week after treatment ##
  trtTime <- get_attr(dat, "trtTime")
  idsSecRec <- which(active == 1 & syph.stage == 3 & syph.trt == 1 &
                       trtTime <= at - 1)
  status[idsSecRec] <- "s"
  syph.stage[idsSecRec] <- 0
  syph.trt[idsSecRec] <- NA
  syph.symp[idsSecRec] <- NA


  ## Tertiary Symptomatic Patients are put on treatment ##
  terTime <- get_attr(dat, "terTime")
  idsTerTrt <- which(active == 1 & syph.stage == 6 & syph.symp == 1 &
                       is.na(syph.trt) & terTime < at)
  nTerTrt <- length(idsTerTrt)

  if (nTerTrt > 0) {
    vecTerTrt <- which(rbinom(nTerTrt, 1, late.trt) == 1)
    if (length(vecTerTrt) > 0) {
      idsTerTrt <- idsTerTrt[vecTerTrt]
      syph.trt[idsTerTrt] <- 1
      trtTime[idsTerTrt] <- at
      dat <- set_attr(dat, "trtTime", trtTime)
    }
  }

  ## Recovery after treatment in tertiary stage
  idsTerRec <- which(active == 1 & syph.stage == 6 & syph.trt == 1 &
                       trtTime <= at - 3)
  status[idsTerRec] <- "s"
  syph.stage[idsTerRec] <- 0
  syph.symp[idsTerRec] <- NA
  syph.trt[idsTerRec] <- NA


  ## Screening of asymptomatic population who are not on treatment##
  nScr <- 0
  idsEligScr <- which(active == 1 & (syph.symp == 0 | is.na(syph.symp)) &
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
      dat <- set_attr(dat, "trtTime", trtTime)
      scrTime[idsScr] <- at
      dat <- set_attr(dat, "scrTime", scrTime)
    }
  }


  idsRec1 <- which(active == 1 & syph.stage < 4 & syph.trt == 1 &
                     trtTime < at - 1)
  status[idsRec1] <- "s"
  syph.stage[idsRec1] <- 0
  syph.symp[idsRec1] <- 0

  idsRec2 <- which(active == 1 & syph.stage == 6 & syph.trt == 1 &
                     trtTime <= at - 3)
  status[idsRec2] <- "s"
  syph.stage[idsRec2] <- 0
  syph.symp[idsRec2] <- 0

  dat <- set_epi(dat, "scr.flow", at, nScr)
  dat <- set_epi(dat, "scr.num", at, sum(active == 1 & syph.scr == 1))
  dat <- set_epi(dat, "trt.num", at, sum(active == 1 & syph.trt == 1))

  dat <- set_attr(dat, "syph.stage", syph.stage)
  dat <- set_attr(dat, "syph.symp", syph.symp)

  return(dat)
}
