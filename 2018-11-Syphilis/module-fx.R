##
## SEIR Model extension: Syphilis Progression Model
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##
## Authors: Samuel M. Jenness, Yuan Zhao
## Date: November 2018
##

## Natural history of syphilis with testing and treatment
# Replacement infection/transmission module -------------------------------

infect <- function(dat, at) {

  ## Uncomment this to function environment interactively
  ## browser()
  ## Attributes ##
  active <- dat$attr$active
  status <- dat$attr$status
  
  ## Initiating an indicator of syphilis status to indicate stage of syphilis##
  if (at == 2) {
    dat$attr$syph.stage <- rep(0, length(active))
    dat$attr$syph.stage <- ifelse(dat$attr$status == "i",1,0)
    dat$attr$syph.symp <- rep(0, length(active))
    
    dat$attr$infTime <- rep(NA, length(active))
    dat$attr$infTime <- ifelse(dat$attr$status == "i",1,dat$attr$infTime)
    dat$attr$priTime <- rep(NA, length(active))
    dat$attr$secTime <- rep(NA, length(active))
    dat$attr$elTime <- rep(NA, length(active))
    dat$attr$llTime <- rep(NA, length(active))
    dat$attr$terTime <- rep(NA, length(active))
    
    dat$attr$syph.dur <- rep(NA, length(active))
    dat$attr$syph2.dur <-  rep(NA, length(active))
    dat$attr$syph3.dur <-  rep(NA, length(active))
    dat$attr$syph4.dur <-  rep(NA, length(active))
    dat$attr$syph5.dur <-  rep(NA, length(active))
    dat$attr$syph6.dur <-  rep(NA, length(active))
    
  }
  
  syph.stage <- dat$attr$syph.stage
  
  ## Parameters ##
  inf.prob1 <- dat$param$inf.prob1
  inf.prob2 <- dat$param$inf.prob2
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
      del$transProb <- ifelse(syph.stage[del$inf] < 4,inf.prob1,inf.prob2)
      del$transProb <- ifelse(syph.stage[del$inf] > 4,0,del$transProb)
      
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
        dat$attr$status[idsNewInf] <- "i"
        dat$attr$infTime[idsNewInf] <- at
        
      }
      
      
    }
    
  }
  ## Save summary statistic for S->i flow
  dat$epi$si.flow[at] <- nInf
  dat$attr$syph.stage <- syph.stage
  dat$attr$syph.dur <- ifelse(syph.stage == 1,(at - dat$attr$infTime),
                              dat$attr$syph.dur)
  dat$epi$syph.dur[at] <- mean(dat$attr$syph.dur,na.rm = TRUE)
  return(dat)
}


# Syphilis progression module ------------------------------------------
# (Replaces the recovery module)

progress <- function(dat, at) {

  ## Uncomment this to function environment interactively
  ## browser()

  ## Attributes: syphilis stage, symptomatic or not, treatment and 
  ## screening inidicators ##
  active <- dat$attr$active
  syph.stage <- dat$attr$syph.stage
  syph.symp <- dat$attr$syph.symp

  ## Parameters of stage transition ##
  ipr.rate <- dat$param$ipr.rate
  prse.rate <- dat$param$prse.rate
  seel.rate <- dat$param$seel.rate
  elll.rate <- dat$param$elll.rate
  llter.rate <- dat$param$llter.rate
  
  ## Parameters of symptomatic ##
  pri.sym <- dat$param$pri.sym 
  sec.sym <- dat$param$sec.sym
  
  
  ## Incubation to primary stage progression process ##
  ## browser()
  nPrim <- 0
  idsEligPri <- which(active == 1 & syph.stage == 1)
  nEligPri <- length(idsEligPri)

  if (nEligPri > 0) {
    vecPri <- which(rbinom(nEligPri, 1, ipr.rate) == 1)
    if (length(vecPri) > 0) {
      idsPri <- idsEligPri[vecPri]
      nPrim  <- length(idsPri)
      syph.stage[idsPri] <- 2
      dat$attr$priTime[idsPri] <- at
      ## Primary Stage Symptomatic Progression ##
      syph.symp[idsPri] <- sample(0:1, size = length(vecPri),
                        prob = c(1 - pri.sym, pri.sym),replace = TRUE)
    }
  }
  dat$attr$syph2.dur <- ifelse(syph.stage == 2,
                               (at - dat$attr$priTime),dat$attr$syph2.dur)
  
 
  ## Primary to secondary stage progression ##
  nSec <- 0
  idsEligSec <- which(active == 1 & syph.stage == 2 & dat$attr$priTime < at)
  nEligSec <- length(idsEligSec)

  if (nEligSec > 0) {
    vecSec <- which(rbinom(nEligSec, 1, prse.rate) == 1)
    if (length(vecSec) > 0) {
      idsSec <- idsEligSec[vecSec]
      nSec <- length(idsSec)
      syph.stage[idsSec] <- 3
      syph.symp[idsSec] <- 0
      dat$attr$secTime[idsSec] <- at
      ## Secondary Stage Symptomatic Progression ##
      syph.symp[idsSec] <- sample(0:1, size = length(vecSec),
                          prob = c(1 - sec.sym,sec.sym),replace = TRUE)
    }
  }
  dat$attr$syph3.dur <- ifelse(syph.stage == 3,
                               (at - dat$attr$secTime),dat$attr$syph3.dur)
  
  
  ## Secondary to early latent progression ##
  nEarL <- 0
  idsEligLat <- which(active == 1 & syph.stage == 3 & dat$attr$secTime < at)
  nEligLat <- length(idsEligLat)
  
  if (nEligLat > 0) {
    vecLat <- which(rbinom(nEligLat, 1, seel.rate) == 1)
    if (length(vecLat) > 0) {
      idsLat <- idsEligLat[vecLat]
      nEarL <- length(idsLat)
      syph.stage[idsLat] <- 4
      syph.symp[idsLat] <- 0
      dat$attr$elTime[idsLat] <- at
    }
  }
  dat$attr$syph4.dur <- ifelse(syph.stage == 4,
                               (at - dat$attr$elTime),dat$attr$syph4.dur)
  
  ## Early latent to late latent progression ##
  nLaL <- 0
  idsEligLaL <- which(active == 1 & syph.stage == 4 & dat$attr$elTime < at)
  nEligLaL <- length(idsEligLaL)
  if (nEligLaL > 0) {
    vecLal <- which(rbinom(nEligLaL, 1, elll.rate) == 1)
    if (length(vecLal) > 0) {
      idsLal <- idsEligLaL[vecLal]
      nLaL <- length(idsLal)
      syph.stage[idsLal] <- 5
      syph.symp[idsLal] <- 0
      dat$attr$llTime[idsLal] <- at
    }
  }
  dat$attr$syph5.dur <- ifelse(syph.stage == 5,
                               (at - dat$attr$llTime),dat$attr$syph5.dur)
  
  ## Late latent to Tertiary progression ##
  nTer <- 0
  idsEligTer <- which(active == 1 & syph.stage == 5 & dat$attr$llTime < at)
  nEligTer <- length(idsEligTer)
  
  if (nEligTer > 0) {
    vecTer <- which(rbinom(nEligTer, 1, llter.rate) == 1)
    if (length(vecTer) > 0) {
      idsTer <- idsEligTer[vecTer]
      nTer <- length(idsTer)
      syph.stage[idsTer] <- 6
      syph.symp[idsTer] <- 1
      dat$attr$terTime[idsTer] <- at
    }
  }
  dat$attr$syph6.dur <- ifelse(syph.stage == 6,
                               (at - dat$attr$terTime),dat$attr$syph6.dur)
  
  dat$attr$syph.stage <- syph.stage
  dat$attr$syph.symp <- syph.symp
  
  ## Save summary statistics ##
  dat$epi$ipr.flow[at] <- nPrim
  dat$epi$prse.flow[at] <- nSec
  dat$epi$seel.flow[at] <- nEarL
  dat$epi$elll.flow[at] <- nLaL
  dat$epi$llter.flow[at] <- nTer
  
 
  dat$epi$inc.num[at] <- sum(active == 1 & syph.stage == 1)
  dat$epi$pr.num[at] <- sum(active == 1 & syph.stage == 2)
  dat$epi$se.num[at] <- sum(active == 1 & syph.stage == 3)
  dat$epi$el.num[at] <- sum(active == 1 & syph.stage == 4)
  dat$epi$ll.num[at] <- sum(active == 1 & syph.stage == 5)
  dat$epi$ter.num[at] <- sum(active == 1 & syph.stage == 6)
  dat$epi$sym.num[at] <- sum(active == 1 & syph.symp == 1)
  
  dat$epi$syph2.dur[at] <- mean(dat$attr$syph2.dur,na.rm = TRUE)
  dat$epi$syph3.dur[at] <- mean(dat$attr$syph3.dur,na.rm = TRUE)
  dat$epi$syph4.dur[at] <- mean(dat$attr$syph4.dur,na.rm = TRUE)
  dat$epi$syph5.dur[at] <- mean(dat$attr$syph5.dur,na.rm = TRUE)
  dat$epi$syph6.dur[at] <- mean(dat$attr$syph6.dur,na.rm = TRUE)
  
  return(dat)
}

# Treatment and testing module ------------------------------------------
tnt <- function(dat, at) {
  
  ## browser()
  
  ## Attributes: syphilis stage, symptomatic or not, treatment and 
  ## screening inidicators ##
  active <- dat$attr$active
  syph.stage <- dat$attr$syph.stage
  syph.symp <- dat$attr$syph.symp
  
  if (at == 2) {
    dat$attr$syph.trt <- rep(NA, length(active))
    dat$attr$syph.scr <- rep(0, length(active))
    dat$attr$trtTime <- rep(NA, length(active))
    dat$attr$scrTime <- rep(NA, length(active))
  }
  
  syph.trt <- dat$attr$syph.trt
  syph.scr <- dat$attr$syph.scr
  
  ## Parameters of treatment and screening ##
  early.trt <- dat$param$early.trt 
  scr.rate <- dat$param$scr.rate
  
  ## Primary symptomatic patients receiving treatment
  nPrim.trt <- 0
  idsPriTrt <- which(active == 1 & syph.stage == 2 & syph.symp == 1 & 
                        is.na(syph.trt) & dat$attr$priTime < at)
  nPriTrt <- length(idsPriTrt)
  
  if (nPriTrt > 0) {
    vecPriTrt <- which(rbinom(nPriTrt, 1, early.trt) == 1)
    if (length(vecPriTrt) > 0) {
      idsPriTrt <- idsPriTrt[vecPriTrt]
      nPrim.trt  <- length(idsPriTrt)
      syph.trt[idsPriTrt] <- 1
      dat$attr$trtTime[idsPriTrt] <- at
    }
  }

  ## Recover in 1 week after treatment ##
  idsRec <- which(active == 1 & syph.stage == 2 & syph.trt == 1 & 
                    dat$attr$trtTime < at - 1)
  dat$attr$status[idsRec] <- "s"
  syph.stage[idsRec] <- 0
  syph.trt[idsRec] <- NA
  syph.symp[idsRec] <- NA
  
  ## Secondary symptomatic patients receiving treatment ##
  nSec.trt <- 0
  idsSecTrt <- which(active == 1 & syph.stage == 3 & syph.symp == 1 & 
                        is.na(syph.trt) & dat$attr$secTime < at)
  nSecTrt <- length(idsSecTrt)
  
  if (nSecTrt > 0) {
    vecSecTrt <- which(rbinom(nSecTrt, 1, early.trt) == 1)
    if (length(vecSecTrt) > 0) {
      idsSecTrt <- idsSecTrt[vecSecTrt]
      nSec.trt  <- length(idsSecTrt)
      syph.trt[idsSecTrt] <- 1
      dat$attr$trtTime[idsSecTrt] <- at
    }
  }
  ## Recover in 1 week after treatment ##
  idsSecRec <- which(active == 1 & syph.stage == 3 & syph.trt == 1 & 
                    dat$attr$trtTime < at - 1)
  {
    dat$attr$status[idsSecRec] <- "s"
    syph.stage[idsSecRec] <- 0
    syph.trt[idsSecRec] <- NA
    syph.symp[idsSecRec] <- NA
  }
  
  ## Tertiary Symptomatic Patients are put on treatment ##
  idsTerTrt <- which(active == 1 & syph.stage == 6 & syph.symp == 1 & 
                        is.na(syph.trt) & dat$attr$terTime < at)
  nTer.trt <- length(idsTerTrt)
  syph.trt[idsTerTrt] <- 1
  dat$attr$trtTime[idsTerTrt] <- at
  idsTerRec <- which(active == 1 & syph.stage == 6 & syph.trt == 1 & 
                    dat$attr$trtTime < at - 3)
  {
    dat$attr$status[idsTerRec] <- "s"
    syph.stage[idsTerRec] <- 0
    syph.symp[idsTerRec] <- NA
    syph.trt[idsTerRec] <- NA
  }
  
  ## Screening of asymptomatic population who are not on treatment##
  nScr <- 0
  idsEligScr <- which(active == 1 & (syph.symp == 0|is.na(syph.symp)) & 
                        is.na(syph.trt))
  nEligScr <- length(idsEligScr)
  
  if (nEligScr > 0) {
    vecScr <- which(rbinom(nEligScr, 1, scr.rate) == 1)
    if (length(vecScr) > 0) {
      idsScr <- idsEligScr[vecScr]
      nScr <- length(idsScr)
      syph.scr[idsScr] <- 1
      syph.trt[idsScr] <- 1
      dat$attr$trtTime[idsScr] <- at
      dat$attr$scrTime[idsScr] <- at
    }
  }
  
  
  idsRec1 <- which(active == 1 & syph.stage < 4 & syph.trt == 1 & 
                    dat$attr$trtTime < at - 1)
    dat$attr$status[idsRec1] <- "s"
    syph.stage[idsRec1] <- 0
    syph.symp[idsRec1] <- 0
  
  idsRec2 <- which(active == 1 & syph.stage == 6 & syph.trt == 1 & 
                    dat$attr$trtTime < at - 3)
    dat$attr$stage[idsRec2] <- "s"
    syph.stage[idsRec2] <- 0
    syph.symp[idsRec2] <- 0
  
  dat$epi$scr.flow[at] <- nScr
  dat$epi$scr.num[at] <- sum(active == 1 & syph.scr == 1)
  dat$epi$trt.num[at] <- sum(active == 1 & syph.trt == 1)
  return(dat)
}
