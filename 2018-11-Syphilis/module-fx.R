##
## SEIR Model extension: Syphilis Progression Model
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##
## Authors: Samuel M. Jenness, Yuan Zhao
## Date: November 2018
##


## Natural history of syphilis with testing and treatment
# Replacement infection/transmission module -------------------------------

infect_natural <- function(dat, at) {

  ## Uncomment this to function environment interactively
  # browser()

  ## Attributes ##
  active <- dat$attr$active
  status <- dat$attr$status
  ## Initiating a indicator of syphilis status to indicate stage of syphilis##
  
  if (at == 2) {
    dat$attr$syph.status <- rep(0, length(active))
    dat$attr$syph.status <- ifelse(dat$attr$status == "i",1,0)
    
    dat$attr$syph.symp <- rep(0, length(active))
    dat$attr$syph.trt <- rep(0, length(active))
    dat$attr$syph.scr <- rep(0, length(active))
    
    dat$attr$infTime <- rep(NA, length(active))
    dat$attr$infTime <- ifelse(dat$attr$status == "i",1,dat$attr$infTime)
    dat$attr$priTime <- rep(NA, length(active))
    dat$attr$secTime <- rep(NA, length(active))
    dat$attr$elTime <- rep(NA, length(active))
    dat$attr$llTime <- rep(NA, length(active))
    dat$attr$terTime <- rep(NA, length(active))
    dat$attr$trtTime <- rep(NA, length(active))
    dat$attr$scrTime <- rep(NA, length(active))
    
    dat$attr$syph.dur <-  rep(NA, length(active))
    dat$attr$syph2.dur <-  rep(NA, length(active))
    dat$attr$syph3.dur <-  rep(NA, length(active))
    dat$attr$syph4.dur <-  rep(NA, length(active))
    dat$attr$syph5.dur <-  rep(NA, length(active))
    dat$attr$syph6.dur <-  rep(NA, length(active))
    
  }
  
  syph.status <- dat$attr$syph.status
  
  
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
      del$transProb <- ifelse(syph.status[del$inf] < 4,inf.prob1,inf.prob2)
      del$transProb <- ifelse(syph.status[del$inf] > 4,0,del$transProb)
      
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
        syph.status[idsNewInf] <- 1
        dat$attr$status[idsNewInf] <- "i"
        dat$attr$infTime[idsNewInf] <- at
        
      }
      
      
    }
    
  }
  ## Save summary statistic for S->i flow
  dat$epi$si.flow[at] <- nInf
  dat$attr$syph.status <- syph.status
  dat$attr$syph.dur <- ifelse(dat$attr$syph.status == 1,(at - dat$attr$infTime),dat$attr$syph.dur)
  dat$epi$syph.dur[at] <- mean(dat$attr$syph.dur,na.rm = TRUE)

  return(dat)
}


# Syphilis progression module ------------------------------------------
# (Replaces the recovery module)

progress <- function(dat, at) {

  ## Uncomment this to function environment interactively
  ## browser()

  ## Attributes: syphilis stage, symptomatic or not, treatment and screening inidicators ##
  active <- dat$attr$active
  syph.status <- dat$attr$syph.status
  syph.symp <- dat$attr$syph.symp
  syph.trt <- dat$attr$syph.trt
  syph.scr <- dat$attr$syph.scr

  ## Parameters of stage transition ##
  ipr.rate <- dat$param$ipr.rate
  prse.rate <- dat$param$prse.rate
  seel.rate <- dat$param$seel.rate
  elll.rate <- dat$param$elll.rate
  llter.rate <- dat$param$llter.rate
  
  ## Parameters of symptomatic ##
  pri.sym <- dat$param$pri.sym 
  sec.sym <- dat$param$sec.sym
  
  ## Parameters of treatment and screening ##
  early.trt <- dat$param$early.trt 
  scr.rate <- dat$param$scr.rate
  
  ## Incubation to primary stage progression process ##
  ## browser()
  nPrim <- 0
  idsEligInf <- which(active == 1 & syph.status == 1)
  nEligInf <- length(idsEligInf)

  if (nEligInf > 0) {
    vecInf <- which(rbinom(nEligInf, 1, ipr.rate) == 1)
    if (length(vecInf) > 0) {
      idsInf <- idsEligInf[vecInf]
      nPrim  <- length(idsInf)
      syph.status[idsInf] <- 2
      dat$attr$priTime[idsInf] <- at
    }
  }
  dat$attr$syph.status <- syph.status
  dat$attr$syph2.dur <- ifelse(dat$attr$syph.status == 2,(at - dat$attr$priTime),dat$attr$syph2.dur)
  
  ## Primary stage symptomatic progression ##
  nPrim.sym <- 0
  idsEligSym <- which(active == 1 & syph.status == 2 & dat$attr$syph2.dur == 0)
  nEligSym <- length(idsEligSym)
  
  if (nEligSym > 0) {
    vecSym <- which(rbinom(nEligSym, 1, pri.sym) == 1)
    if (length(vecSym) > 0) {
      idsSym <- idsEligSym[vecSym]
      nPrim.sym <- length(idsSym)
      syph.symp[idsSym] <- 1
    }
  }
  dat$attr$syph.symp <- syph.symp
  
  ## Primary symptomatic patients receiving treatment
  nPrim.trt <- 0
  idsEligTrt <- which(active == 1 & syph.status == 2 & syph.symp == 1 & syph.trt == 0 & dat$attr$priTime < at)
  nEligTrt <- length(idsEligTrt)
  
  if (nEligTrt > 0) {
    vecTrt <- which(rbinom(nEligTrt, 1, early.trt) == 1)
    if (length(vecTrt) > 0) {
      idsTrt <- idsEligInf[vecTrt]
      nPrim.trt  <- length(idsTrt)
      syph.trt[idsTrt] <- 1
      dat$attr$trtTime[idsTrt] <- at
    }
  }
  dat$attr$syph.trt <- syph.trt
  
  idsRec <- which(active == 1 & syph.status == 2 & syph.trt == 1 & dat$attr$trtTime < at - 1)
  {
    dat$attr$status[idsRec] <- "s"
    syph.status[idsRec] <- 0
    syph.symp[idsRec] <- 0
    }
  
  

  ## Primary to secondary stage progression ##
  nSec <- 0
  idsEligRec <- which(active == 1 & syph.status == 2 & dat$attr$priTime < at)
  nEligRec <- length(idsEligRec)

  if (nEligRec > 0) {
    vecRec <- which(rbinom(nEligRec, 1, prse.rate) == 1)
    if (length(vecRec) > 0) {
      idsRec <- idsEligRec[vecRec]
      nSec <- length(idsRec)
      syph.status[idsRec] <- 3
      syph.symp[idsRec] <- 0
      dat$attr$secTime[idsRec] <- at
    }
  }
  dat$attr$syph.status <- syph.status
  dat$attr$syph3.dur <- ifelse(dat$attr$syph.status == 3,(at - dat$attr$secTime),dat$attr$syph3.dur)
  
  ## Secondary stage symptomatic progression
  nSec.sym <- 0
  idsEligSym2 <- which(active == 1 & syph.status == 3 & dat$attr$syph3.dur == 0)
  nEligSym2 <- length(idsEligSym2)
  
  if (nEligSym2 > 0) {
    vecSym2 <- which(rbinom(nEligSym2, 1, sec.sym) == 1)
    if (length(vecSym2) > 0) {
      idsSym2 <- idsEligSym2[vecSym2]
      nSec.sym <- length(idsSym2)
      syph.symp[idsSym2] <- 1
    }
  }
  dat$attr$syph.symp <- syph.symp
  
  ## Secondary symptomatic receiving treatment
  nSec.trt <- 0
  idsEligTrt <- which(active == 1 & syph.status == 3 & syph.symp == 1 & syph.trt == 0 & dat$attr$secTime < at)
  nEligTrt <- length(idsEligTrt)
  
  if (nEligTrt > 0) {
    vecTrt <- which(rbinom(nEligTrt, 1, early.trt) == 1)
    if (length(vecTrt) > 0) {
      idsTrt <- idsEligInf[vecTrt]
      nSec.trt  <- length(idsTrt)
      syph.trt[idsTrt] <- 1
      dat$attr$trtTime[idsTrt] <- at
    }
  }
  dat$attr$syph.trt <- syph.trt
  
  idsRec <- which(active == 1 & syph.status == 3 & syph.trt == 1 & dat$attr$trtTime < at - 1)
  {
    dat$attr$status[idsRec] <- "s"
    syph.status[idsRec] <- 0
    syph.symp[idsRec] <- 0
  }
  
  
  ## Secondary to early latent progression ##
  nEarL <- 0
  idsEligRec <- which(active == 1 & syph.status == 3 & dat$attr$secTime < at)
  nEligRec <- length(idsEligRec)
  
  if (nEligRec > 0) {
    vecRec <- which(rbinom(nEligRec, 1, seel.rate) == 1)
    if (length(vecRec) > 0) {
      idsRec <- idsEligRec[vecRec]
      nEarL <- length(idsRec)
      syph.status[idsRec] <- 4
      syph.symp[idsRec] <- 0
      dat$attr$elTime[idsRec] <- at
    }
  }
  dat$attr$syph.status <- syph.status
  dat$attr$syph.symp <- syph.symp
  dat$attr$syph4.dur <- ifelse(dat$attr$syph.status == 4,(at - dat$attr$elTime),dat$attr$syph4.dur)
  
  ## Early latent to late latent progression ##
  nLaL <- 0
  idsEligRec <- which(active == 1 & syph.status == 4 & dat$attr$elTime < at)
  nEligRec <- length(idsEligRec)
  
  if (nEligRec > 0) {
    vecRec <- which(rbinom(nEligRec, 1, elll.rate) == 1)
    if (length(vecRec) > 0) {
      idsRec <- idsEligRec[vecRec]
      nLaL <- length(idsRec)
      syph.status[idsRec] <- 5
      syph.symp[idsRec] <- 0
      dat$attr$llTime[idsRec] <- at
    }
  }
  dat$attr$syph.status <- syph.status
  dat$attr$syph5.dur <- ifelse(dat$attr$syph.status == 5,(at - dat$attr$llTime),dat$attr$syph5.dur)
  
  ## Late latent to Tertiary progression ##
  nTer <- 0
  idsEligRec <- which(active == 1 & syph.status == 5 & dat$attr$llTime < at)
  nEligRec <- length(idsEligRec)
  
  if (nEligRec > 0) {
    vecRec <- which(rbinom(nEligRec, 1, llter.rate) == 1)
    if (length(vecRec) > 0) {
      idsRec <- idsEligRec[vecRec]
      nTer <- length(idsRec)
      syph.status[idsRec] <- 6
      syph.symp[idsRec] <- 1
      syph.trt[idsRec] <- 1
      dat$attr$terTime[idsRec] <- at
      dat$attr$trtTime[idsRec] <- at
    }
  }

  dat$attr$syph.status <- syph.status
  dat$attr$syph.symp <- syph.symp
  dat$attr$syph6.dur <- ifelse(dat$attr$syph.status == 6,(at - dat$attr$terTime),dat$attr$syph6.dur)
  
  idsRec <- which(active == 1 & syph.status == 6 & syph.trt == 1 & dat$attr$trtTime < at - 3)
  {
    dat$attr$status[idsRec] <- "s"
    syph.status[idsRec] <- 0
    syph.symp[idsRec] <- 0
  }
  dat$attr$syph.status <- syph.status
  dat$attr$syph.symp <- syph.symp
  
  ## Screening of asymptomatic population
  nScr <- 0
  #idsEligScr <- which((active == 1 & syph.scr ==0 & syph.symp == 0)|(active == 1 & dat$attr$scrTime == at-13))
  idsEligScr <- which(active == 1 & syph.scr ==0 & syph.symp == 0)
  nEligScr <- length(idsEligScr)
  
  if (nEligScr > 0) {
    vecScr <- which(rbinom(nEligScr, 1, scr.rate) == 1)
    if (length(vecScr) > 0) {
      idsScr <- idsEligRec[vecScr]
      nScr <- length(idsScr)
      syph.scr[idsScr] <- 1
      syph.trt[idsScr] <- 1
      dat$attr$trtTime[idsScr] <- at
      dat$attr$scrTime[idsScr] <- at
    }
  }
  
  idsRec <- which(active == 1 & syph.status == 6 & syph.trt == 1 & dat$attr$trtTime < at - 3)
  {
    dat$attr$status[idsRec] <- "s"
    syph.status[idsRec] <- 0
    syph.symp[idsRec] <- 0
  }
  
  idsRec <- which(active == 1 & syph.status < 4 & syph.trt == 1 & dat$attr$trtTime < at - 1)
  {
    dat$attr$status[idsRec] <- "s"
    syph.status[idsRec] <- 0
    syph.symp[idsRec] <- 0
  }
  
  dat$attr$syph.scr <- syph.scr
  dat$attr$syph.trt <- syph.trt
  dat$attr$syph.status <- syph.status
  dat$attr$syph.symp <- syph.symp
  

  ## Save summary statistics ##
  dat$epi$ipr.flow[at] <- nPrim
  dat$epi$prse.flow[at] <- nSec
  dat$epi$seel.flow[at] <- nEarL
  dat$epi$elll.flow[at] <- nLaL
  dat$epi$llter.flow[at] <- nTer
  dat$epi$scr.flow[at] <- nScr
 
  dat$epi$inc.num[at] <- sum(active == 1 & syph.status == 1)
  dat$epi$pr.num[at] <- sum(active == 1 & syph.status == 2)
  dat$epi$se.num[at] <- sum(active == 1 & syph.status == 3)
  dat$epi$el.num[at] <- sum(active == 1 & syph.status == 4)
  dat$epi$ll.num[at] <- sum(active == 1 & syph.status == 5)
  dat$epi$ter.num[at] <- sum(active == 1 & syph.status == 6)
  dat$epi$sym.num[at] <- sum(active == 1 & syph.symp == 1)
  dat$epi$scr.num[at] <- sum(active == 1 & syph.scr == 1)
  dat$epi$trt.num[at] <- sum(active == 1 & syph.trt == 1)
  
  dat$epi$syph2.dur[at] <- mean(dat$attr$syph2.dur,na.rm = TRUE)
  dat$epi$syph3.dur[at] <- mean(dat$attr$syph3.dur,na.rm = TRUE)
  dat$epi$syph4.dur[at] <- mean(dat$attr$syph4.dur,na.rm = TRUE)
  dat$epi$syph5.dur[at] <- mean(dat$attr$syph5.dur,na.rm = TRUE)
  dat$epi$syph6.dur[at] <- mean(dat$attr$syph6.dur,na.rm = TRUE)
  
  
  
  return(dat)
}


