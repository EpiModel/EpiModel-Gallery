
##
## SEIR Model extension: Syphilis Progression Model
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##
## Authors: Samuel M. Jenness, Yuan Zhao
## Date: November 2018
##


# Replacement infection/transmission module -------------------------------

infect <- function(dat, at) {

  ## Uncomment this to function environment interactively
  # browser()

  ## Attributes ##
  active <- dat$attr$active
  status <- dat$attr$status
  ## Initiating a indicator of syphilis status ##
  if (at == 2) {
    dat$attr$syph.status <- rep(0, length(active))
    dat$attr$syph.status<-ifelse(dat$attr$status=="i",1,0)
  }
  syph.status<- dat$attr$syph.status
  
  
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
      del$transProb <- ifelse(syph.status[del$inf]<4,inf.prob1,inf.prob2)
      del$transProb <- ifelse(syph.status[del$inf]>4,0,del$transProb)
      
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
        syph.status[idsNewInf]<-1
        dat$attr$status[idsNewInf] <- "i"
        dat$attr$infTime[idsNewInf] <- at
        
        
      }
      
      
    }
    
  }
  ## Save summary statistic for S->i flow
  dat$epi$si.flow[at] <- nInf
  dat$attr$syph.status<-syph.status

  return(dat)
}


# New disease progression module ------------------------------------------
# (Replaces the recovery module)

progress <- function(dat, at) {

  ## Uncomment this to function environment interactively
  # browser()

  ## Attributes ##
  active <- dat$attr$active
  syph.status <- dat$attr$syph.status

  ## Parameters ##
  ipr.rate <- dat$param$ipr.rate
  prse.rate <- dat$param$prse.rate
  seel.rate <- dat$param$seel.rate
  elll.rate <- dat$param$elll.rate
  llter.rate <- dat$param$llter.rate
  
  ## Incubation to primary stage progression process ##
  nPrim <- 0
  idsEligInf <- which(active == 1 & syph.status == 1)
  nEligInf <- length(idsEligInf)

  if (nEligInf > 0) {
    vecInf <- which(rbinom(nEligInf, 1, ipr.rate) == 1)
    if (length(vecInf) > 0) {
      idsInf <- idsEligInf[vecInf]
      nPrim  <- length(idsInf)
      syph.status[idsInf] <- 2
    }
  }

  ## Primary to secondary stage progression ##
  nSec <- 0
  idsEligRec <- which(active == 1 & syph.status ==2)
  nEligRec <- length(idsEligRec)

  if (nEligRec > 0) {
    vecRec <- which(rbinom(nEligRec, 1, prse.rate) == 1)
    if (length(vecRec) > 0) {
      idsRec <- idsEligRec[vecRec]
      nSec <- length(idsRec)
      syph.status[idsRec] <- 3
    }
  }
  
  ## Secondary to early latent progression ##
  nEarL <- 0
  idsEligRec <- which(active == 1 & syph.status == 3)
  nEligRec <- length(idsEligRec)
  
  if (nEligRec > 0) {
    vecRec <- which(rbinom(nEligRec, 1, seel.rate) == 1)
    if (length(vecRec) > 0) {
      idsRec <- idsEligRec[vecRec]
      nEarL <- length(idsRec)
      syph.status[idsRec] <- 4
    }
  }
  
  ## Early latent to late latent progression ##
  nLaL <- 0
  idsEligRec <- which(active == 1 & syph.status == 4)
  nEligRec <- length(idsEligRec)
  
  if (nEligRec > 0) {
    vecRec <- which(rbinom(nEligRec, 1, elll.rate) == 1)
    if (length(vecRec) > 0) {
      idsRec <- idsEligRec[vecRec]
      nLaL <- length(idsRec)
      syph.status[idsRec] <- 5
    }
  }
  
  ## Late latent to Tertiary progression ##
  nTer <- 0
  idsEligRec <- which(active == 1 & syph.status == 5)
  nEligRec <- length(idsEligRec)
  
  if (nEligRec > 0) {
    vecRec <- which(rbinom(nEligRec, 1, llter.rate) == 1)
    if (length(vecRec) > 0) {
      idsRec <- idsEligRec[vecRec]
      nTer <- length(idsRec)
      syph.status[idsRec] <- 6
    }
  }

  ## Write out updated status attribute ##
  dat$attr$syph.status <- syph.status

  ## Save summary statistics ##
  dat$epi$ipr.flow[at] <- nPrim
  dat$epi$prse.flow[at] <- nSec
  dat$epi$seel.flow[at] <- nEarL
  dat$epi$elll.flow[at] <- nLaL
  dat$epi$llter.flow[at] <- nTer
  
  dat$epi$inc.num[at] <- sum(active == 1 & syph.status == 1)
  dat$epi$pr.num[at] <- sum(active == 1 & syph.status == 2)
  dat$epi$se.num[at] <- sum(active == 1 & syph.status == 3)
  dat$epi$el.num[at] <- sum(active == 1 & syph.status == 4)
  dat$epi$ll.num[at] <- sum(active == 1 & syph.status == 5)
  dat$epi$ter.num[at] <- sum(active == 1 & syph.status ==6)
  
  return(dat)
}


