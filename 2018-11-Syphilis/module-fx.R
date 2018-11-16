
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

  ## Parameters ##
  inf.prob <- dat$param$inf.prob
  act.rate <- dat$param$act.rate

  ## Find infected nodes ##
  idsInf <- which(active == 1 & (status == "i"|status == "pr"|status == "se"|status == "el"))
  nActive <- sum(active == 1)
  nElig <- length(idsInf)

  ## Initialize default incidence at 0 ##
  nInf <- 0

  ## If any infected nodes, proceed with transmission ##
  if (nElig > 0 && nElig < nActive) {
    browser()
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
        dat$attr$status[idsNewInf] <- "i"
        dat$attr$infTime[idsNewInf] <- at
      }
    }
  }

  ## Save summary statistic for S->i flow
  dat$epi$si.flow[at] <- nInf

  return(dat)
}


# New disease progression module ------------------------------------------
# (Replaces the recovery module)

progress <- function(dat, at) {

  ## Uncomment this to function environment interactively
  # browser()

  ## Attributes ##
  active <- dat$attr$active
  status <- dat$attr$status

  ## Parameters ##
  ipr.rate <- dat$param$ipr.rate
  prse.rate <- dat$param$prse.rate
  seel.rate <- dat$param$seel.rate
  elll.rate <- dat$param$elll.rate
  llter.rate <- dat$param$llter.rate
  

  ## Incubation to primary stage progression process ##
  nPrim <- 0
  idsEligInf <- which(active == 1 & status == "i")
  nEligInf <- length(idsEligInf)

  if (nEligInf > 0) {
    vecInf <- which(rbinom(nEligInf, 1, ipr.rate) == 1)
    if (length(vecInf) > 0) {
      idsInf <- idsEligInf[vecInf]
      nPrim  <- length(idsInf)
      status[idsInf] <- "pr"
    }
  }

  ## Primary to secondary stage progression ##
  nSec <- 0
  idsEligRec <- which(active == 1 & status == "pr")
  nEligRec <- length(idsEligRec)

  if (nEligRec > 0) {
    vecRec <- which(rbinom(nEligRec, 1, prse.rate) == 1)
    if (length(vecRec) > 0) {
      idsRec <- idsEligRec[vecRec]
      nSec <- length(idsRec)
      status[idsRec] <- "se"
    }
  }
  
  ## Secondary to early latent progression ##
  nEarL <- 0
  idsEligRec <- which(active == 1 & status == "se")
  nEligRec <- length(idsEligRec)
  
  if (nEligRec > 0) {
    vecRec <- which(rbinom(nEligRec, 1, seel.rate) == 1)
    if (length(vecRec) > 0) {
      idsRec <- idsEligRec[vecRec]
      nEarL <- length(idsRec)
      status[idsRec] <- "el"
    }
  }
  
  ## Early latent to late latent progression ##
  nLaL <- 0
  idsEligRec <- which(active == 1 & status == "el")
  nEligRec <- length(idsEligRec)
  
  if (nEligRec > 0) {
    vecRec <- which(rbinom(nEligRec, 1, elll.rate) == 1)
    if (length(vecRec) > 0) {
      idsRec <- idsEligRec[vecRec]
      nLaL <- length(idsRec)
      status[idsRec] <- "ll"
    }
  }
  
  ## Late latent to Tertiary progression ##
  nTer <- 0
  idsEligRec <- which(active == 1 & status == "ll")
  nEligRec <- length(idsEligRec)
  
  if (nEligRec > 0) {
    vecRec <- which(rbinom(nEligRec, 1, ir.rate) == 1)
    if (length(vecRec) > 0) {
      idsRec <- idsEligRec[vecRec]
      nTer <- length(idsRec)
      status[idsRec] <- "ter"
    }
  }

  ## Write out updated status attribute ##
  dat$attr$status <- status

  ## Save summary statistics ##
  dat$epi$ipr.flow[at] <- nInf
  dat$epi$prse.flow[at] <- nSec
  dat$epi$seel.flow[at] <- nEarL
  dat$epi$elll.flow[at] <- nLaL
  dat$epi$llter.flow[at] <- nTer
  
  dat$epi$pr.num[at] <- sum(active == 1 & status == "pr")
  dat$epi$se.num[at] <- sum(active == 1 & status == "se")
  dat$epi$el.num[at] <- sum(active == 1 & status == "el")
  dat$epi$ll.num[at] <- sum(active == 1 & status == "ll")
  dat$epi$ter.num[at] <- sum(active == 1 & status == "ter")
  
  return(dat)
}



# Extension #1: Adding an R --> S Transition (SEIRS) ----------------------

progress2 <- function(dat, at) {

  ## Uncomment this to function environment interactively
  # browser()

  ## Attributes ##
  active <- dat$attr$active
  status <- dat$attr$status

  ## Parameters ##
  ei.rate <- dat$param$ei.rate
  ir.rate <- dat$param$ir.rate
  rs.rate <- dat$param$rs.rate

  ## E to I progression process ##
  nInf <- 0
  idsEligInf <- which(active == 1 & status == "e")
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

  # ## R to S progression process ##
  nSus <- 0
  idsEligSus <- which(active == 1 & status == "r")
  nEligSus <- length(idsEligSus)

  if (nEligSus > 0) {
    vecSus <- which(rbinom(nEligSus, 1, rs.rate) == 1)
    if (length(vecSus) > 0) {
      idsSus <- idsEligSus[vecSus]
      nSus <- length(idsSus)
      status[idsSus] <- "s"
    }
  }

  ## Write out updated status attribute ##
  dat$attr$status <- status

  ## Save summary statistics ##
  dat$epi$ei.flow[at] <- nInf
  dat$epi$ir.flow[at] <- nRec
  dat$epi$rs.flow[at] <- nSus
  dat$epi$e.num[at] <- sum(active == 1 & status == "e")
  dat$epi$r.num[at] <- sum(active == 1 & status == "r")
  dat$epi$s.num[at] <- sum(active == 1 & status == "s")

  return(dat)
}
