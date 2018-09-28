# Update Transmission Module ----------------------------------------------

new_infect_mod <- function(dat, at) {
  
  ## Attributes ##
  active <- dat$attr$active
  status <- dat$attr$status
  
  ## Parameters ##
  inf.prob <- dat$param$inf.prob
  act.rate <- dat$param$act.rate
  
  # Vector of infected and susceptible IDs
  idsInf <- which(active == 1 & status == "i")
  nActive <- sum(active == 1)
  nElig <- length(idsInf)
  
  # Initialize vectors at 0
  nInf <- totInf <- 0
  
  ## Processes ##
  # If some infected AND some susceptible, then proceed
  if (nElig > 0 && nElig < nActive) {
    
    # Get discordant edgelist
    del <- discord_edgelist(dat, at)
    
    # If some discordant edges, then proceed
    if (!(is.null(del))&(get_degree(dat$nw)>=dat$param$minimun.degree)) {
      
      # Infection probabilities
      del<-nw[which(get_degree(nw)>3),]
      del$transProb <- inf.prob
      
      # Act rates
      del$actRate <- act.rate
      
      # Calculate final transmission probability per timestep
      del$finalProb <- 1 - (1 - del$transProb) ^ del$actRate
      
      # Randomize transmissions and subset df
      transmit <- rbinom(nrow(del), 1, del$finalProb)
      del <- del[which(transmit == 1), ]
      
      # Set new infections vector
      idsNewInf <- unique(del$sus)
      totInf <- length(idsNewInf)
      
      # Update attributes for newly infected
      if (totInf > 0) {
        dat$attr$status[idsNewInf] <- "i"
        dat$attr$infTime[idsNewInf] <- at
      }
      
    }
  }
  
  ## Summary statistics ##
  dat$epi$si.flow[at] <- totInf
  
  return(dat)
}