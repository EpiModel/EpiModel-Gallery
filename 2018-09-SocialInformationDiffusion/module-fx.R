# Update Transmission Module ----------------------------------------------

infect_mod <- function(dat, at) {
      
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
    #browser()
    # Get discordant edgelist
    del <- discord_edgelist(dat, at)
    
    # If some discordant edges, then proceed
    if (!(is.null(del))) {
      # Infection probabilities
      
      del$degree<-1
      subdel<-aggregate(degree~sus,FUN=length,data=del)
      del<-merge(del[, 1:3], subdel, by = "sus", all=TRUE)
      
      #If some susceptible nodes have more than minimum degree with infected person, 
      #then set their transmission prob as transmission probility
      #browser()
      ##Test if work for no one having more than min degree
      del$transProb <- ifelse(del$degree >= dat$param$min.degree, inf.prob, 0)
      
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


infect_newmod <- function(dat, at) {
  
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
    #browser()
    # Get discordant edgelist
    del <- discord_edgelist(dat, at)
    
    # If some discordant edges, then proceed
    if (!(is.null(del))) {
      # get discordant degree for susceptible individuals
      del$degree<-1
      subdel<-aggregate(degree~sus,FUN=length,data=del)
      del<-merge(del[, 1:3], subdel, by = "sus", all=TRUE)
      #The probability of infection is logistic function of susceptible nodes' degree with infected nodes
      #With parameters of a and b
      browser()
      # Two ways of modeling prob, one is binomial logistic, more scientific sense?
      # Another is the direct logistic model easier?
      # del$transProb <- 1/(1+dat$param$log_a*dat$param$log_b^del$degree)
      del$transProb <- exp(dat$param$log_a*dal$degree-dat$param$log_b)/(1+exp(dat$param$log_a*dal$degree-dat$param$log_b))
      
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