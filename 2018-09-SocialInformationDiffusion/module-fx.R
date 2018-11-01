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
  
    # Get discordant edgelist
    del <- discord_edgelist(dat, at)
    
    # If some discordant edges, then proceed
    if (!(is.null(del))) {
      # Infection probabilities
      del$count<-1
      subdel<-aggregate(count~sus,FUN=length,data=del)
      del<-merge(del, subdel, by = "sus", all=TRUE)
      ##First set everyone's transmission probability as 0
      del$transProb <- 0
      #If some susceptible nodes have more than minimum degree of edge 
      #with infected person, then set their transmission prob as transmission probility
      if (dim(del[which(del$count.y>=dat$param$mini.degree),])[1]!=0)
      del[which(del$count.y>=dat$param$mini.degree),]$transProb <- inf.prob
      
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