##
## Modeling Two Competing Strains in an SIS Epidemic
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##
## Author: Steven M. Goodreau (University of Washington)
## Date: September 2018
##

# Revised transmission function : 
        # any lines changed are commented with ## EDITED
        # any lines deleted are commented out and marked with ## DELETED
        # any lines added are commented with ## ADDED

infection.2strains <- function(dat, at) {
  
  ## Uncomment this to function environment interactively
  #browser()
  
  # Note: strain for the initial population are assigned in the 
  #   recovery module, since that is run first
  
  # Variables ---------------------------------------------------------------
  active <- dat$attr$active
  status <- dat$attr$status
  modes <- dat$param$modes
  mode <- idmode(dat$nw)
  
  inf.prob <- dat$param$inf.prob
  inf.prob.m2 <- dat$param$inf.prob.m2
  act.rate <- dat$param$act.rate
  inf.prob.st2 <- dat$param$inf.prob.st2  ## ADDED
  
  nw <- dat$nw
  tea.status <- dat$control$tea.status
  
  # Vector of infected and susceptible IDs
  idsInf <- which(active == 1 & status == "i")
  nActive <- sum(active == 1)
  nElig <- length(idsInf)
  
  # Initialize vectors
  nInf <- nInfST2 <- totInf <- 0  ## EDITED from m2 to st2
  
  ## ADDED: for this example, allowing for different transmission probabilities by 
  ## strain means suppressing heterogeniety by mode in transmission probabilities
  if (modes==2) stop("Error: modeling two pathogen strains is only impemented for 
                     one-mode networks")
  if (!is.null(inf.prob.m2)) stop("Error: inf.prob.m2 cannot be set when modeling 
                     two pathogen strains")
  
  # Process -----------------------------------------------------------------
  # If some infected AND some susceptible, then proceed
  if (nElig > 0 && nElig < nActive) {
    
    # Get discordant edgelist
    del <- discord_edgelist(dat, at)
    
    # If some discordant edges, then proceed
    if (!(is.null(del))) {
      
      # Infection duration to at
      del$infDur <- at - dat$attr$infTime[del$inf]
      del$infDur[del$infDur == 0] <- 1
      
      # Calculate infection-stage transmission rates
      linf.prob <- length(inf.prob)
      linf.prob.st2 <- length(inf.prob.st2)  ## ADDED
      
      ## DELETED
      #if (is.null(inf.prob.m2)) {
      #  del$transProb <- ifelse(del$infDur <= linf.prob,
      #                          inf.prob[del$infDur],
      #                          inf.prob[linf.prob])
      #} else {
      #  del$transProb <- ifelse(del$sus <= nw %n% "bipartite",
      #                          ifelse(del$infDur <= linf.prob,
      #                                 inf.prob[del$infDur],
      #                                 inf.prob[linf.prob]),
      #                          ifelse(del$infDur <= linf.prob,
      #                                 inf.prob.m2[del$infDur],
      #                                 inf.prob.m2[linf.prob]))
      #}
    
      ## ADDED
      del$transProb <- ifelse(dat$attr$strain[del$inf] == 1,
                                ifelse(del$infDur <= linf.prob,
                                       inf.prob[del$infDur],
                                       inf.prob[linf.prob]),
                                ifelse(del$infDur <= linf.prob.st2,
                                       inf.prob.st2[del$infDur],
                                       inf.prob.st2[linf.prob]))
      
      # Interventions
      if (!is.null(dat$param$inter.eff) && at >= dat$param$inter.start) {
        del$transProb <- del$transProb * (1 - dat$param$inter.eff)
      }
      
      # Calculate infection-stage act/contact rates
      lact.rate <- length(act.rate)
      del$actRate <- ifelse(del$infDur <= lact.rate,
                            act.rate[del$infDur],
                            act.rate[lact.rate])
      
      # Calculate final transmission probability per timestep
      del$finalProb <- 1 - (1 - del$transProb) ^ del$actRate
      
      # Randomize transmissions and subset df
      transmit <- rbinom(nrow(del), 1, del$finalProb)
      del <- del[which(transmit == 1), ]
      
      # Set new infections vector
      idsNewInf <- unique(del$sus)
      
      ## ADDED : determine which pairs actually involved infection
      ## This line is needed because it is possible for someone to be flagged for infection
      ## by more than one partner in the same time step; we need to make sure to only assign 
      ## them the strain of the first person who infects them
      infpairs <- sapply(idsNewInf, function(x) min(which(del$sus==x)))
      
      if(length(infpairs)>0) { 
        ## ADDED : assign the strain from the infecting partner to the newly infected partner
        infectors <- del$inf[infpairs]
        infectors_strain <- dat$attr$strain[infectors]
         
        ## EDITED from distinguishing between m1 and m2 to st1 and st2
        nInf <- sum(infectors_strain == 1)
        nInfST2 <- sum(infectors_strain == 2)
        totInf <- nInf + nInfST2
      } else {
        nInf <- nInfST2 <- totInf <- 0
      }
          
      # Update nw attributes
      if (totInf > 0) {
        if (tea.status == TRUE) {
          nw <- activate.vertex.attribute(nw,
                                          prefix = "testatus",
                                          value = "i",
                                          onset = at,
                                          terminus = Inf,
                                          v = idsNewInf)
        }
        dat$attr$status[idsNewInf] <- "i"
        dat$attr$infTime[idsNewInf] <- at
        
        
        dat$attr$strain[idsNewInf] <- dat$attr$strain[infectors] ##ADDED
        
        form <- get_nwparam(dat)$formation
        fterms <- get_formula_terms(form)
        if ("status" %in% fterms) {
          nw <- set.vertex.attribute(nw, "status", dat$attr$status)
        }
      }
      
      # Substitute PIDs for vital bipartite sims
      if (any(names(nw$gal) %in% "vertex.pid")) {
        del$sus <- get.vertex.pid(nw, del$sus)
        del$inf <- get.vertex.pid(nw, del$inf)
      }
      
    } # end some discordant edges condition
  } # end some active discordant nodes condition
  
  
  # Output ------------------------------------------------------------------
  
  # Save transmission matrix
  if (totInf > 0) {
    del <- del[!duplicated(del$sus), ]
    if (at == 2) {
      dat$stats$transmat <- del
    } else {
      dat$stats$transmat <- rbind(dat$stats$transmat, del)
    }
  }
  
  ## Save incidence vector
  if (at == 2) {
    dat$epi$si.flow <- c(0, nInf)
    dat$epi$si.flow.st2 <- c(0,nInfST2)       ## ADDED
    #if (modes == 2) {                        ## DELETED
    #  dat$epi$si.flow.m2 <- c(0, nInfM2)
    #}
  } else {
    dat$epi$si.flow[at] <- nInf
    dat$epi$si.flow.st2[at] <- nInfST2        ## ADDED
    #if (modes == 2) {                        ## DELETED
    #  dat$epi$si.flow.m2[at] <- nInfM2
    #}
  }

  ## ADDED
  dat$epi$i.num.st1[at] <- sum(dat$attr$status == "i" & dat$attr$strain == 1)
  dat$epi$i.num.st2[at] <- sum(dat$attr$status == "i" & dat$attr$strain == 2)
  
  dat$nw <- nw

  return(dat)
}



# Updated Recovery Module --------------------------------------------------

recov.2strains <- function(dat, at) {

  ## Adapted from the recov module in the EpiModel Gallery TestAndTreat example
  ## Uncomment this to function environment interactively
  #browser()
  
  active <- dat$attr$active
  status <- dat$attr$status

  ## ADDED
  # Initialize strain attribute ---------------------------------------------
  strain <- dat$attr$strain
  pct.st2 <- dat$param$pct.st2
  idsInf <- which(active == 1 & status == "i")
  nActive <- sum(active == 1)
  nElig <- length(idsInf)
  if(at == 2) {
    strain <- rep(NA, nActive)
    strain[idsInf] <- rbinom(nElig, 1, pct.st2) + 1  # Strains are labeled 1 and 2
    dat$attr$strain <- strain
  }
  
  ## Parameters ##
  rec.rate <- dat$param$rec.rate
  rec.rate.st2 <- dat$param$rec.rate.st2      ## EDITED

  ## Determine eligible to recover ##
  idsElig <- which(active == 1 & status == "i")
  nElig <- length(idsElig)

  ## Determine strain of eligible ##
  strain.elig <- strain[idsElig]              ## EDITED 

  ## Recovery rates dependent on strain ##
  ratesElig <- ifelse(strain.elig == 1, rec.rate, rec.rate.st2)   ## EDITED
  
  ## Vector of recovered IDs after stochastic process
  vecRecov <- which(rbinom(nElig, 1, ratesElig) == 1)
  idsRecov <- idsElig[vecRecov]
  nRecov <- length(idsRecov)
  nRecov.st1 <- sum(strain[idsRecov]==1)            ## ADDED
  nRecov.st2 <- sum(strain[idsRecov]==2)            ## ADDED
  
  ## Update attributes if any recovered ##
  status[idsRecov] <- "s"
  strain[idsRecov] <- NA                            ## EDITED
  
  ## Write out updated attributes ##
  dat$attr$status <- status
  dat$attr$strain <- strain                         ## EDITED
  
  ## Write out summary statistics ##
  dat$epi$is.flow[at] <- nRecov
  dat$epi$is.flow.st1[at] <- nRecov.st1             ## ADDED
  dat$epi$is.flow.st2[at] <- nRecov.st2             ## ADDED

  return(dat)
}
