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
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")

  inf.prob <- get_param(dat, "inf.prob")
  act.rate <- get_param(dat, "act.rate")
  inf.prob.st2 <- get_param(dat, "inf.prob.st2")

  # Vector of infected and susceptible IDs
  idsInf <- which(active == 1 & status == "i")
  nActive <- sum(active == 1)
  nElig <- length(idsInf)

  # Initialize vectors
  nInf <- nInfST2 <- totInf <- 0

  ## Process ##

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
      infpairs <- sapply(idsNewInf, function(x) min(which(del$sus == x)))

      if (length(infpairs) > 0) {
        ## ADDED : assign the strain from the infecting partner to the newly infected partner
        infectors <- del$inf[infpairs]
        strain <- get_attr(dat, "strain")
        infectors_strain <- strain[infectors]

        ## EDITED from distinguishing between m1 and m2 to st1 and st2
        nInf <- sum(infectors_strain == 1)
        nInfST2 <- sum(infectors_strain == 2)
        totInf <- nInf + nInfST2
      } else {
        nInf <- nInfST2 <- totInf <- 0
      }

      # Update nw attributes
      if (totInf > 0) {
        status[idsNewInf] <- "i"
        dat <- set_attr(dat, "status", status)
        infTime[idsNewInf] <- at
        dat <- set_attr(dat, "infTime", infTime)
        strain[idsNewInf] <- strain[infectors]
        dat <- set_attr(dat, "strain", strain)
      }

    }
  }

  ## Output ##

  ## Save incidence vector
  dat$epi$si.flow[at] <- nInf
  dat$epi$si.flow.st2[at] <- nInfST2

  ## Save prevalence vector
  dat$epi$i.num.st1[at] <- sum(dat$attr$status == "i" & dat$attr$strain == 1)
  dat$epi$i.num.st2[at] <- sum(dat$attr$status == "i" & dat$attr$strain == 2)

  return(dat)
}



# Updated Recovery Module --------------------------------------------------

recov.2strains <- function(dat, at) {

  ## Uncomment this to function environment interactively
  #browser()

  active <- dat$attr$active
  status <- dat$attr$status

  ## Initialize strain attribute
  strain <- dat$attr$strain
  pct.st2 <- dat$param$pct.st2
  idsInf <- which(active == 1 & status == "i")
  nActive <- sum(active == 1)
  nElig <- length(idsInf)
  if (at == 2) {
    strain <- rep(NA, nActive)
    strain[idsInf] <- rbinom(nElig, 1, pct.st2) + 1  # Strains are labeled 1 and 2
    dat$attr$strain <- strain
  }

  ## Parameters ##
  rec.rate <- dat$param$rec.rate
  rec.rate.st2 <- dat$param$rec.rate.st2

  ## Determine eligible to recover ##
  idsElig <- which(active == 1 & status == "i")
  nElig <- length(idsElig)

  ## Determine strain of eligible ##
  strain.elig <- strain[idsElig]

  ## Recovery rates dependent on strain ##
  ratesElig <- ifelse(strain.elig == 1, rec.rate, rec.rate.st2)

  ## Vector of recovered IDs after stochastic process
  vecRecov <- which(rbinom(nElig, 1, ratesElig) == 1)
  idsRecov <- idsElig[vecRecov]
  nRecov <- length(idsRecov)
  nRecov.st1 <- sum(strain[idsRecov] == 1)
  nRecov.st2 <- sum(strain[idsRecov] == 2)

  ## Update attributes if any recovered ##
  status[idsRecov] <- "s"
  strain[idsRecov] <- NA

  ## Write out updated attributes ##
  dat$attr$status <- status
  dat$attr$strain <- strain

  ## Write out summary statistics ##
  dat$epi$is.flow[at] <- nRecov
  dat$epi$is.flow.st1[at] <- nRecov.st1
  dat$epi$is.flow.st2[at] <- nRecov.st2

  return(dat)
}
