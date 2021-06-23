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

init_strain <- function(dat, at) {

  # Module only runs at initial time step
  if (at == 2) {

    # Pull attributes
    active <- get_attr(dat, "active")
    status <- get_attr(dat, "status")

    nActive <- sum(active == 1)
    idsInf <- which(active == 1 & status == "i")

    # Pull parameter
    pct.st2 <- get_param(dat, "pct.st2")

    # Set up strain attr, with NA as default
    strain <- rep(NA, nActive)

    # Strains are labeled 1 and 2
    strain[idsInf] <- rbinom(length(idsInf), 1, pct.st2) + 1

    # Set new attr on dat
    dat <- set_attr(dat, "strain", strain)
  }

  return(dat)
}

infection_2strains <- function(dat, at) {

  ## Uncomment this to function environment interactively
  #browser()

  # Note: strain for the initial population are assigned in the
  #   recovery module, since that is run first

  # Variables ---------------------------------------------------------------
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  infTime <- get_attr(dat, "infTime")
  strain <- get_attr(dat, "strain")

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
      del$infDur <- at - infTime[del$inf]
      del$infDur[del$infDur == 0] <- 1

      # Calculate infection-stage transmission rates
      linf.prob <- length(inf.prob)
      linf.prob.st2 <- length(inf.prob.st2)  ## ADDED

      ## ADDED
      del$transProb <- ifelse(strain[del$inf] == 1,
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
  dat <- set_epi(dat, "si.flow", at, nInf)
  dat <- set_epi(dat, "si.flow.st2", at,nInfST2)

  ## Save prevalence vector
  dat <- set_epi(dat, "i.num.st1", at, sum(status == "i" & strain == 1))
  dat <- set_epi(dat, "i.num.st2", at, sum(status == "i" & strain == 2))

  return(dat)
}



# Updated Recovery Module --------------------------------------------------

recov_2strains <- function(dat, at) {

  ## Uncomment this to function environment interactively
  #browser()

  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")

  ## Initialize strain attribute
  pct.st2 <- get_param(dat, "pct.st2")
  idsInf <- which(active == 1 & status == "i")
  nElig <- length(idsInf)
  strain <- get_attr(dat, "strain")

  ## Parameters ##
  rec.rate <- get_param(dat, "rec.rate")
  rec.rate.st2 <- get_param(dat,"rec.rate.st2")

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
  dat <- set_attr(dat, "status", status)
  dat <- set_attr(dat, "strain", strain)

  ## Write out summary statistics ##
  dat <- set_epi(dat, "is.flow", at, nRecov)
  dat <- set_epi(dat, "is.flow.st1", at, nRecov.st1)
  dat <- set_epi(dat, "is.flow.st2", at, nRecov.st2)

  return(dat)
}
