
##
## Modeling Epidemics over Observed Networks
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##
## Author: Samuel M. Jenness (Emory University)
## Date: August 2018
##


# Updated Initialization Module -------------------------------------------

new_init_mod <- function(x, param, init, control, s) {

  # Master Data List
  dat <- list()
  dat$param <- param
  dat$init <- init
  dat$control <- control

  dat$attr <- list()
  dat$stats <- list()
  dat$temp <- list()

  # Network parameters
  dat$nw <- x
  dat$param$modes <- 1

  ## Infection status and time attributes
  n <- network.size(dat$nw)
  dat$attr$status <- rep("s", n)
  dat$attr$status[sample(1:n, init$i.num)] <- "i"

  dat$attr$active <- rep(1, n)
  dat$attr$entrTime <- rep(1, n)
  dat$attr$exitTime <- rep(NA, n)

  dat$attr$infTime <- rep(NA, n)
  dat$attr$infTime[dat$attr$status == "i"] <- 1

  ## Get initial prevalence
  dat <- get_prev.net(dat, at = 1)

  return(dat)
}


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
    if (!(is.null(del))) {

      # Infection probabilities
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
