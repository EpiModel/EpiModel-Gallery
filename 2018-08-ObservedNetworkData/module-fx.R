
##
## Modeling Epidemics over Observed Networks
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##
## Author: Samuel M. Jenness (Emory University)
## Date: August 2018
##


# Example 1: Base Model ---------------------------------------------------

## Updated Initialization Module ##

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


## Update Transmission Module ##

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



# Example 2: Adding Network Tracking and Time-Varying Risk ----------------


## Updated Initialization Module ##

new_init_mod2 <- function(x, param, init, control, s) {

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

  # Set time-varying status attribute on network (for plotting)
  dat$nw <- activate.vertex.attribute(dat$nw,
                                      prefix = "testatus",
                                      value = dat$attr$status,
                                      onset = 1,
                                      terminus = Inf)

  ## Get initial prevalence
  dat <- get_prev.net(dat, at = 1)

  return(dat)
}


## Update Transmission Module ##

new_infect_mod2 <- function(dat, at) {

  ## Attributes ##
  active <- dat$attr$active
  status <- dat$attr$status
  infTime <- dat$attr$infTime

  ## Parameters ##
  inf.prob.stage1 <- dat$param$inf.prob.stage1
  inf.prob.stage2 <- dat$param$inf.prob.stage2
  dur.stage1 <- dat$param$dur.stage1
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
      infDur.del <- at - infTime[del$inf]
      inf.prob.del <- ifelse(infDur.del <= dur.stage1, inf.prob.stage1, inf.prob.stage2)
      del$transProb <- inf.prob.del

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

        # Update time-varying status on network (for plotting)
        dat$nw <- activate.vertex.attribute(dat$nw,
                                            prefix = "testatus",
                                            value = "i",
                                            onset = at,
                                            terminus = Inf,
                                            v = idsNewInf)
      }

    }
  }

  ## Summary statistics ##
  dat$epi$si.flow[at] <- totInf

  return(dat)
}
