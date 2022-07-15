
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
  dat$nw[[1]] <- x
  dat <- set_param(dat, "groups", 1)

  # Epidemic parameters
  i.num <- get_init(dat, "i.num")

  ## Core attributes and Infection status attributes
  n <- network.size(dat$nw[[1]])
  dat <- append_core_attr(dat, 1, n)

  status <- rep("s", n)
  status[sample(1:n, i.num)] <- "i"
  dat <- set_attr(dat, "status", status)

  infTime <- rep(NA, n)
  infTime[which(status == "i")] <- 1
  dat <- set_attr(dat, "infTime", infTime)

  dat <- prevalence.net(dat, 1)
  return(dat)
}


## Update Transmission Module ##

new_infect_mod <- function(dat, at) {

  ## Attributes ##
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  infTime <- get_attr(dat, "infTime")

  ## Parameters ##
  inf.prob <- get_param(dat, "inf.prob")
  act.rate <- get_param(dat, "act.rate")

  # Vector of infected and susceptible IDs
  idsInf <- which(active == 1 & status == "i")
  nActive <- sum(active == 1)
  nElig <- length(idsInf)

  # Initialize vectors at 0
  totInf <- 0

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
        status[idsNewInf] <- "i"
        infTime[idsNewInf] <- at
        dat <- set_attr(dat, "status", status)
        dat <- set_attr(dat, "infTime", infTime)
      }

    }
  }

  ## Summary statistics ##
  dat <- set_epi(dat, "si.flow", at, totInf)

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
  dat$nw[[1]] <- x
  dat <- set_param(dat, "groups", 1)

  # Epidemic parameters
  i.num <- get_init(dat, "i.num")

  ## Core attributes and Infection status attributes
  n <- network.size(dat$nw[[1]])
  dat <- append_core_attr(dat, 1, n)

  status <- rep("s", n)
  status[sample(1:n, i.num)] <- "i"
  dat <- set_attr(dat, "status", status)

  infTime <- rep(NA, n)
  infTime[which(status == "i")] <- 1
  dat <- set_attr(dat, "infTime", infTime)

  # Set time-varying status attribute on network (for plotting)
  dat$nw[[1]] <- activate.vertex.attribute(dat$nw[[1]],
                                           prefix = "testatus",
                                           value = get_attr(dat, "status"),
                                           onset = 1,
                                           terminus = Inf)

  dat <- prevalence.net(dat, 1)

  return(dat)
}


## Update Transmission Module ##

new_infect_mod2 <- function(dat, at) {

  ## Attributes ##
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  infTime <- get_attr(dat, "infTime")

  ## Parameters ##
  inf.prob.stage1 <- get_param(dat, "inf.prob.stage1")
  inf.prob.stage2 <- get_param(dat, "inf.prob.stage2")
  dur.stage1 <- get_param(dat, "dur.stage1")
  act.rate <- get_param(dat, "act.rate")

  # Vector of infected and susceptible IDs
  idsInf <- which(active == 1 & status == "i")
  nActive <- sum(active == 1)
  nElig <- length(idsInf)

  # Initialize vectors at 0
  totInf <- 0

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
        status[idsNewInf] <- "i"
        infTime[idsNewInf] <- at
        dat <- set_attr(dat, "status", status)
        dat <- set_attr(dat, "infTime", infTime)

        # Update time-varying status on network (for plotting)
        dat$nw[[1]] <- activate.vertex.attribute(dat$nw[[1]],
                                                 prefix = "testatus",
                                                 value = "i",
                                                 onset = at,
                                                 terminus = Inf,
                                                 v = idsNewInf)
      }

    }
  }

  ## Summary statistics ##
  dat <- set_epi(dat, "si.flow", at, totInf)

  return(dat)
}
