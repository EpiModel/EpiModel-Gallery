# Scenario 1: Minimum degree transmission pattern
# Update Transmission Module ----------------------------------------------

diffuse_mod <- function(dat, at) {

  ## Attributes ##
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  if (at == 2) {
    infTime <- rep(NA, length(active))
    infTime[which(status == "i")] <- 1
    set_attr(dat, "infTime", infTime)
  } else {
    infTime <- get_attr(dat, "infTime")
  }


  ## Parameters ##
  inf.prob <- get_param(dat, "inf.prob")
  act.rate <- get_param(dat, "act.rate")
  min.degree <- get_param(dat, "min.degree")

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

      # Initialize degree variable on edgelist
      del$degree <- 1

      # Count degree in the del
      subdel <- aggregate(degree ~ sus, FUN = length, data = del)

      # Merge new degree count in
      del <- merge(del[, 1:3], subdel, by = "sus", all = TRUE)

      # If some susceptible nodes have more than minimum degree with infected person,
      # then set their transmission prob as transmission probility
      ## Test if work for no one having more than min degree
      del$transProb <- ifelse(del$degree >= min.degree, inf.prob, 0)

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
        dat <- set_attr(dat, "status", status)
        infTime[idsNewInf] <- at
        dat <- set_attr(dat, "infTime", infTime)
      }

    }
  }

  ## Summary statistics ##
  dat <- set_epi(dat, "si.flow", at, totInf)

  return(dat)
}

# Scenario 2: Transmission probability as a function of degree of discordant edgelist

diffuse_mod2 <- function(dat, at) {

  ## Attributes ##
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  if (at == 2) {
    infTime <- rep(NA, length(active))
    infTime[which(status == "i")] <- 1
    dat <- set_attr(dat, "infTime", infTime)
  } else {
    infTime <- get_attr(dat, "infTime")
  }


  ## Parameters ##
  beta0 <- get_param(dat, "beta0")
  beta1 <- get_param(dat, "beta1")
  act.rate <- get_param(dat, "act.rate")

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

      # Initialize degree variable on edgelist
      del$degree <- 1

      # Count degree in the del
      subdel <- aggregate(degree ~ sus, FUN = length, data = del)

      # Merge new degree count in
      del <- merge(del[, 1:3], subdel, by = "sus", all = TRUE)

      # The probability of infection is logistic function of susceptible nodes' degree with infected nodes
      # With parameters of beta0 and beta1
      ## beta1: log odds ratio with 1 degree increase of discordant relationship
      ## beta0: baseline log odds of transmission when degree of discordant edgelist is 0
      del$transProb <- plogis(beta0 + beta1*del$degree)

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
        dat <- set_attr(dat, "status", status)
        infTime[idsNewInf] <- at
        dat <- set_attr(dat, "infTime", infTime)
      }

    }
  }

  ## Summary statistics ##
  dat <- set_epi(dat, "si.flow", at, totInf)

  return(dat)
}
