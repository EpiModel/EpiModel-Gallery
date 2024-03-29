
##
## Simple Cost-effectiveness Model (Simple SI model with cost and utility tracking)
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##
## Author: Gregory Knowlton (University of Minnesota)
## Date: October 2021
##

# Cost/Utility Tracking Module ------------------------------------------------

costeffect <- function(dat, at) {

  # If the start of the analytic time horizon has not been reached, skip module
  cea.start <- get_param(dat, "cea.start")
  if (at < cea.start) {
    return(dat)
  }

  # Import attributes
  active <- get_attr(dat, "active")
  age <- get_attr(dat, "age")
  status <- get_attr(dat, "status")

  # Import parameters
  sus.cost <- get_param(dat, "sus.cost")
  inf.cost <- get_param(dat, "inf.cost")
  sus.qaly <- get_param(dat, "sus.qaly")
  inf.qaly <- get_param(dat, "inf.qaly")
  age.decrement <- get_param(dat, "age.decrement")
  disc.rate <- get_param(dat, "disc.rate")
  inter.cost <- get_param(dat, "inter.cost", override.null.error = TRUE)
  inter.start <- get_param(dat, "inter.start", override.null.error = TRUE)
  end.horizon <- get_param(dat, "end.horizon", override.null.error = TRUE)

  # Identify relevant sub-populations
  idsSus <- which(active == 1 & status == "s")
  idsInf <- which(active == 1 & status == "i")

  # Account for costs of each sub-population
  pop.sus.cost <- length(idsSus) * sus.cost
  pop.inf.cost <- length(idsInf) * inf.cost

  # Account for intervention costs
  # Costs accrue uniformly from intervention start to end horizon
  if (!is.null(inter.start) && at < end.horizon) {
    inter.cost <- inter.cost / (end.horizon - inter.start)
  } else {
    inter.cost <- 0
  }

  # Account for effects (QALYs) of each sub-population
  # QALY parameters by health state and age decrements must be converted to weekly time-step
  pop.sus.qaly <- sum((age[idsSus] * age.decrement + sus.qaly) / 52, na.rm = TRUE)
  pop.inf.qaly <- sum((age[idsInf] * age.decrement + inf.qaly) / 52, na.rm = TRUE)

  # Aggregate costs and effects
  pop.cost <- pop.sus.cost + pop.inf.cost + inter.cost
  pop.qaly <- pop.sus.qaly + pop.inf.qaly

  # Discount aggregated costs and effects
  t <- at - cea.start
  pop.cost.disc <- pop.cost * (1 - disc.rate) ^ (t / 52)
  pop.qaly.disc <- pop.qaly * (1 - disc.rate) ^ (t / 52)

  ## Summary statistics ##
  dat <- set_epi(dat, "cost", at, pop.cost)
  dat <- set_epi(dat, "qaly", at, pop.qaly)
  dat <- set_epi(dat, "cost.disc", at, pop.cost.disc)
  dat <- set_epi(dat, "qaly.disc", at, pop.qaly.disc)

  return(dat)
}

# Updated Aging Module --------------------------------------------------------

aging <- function(dat, at) {

  # Update age on attr and also the network
  active <- get_attr(dat, "active")
  idsActive <- which(active == 1)
  age <- get_attr(dat, "age")
  age[idsActive] <- age[idsActive] + 1 / 52
  dat <- set_attr(dat, "age", age)

  ## Summary statistics ##
  dat <- set_epi(dat, "meanAge", at, mean(age[idsActive], na.rm = TRUE))

  return(dat)
}


# Updated Departure Module -----------------------------------------------------

dfunc <- function(dat, at) {

  ## Attributes
  active <- get_attr(dat, "active")
  active.s <- get_attr(dat, "active.s")
  exitTime <- get_attr(dat, "exitTime")
  age <- get_attr(dat, "age")

  ## Parameters
  death.rates <- get_param(dat, "death.rates")

  ## Query active
  idsElig <- which(active == 1)
  nElig <- length(idsElig)
  nDeaths <- 0

  if (nElig > 0) {

    ## Calculate age-specific departure rates for each eligible node
    ## Everyone older than 85 gets the final mortality
    whole_ages_of_elig <- pmin(ceiling(age[idsElig]), 101)
    death_rates_of_elig <- death.rates[whole_ages_of_elig]

    ## Simulate departure process
    vecDeaths <- which(rbinom(nElig, 1, 1 - exp(-death_rates_of_elig)) == 1)
    idsDeaths <- idsElig[vecDeaths]
    nDeaths <- length(idsDeaths)

    ## Update nodal attributes
    if (nDeaths > 0) {
      active.s[idsDeaths] <- 0
      active[idsDeaths] <- 0
      exitTime[idsDeaths] <- at
    }

    ## 65+ who did not die this time-step
    idsRetire <- setdiff(which(age >= 65 & active.s == 1), idsDeaths)
    active.s[idsRetire] <- 0
  }

  ## Reset attr
  dat <- set_attr(dat, "active.s", active.s)
  dat <- set_attr(dat, "active", active)
  dat <- set_attr(dat, "exitTime", exitTime)

  ## Summary statistics
  dat <- set_epi(dat, "d.flow", at, nDeaths)

  return(dat)
}


# Updated Arrivals Module ----------------------------------------------------

afunc <- function(dat, at) {

  # If the end horizon has been reached, skip arrival module
  end.horizon <- get_param(dat, "end.horizon")
  if (at >= end.horizon) {
    return(dat)
  }

  ## Parameters
  n <- sum(get_attr(dat, "active") == 1)
  a.rate <- get_param(dat, "arrival.rate")

  ## Process
  nArrivalsExp <- n * a.rate
  nArrivals <- rpois(1, nArrivalsExp)

  ## Update attributes
  if (nArrivals > 0) {
    dat <- append_core_attr(dat, at, nArrivals)
    dat <- append_attr(dat, "status", "s", nArrivals)
    dat <- append_attr(dat, "infTime", NA, nArrivals)
    dat <- append_attr(dat, "age", 16, nArrivals)
    dat <- append_attr(dat, "active.s", 1, nArrivals)
  }

  ## Summary statistics
  dat <- set_epi(dat, "a.flow", at, nArrivals)

  return(dat)
}

# Updated Infection Module ---------------------------------------------------

ifunc <- function(dat, at) {

  # If the end horizon has been reached, skip infection module
  end.horizon <- get_param(dat, "end.horizon", override.null.error = TRUE)
  if (at >= end.horizon) {
    return(dat)
  }

  active.s <- get_attr(dat, "active.s")
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  infTime <- get_attr(dat, "infTime")
  inf.prob <- get_param(dat, "inf.prob")
  act.rate <- get_param(dat, "act.rate")
  inter.eff <- get_param(dat, "inter.eff", override.null.error = TRUE)
  inter.start <- get_param(dat, "inter.start", override.null.error = TRUE)

  # Identify which individuals are still sexually active
  idsActive.s <- which(active.s == 1)
  idsInf <- which(active == 1 & status == "i")
  nActive <- sum(active == 1)
  nElig <- length(idsInf)
  nInf <- 0
  if (nElig > 0 && nElig < nActive) {
    del <- discord_edgelist(dat, at, network = 1)
    if (!(is.null(del))) {
      # Select only rows of discordant edge list with sexually active partners
      del <- del[which(del$sus %in% idsActive.s & del$inf %in% idsActive.s), ]
      del$infDur <- at - infTime[del$inf]
      del$infDur[del$infDur == 0] <- 1
      linf.prob <- length(inf.prob)
      del$transProb <- ifelse(del$infDur <= linf.prob, inf.prob[del$infDur], inf.prob[linf.prob])
      if (!is.null(inter.eff) && at >= inter.start) {
        # Apply reduction in transmission probability due to prophylaxis
        del$transProb <- del$transProb * (1 - inter.eff)
      }
      lact.rate <- length(act.rate)
      del$actRate <- ifelse(del$infDur <= lact.rate, act.rate[del$infDur],
                            act.rate[lact.rate])
      del$finalProb <- 1 - (1 - del$transProb)^del$actRate
      transmit <- rbinom(nrow(del), 1, del$finalProb)
      del <- del[which(transmit == 1), ]
      idsNewInf <- unique(del$sus)
      status <- get_attr(dat, "status")
      status[idsNewInf] <- "i"
      dat <- set_attr(dat, "status", status)
      infTime[idsNewInf] <- at
      dat <- set_attr(dat, "infTime", infTime)
      nInf <- length(idsNewInf)
    }
  }
  if (nInf > 0) {
    dat <- set_transmat(dat, del, at)
  }
  dat <- set_epi(dat, "si.flow", at, nInf)
  return(dat)
}
