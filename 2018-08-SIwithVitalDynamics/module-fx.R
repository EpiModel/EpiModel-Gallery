
##
## Simple SI Model with Vital Dynamics
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##
## Author: Samuel M. Jenness (Emory University)
## Date: August 2018
##


# Aging Module ----------------------------------------------------------------

aging <- function(dat, at) {
  # Increment age by 1/52 (one week) per timestep and track mean age

  age <- get_attr(dat, "age")
  age <- age + 1 / 52
  dat <- set_attr(dat, "age", age)

  ## Summary statistics ##
  dat <- set_epi(dat, "meanAge", at, mean(age, na.rm = TRUE))

  return(dat)
}


# Departure Module ------------------------------------------------------------

dfunc <- function(dat, at) {
  # Simulate age-specific mortality with optional disease-induced excess mortality.
  # Each active node faces a weekly departure probability drawn from an
  # age-specific rate vector; infected nodes have their rate multiplied by
  # departure.disease.mult.

  ## Attributes ##
  active <- get_attr(dat, "active")
  exitTime <- get_attr(dat, "exitTime")
  age <- get_attr(dat, "age")
  status <- get_attr(dat, "status")

  ## Parameters ##
  dep.rates <- get_param(dat, "departure.rates")
  dep.dis.mult <- get_param(dat, "departure.disease.mult")

  ## Query alive ##
  idsElig <- which(active == 1)
  nElig <- length(idsElig)
  nDepts <- 0

  if (nElig > 0) {

    # Map continuous age to the mortality rate vector index.
    # The rate vector has 86 elements: position 1 = age <1, position 2 = age 1,
    # ..., position 86 = age 85+. ceiling(age) converts continuous age to the
    # correct index; pmin caps at 86 so ages >85 get the final (highest) rate.
    whole_ages_of_elig <- pmin(ceiling(age[idsElig]), 86)
    departure_rates_of_elig <- dep.rates[whole_ages_of_elig]

    # Disease-induced excess mortality: infected individuals face a higher
    # departure rate, scaled by the disease multiplier parameter
    idsElig.inf <- which(status[idsElig] == "i")
    departure_rates_of_elig[idsElig.inf] <-
      departure_rates_of_elig[idsElig.inf] * dep.dis.mult

    ## Simulate departure process ##
    vecDepts <- which(rbinom(nElig, 1, departure_rates_of_elig) == 1)
    idsDepts <- idsElig[vecDepts]
    nDepts <- length(idsDepts)

    ## Update nodal attributes ##
    if (nDepts > 0) {
      active[idsDepts] <- 0
      exitTime[idsDepts] <- at
    }
  }

  ## Write out updated attributes ##
  dat <- set_attr(dat, "active", active)
  dat <- set_attr(dat, "exitTime", exitTime)

  ## Summary statistics ##
  dat <- set_epi(dat, "d.flow", at, nDepts)

  return(dat)
}


# Arrival Module ---------------------------------------------------------------

afunc <- function(dat, at) {
  # Simulate births as new susceptible nodes entering at age 0.
  # The arrival rate is calibrated to the mean departure rate so that
  # the population size remains approximately stable over time.

  ## Parameters ##
  n <- sum(get_attr(dat, "active") == 1)
  a.rate <- get_param(dat, "arrival.rate")

  ## Process: Poisson-distributed number of births ##
  nArrivalsExp <- n * a.rate
  nArrivals <- rpois(1, nArrivalsExp)

  ## Update attributes ##
  if (nArrivals > 0) {
    dat <- append_core_attr(dat, at, nArrivals)
    dat <- append_attr(dat, "status", "s", nArrivals)
    dat <- append_attr(dat, "infTime", NA, nArrivals)
    dat <- append_attr(dat, "age", 0, nArrivals)
  }

  ## Summary statistics ##
  dat <- set_epi(dat, "a.flow", at, nArrivals)

  return(dat)
}
