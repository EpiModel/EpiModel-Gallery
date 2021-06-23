
##
## Simple SI Model with Vital Dynamics
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##
## Author: Samuel M. Jenness (Emory University)
## Date: August 2018
##


# New Aging Module --------------------------------------------------------

aging <- function(dat, at) {

  # Update age on attr and also the network
  age <- get_attr(dat, "age")
  age <- age + 1/52
  dat <- set_attr(dat, "age", age)

  ## Summary statistics ##
  dat <- set_epi(dat, "meanAge", at, mean(age, na.rm = TRUE))

  return(dat)
}


# Update Departure Module -----------------------------------------------------

dfunc <- function(dat, at) {

  ## Attributes
  active <- get_attr(dat, "active")
  exitTime <- get_attr(dat, "exitTime")
  age <- get_attr(dat, "age")
  status <- get_attr(dat, "status")

  ## Parameters
  dep.rates <- get_param(dat, "departure.rates")
  dep.dis.mult <- get_param(dat, "departure.disease.mult")

  ## Query alive
  idsElig <- which(active == 1)
  nElig <- length(idsElig)
  nDepts <- 0

  if (nElig > 0) {

    ## Calculate age-specific departure rates for each eligible node
    ## Everyone older than 85 gets the final mortality
    whole_ages_of_elig <- pmin(ceiling(age[idsElig]), 86)
    departure_rates_of_elig <- dep.rates[whole_ages_of_elig]

    ## Multiply departure rates for diseased persons
    idsElig.inf <- which(status[idsElig] == "i")
    departure_rates_of_elig[idsElig.inf] <- departure_rates_of_elig[idsElig.inf] * dep.dis.mult

    ## Simulate departure process
    vecDepts <- which(rbinom(nElig, 1, departure_rates_of_elig) == 1)
    idsDepts <- idsElig[vecDepts]
    nDepts <- length(idsDepts)

    ## Update nodal attributes
    if (nDepts > 0) {
      active[idsDepts] <- 0
      exitTime[idsDepts] <- at
    }
  }

  ## Reset attr
  dat <- set_attr(dat, "active", active)
  dat <- set_attr(dat, "exitTime", exitTime)

  ## Summary statistics
  dat <- set_epi(dat, "d.flow", at, nDepts)

  return(dat)
}


# Updated Arrivals Module ----------------------------------------------------

afunc <- function(dat, at) {

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
    dat <- append_attr(dat, "age", 0, nArrivals)
  }

  ## Summary statistics
  dat <- set_epi(dat, "a.flow", at, nArrivals)

  return(dat)
}
