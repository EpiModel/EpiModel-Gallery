
##
## Test and Treat Intervention for an SIS Epidemic
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##
## Author: Samuel M. Jenness (Emory University)
## Date: August 2018
##

# New Test and Treat Module -----------------------------------------------

tnt <- function(dat, at) {

  ## Uncomment this to function environment interactively
  # browser()

  ## Attributes ##
  active <- get_attr(dat, "active")

  # Initialize new Dx attr at sim start
  if (at == 2) {
    dat <- set_attr(dat, "diag.status", rep(0, length(active)))
    dat <- set_attr(dat, "diag.time", rep(NA, length(active)))
  }
  diag.status <- get_attr(dat, "diag.status")
  diag.time <- get_attr(dat, "diag.time")

  ## Parameters ##
  test.rate <- get_param(dat, "test.rate")
  test.dur <- get_param(dat,"test.dur")

  ## Determine eligible to test ##
  idsElig <- which(active == 1 & diag.status == 0)
  nElig <- length(idsElig)

  ## Vector of recovered IDs after stochastic process ##
  vecTest <- which(rbinom(nElig, 1, test.rate) == 1)
  idsTest <- idsElig[vecTest]
  nTest <- length(idsTest)

  ## Update attribute if any tested ##
  diag.status[idsTest] <- 1
  diag.time[idsTest] <- at

  ## Dx lasts for test.dur weeks, then is reset ##
  ## Probably should be related to 1/rec.rate.tx ##
  idsReset <- which(at - diag.time > (test.dur - 1))
  diag.status[idsReset] <- 0
  diag.time[idsReset] <- NA

  ## Write out updated attributes ##
  dat <- set_attr(dat, "diag.status", diag.status)
  dat <- set_attr(dat, "diag.time", diag.time)

  ## Write out summary statistics ##
  dat <- set_epi(dat, "nTest", at, nTest)
  dat <- set_epi(dat, "nRest", at, length(idsReset))

  return(dat)
}


# Updated Recovery Module --------------------------------------------------

recov <- function(dat, at) {

  ## Uncomment this to function environment interactively
  # browser()

  ## Attributes ##
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  diag.status <- get_attr(dat, "diag.status")
  diag.time <- get_attr(dat, "diag.time")

  ## Parameters ##
  rec.rate <- get_param(dat, "rec.rate")
  rec.rate.tx <- get_param(dat, "rec.rate.tx")

  ## Determine eligible to recover ##
  idsElig <- which(active == 1 & status == "i")
  nElig <- length(idsElig)

  ## Determined Dx status of eligible ##
  diag.status.elig <- diag.status[idsElig]

  ## Recovery rates dependent on Dx status ##
  ratesElig <- ifelse(diag.status.elig == 1, rec.rate.tx, rec.rate)

  ## Vector of recovered IDs after stochastic process
  vecRecov <- which(rbinom(nElig, 1, ratesElig) == 1)
  idsRecov <- idsElig[vecRecov]
  nRecov <- length(idsRecov)

  ## Update attributes if any recovered ##
  status[idsRecov] <- "s"
  diag.status[idsRecov] <- 0
  diag.time[idsRecov] <- NA

  ## Write out updated attributes ##
  dat <- set_attr(dat, "status", status)
  dat <- set_attr(dat, "diag.status", diag.status)
  dat <- set_attr(dat, "diag.time", diag.time)

  ## Write out summary statistics ##
  dat <- set_epi(dat, "is.flow", at, nRecov)

  return(dat)
}
