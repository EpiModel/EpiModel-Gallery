
##
## Test and Treat Intervention for an SIS Epidemic
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##
## Author: Samuel M. Jenness (Emory University)
## Date: August 2018
##


# Test and Treat Module ----------------------------------------------------

tnt <- function(dat, at) {
  # Simulate testing and diagnosis. Undiagnosed individuals test at rate
  # test.rate per timestep. Diagnosis persists for test.dur timesteps
  # (representing the treatment course duration), then resets.

  ## Attributes ##
  active <- get_attr(dat, "active")

  # EpiModel convention: modules first run at timestep 2 (timestep 1 is
  # reserved for initialization via init.net). Custom attributes that are
  # not set in init.net must be initialized here.
  if (at == 2) {
    dat <- set_attr(dat, "diag.status", rep(0, length(active)))
    dat <- set_attr(dat, "diag.time", rep(NA, length(active)))
  }
  diag.status <- get_attr(dat, "diag.status")
  diag.time <- get_attr(dat, "diag.time")

  ## Parameters ##
  test.rate <- get_param(dat, "test.rate")
  test.dur <- get_param(dat, "test.dur")

  ## Testing process ##
  # Eligible: active and not currently diagnosed. Both susceptible and
  # infected individuals test (testing does not require symptoms).
  idsElig <- which(active == 1 & diag.status == 0)
  nElig <- length(idsElig)

  # Stochastic testing: each eligible individual independently tests
  # with probability test.rate per timestep (Bernoulli trial)
  vecTest <- which(rbinom(nElig, 1, test.rate) == 1)
  idsTest <- idsElig[vecTest]
  nTest <- length(idsTest)

  ## Update diagnosis status for newly tested ##
  diag.status[idsTest] <- 1
  diag.time[idsTest] <- at

  ## Diagnosis reset ##
  # Diagnosis persists for exactly test.dur timesteps. After that, the
  # diagnosis status resets to 0, representing the end of the treatment
  # course. If the individual is still infected, they return to the
  # untreated recovery rate until re-tested.
  idsReset <- which(at - diag.time > (test.dur - 1))
  diag.status[idsReset] <- 0
  diag.time[idsReset] <- NA

  ## Write out updated attributes ##
  dat <- set_attr(dat, "diag.status", diag.status)
  dat <- set_attr(dat, "diag.time", diag.time)

  ## Summary statistics ##
  dat <- set_epi(dat, "nTest", at, nTest)
  dat <- set_epi(dat, "nRest", at, length(idsReset))
  dat <- set_epi(dat, "nDiag", at, sum(diag.status == 1, na.rm = TRUE))

  return(dat)
}


# Recovery Module -----------------------------------------------------------

recov <- function(dat, at) {
  # Simulate recovery (I -> S) with diagnosis-dependent treatment rates.
  # Diagnosed individuals recover at rate rec.rate.tx (treatment), while
  # undiagnosed individuals recover at the slower rate rec.rate (natural
  # clearance). On recovery, diagnosis status is also cleared.

  ## Attributes ##
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  diag.status <- get_attr(dat, "diag.status")

  ## Parameters ##
  rec.rate <- get_param(dat, "rec.rate")
  rec.rate.tx <- get_param(dat, "rec.rate.tx")

  ## Determine eligible to recover ##
  idsElig <- which(active == 1 & status == "i")
  nElig <- length(idsElig)

  ## Heterogeneous recovery rates based on diagnosis status ##
  # This is the core mechanism of the test-and-treat intervention:
  # diagnosed individuals receive treatment and recover faster.
  diag.status.elig <- diag.status[idsElig]
  ratesElig <- ifelse(diag.status.elig == 1, rec.rate.tx, rec.rate)

  ## Stochastic recovery process ##
  vecRecov <- which(rbinom(nElig, 1, ratesElig) == 1)
  idsRecov <- idsElig[vecRecov]
  nRecov <- length(idsRecov)

  ## Update attributes for recovered individuals ##
  # Recovery clears both disease status and diagnosis status
  status[idsRecov] <- "s"
  diag.status[idsRecov] <- 0

  ## Write out updated attributes ##
  dat <- set_attr(dat, "status", status)
  dat <- set_attr(dat, "diag.status", diag.status)

  ## Summary statistics ##
  dat <- set_epi(dat, "is.flow", at, nRecov)

  return(dat)
}
