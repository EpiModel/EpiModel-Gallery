
##
## Test and Treat Intervention for an SIS Epidemic
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##
## Author: Samuel M. Jenness (Emory University)
## Date: August 2018
##


# Test and Treat Module ----------------------------------------------------

tnt <- function(dat, at) {
  # Simulate testing and diagnosis. Two attributes track separate aspects of
  # the cascade so that "tested recently" and "diagnosed positive" are not
  # conflated:
  #
  #   tested.status -- 1 during the test.dur window after any test
  #                    (positive or negative). Gates re-testing eligibility.
  #                    Resets to 0 after the window.
  #
  #   diag.status   -- 1 only for an infected individual whose test came
  #                    back positive. Gates the treatment-rate recovery in
  #                    the recovery module. Resets to 0 either when the
  #                    individual recovers or when their treatment course
  #                    (test.dur timesteps) expires without recovery.
  #
  # Universal screening: both susceptibles and infected individuals can
  # test, and only infected testers receive a positive diagnosis.

  ## Attributes ##
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")

  # EpiModel convention: modules first run at timestep 2 (timestep 1 is
  # reserved for initialization via init.net). Custom attributes that are
  # not set in init.net must be initialized here.
  if (at == 2) {
    dat <- set_attr(dat, "tested.status", rep(0, length(active)))
    dat <- set_attr(dat, "tested.time", rep(NA, length(active)))
    dat <- set_attr(dat, "diag.status", rep(0, length(active)))
    dat <- set_attr(dat, "diag.time", rep(NA, length(active)))
  }
  tested.status <- get_attr(dat, "tested.status")
  tested.time   <- get_attr(dat, "tested.time")
  diag.status   <- get_attr(dat, "diag.status")
  diag.time     <- get_attr(dat, "diag.time")

  ## Parameters ##
  test.rate <- get_param(dat, "test.rate")
  test.dur  <- get_param(dat, "test.dur")

  ## Testing process ##
  # Eligible: active, not currently in a test window, not currently in an
  # active treatment course.
  idsElig <- which(active == 1 & tested.status == 0 & diag.status == 0)
  nElig <- length(idsElig)

  # Stochastic testing: each eligible individual independently tests with
  # probability test.rate per timestep (Bernoulli trial).
  vecTest <- which(rbinom(nElig, 1, test.rate) == 1)
  idsTest <- idsElig[vecTest]
  nTest <- length(idsTest)

  # Every tester enters the test-window (gates re-testing).
  tested.status[idsTest] <- 1
  tested.time[idsTest]   <- at

  # Only infected testers receive a positive diagnosis (the rest are
  # tested-negative).
  idsTestPos <- idsTest[status[idsTest] == "i"]
  diag.status[idsTestPos] <- 1
  diag.time[idsTestPos]   <- at
  nDiagNew <- length(idsTestPos)

  ## Test-window reset ##
  # After test.dur timesteps, the tester is again eligible to test.
  idsTestedReset <- which(at - tested.time > (test.dur - 1))
  tested.status[idsTestedReset] <- 0
  tested.time[idsTestedReset]   <- NA

  ## Treatment-course reset ##
  # If a diagnosed individual has not recovered within test.dur timesteps,
  # the treatment course ends. They return to the natural clearance rate
  # until re-tested. (Recovery itself clears diag.status separately, in the
  # recovery module.)
  idsDiagReset <- which(at - diag.time > (test.dur - 1))
  diag.status[idsDiagReset] <- 0
  diag.time[idsDiagReset]   <- NA

  ## Write out updated attributes ##
  dat <- set_attr(dat, "tested.status", tested.status)
  dat <- set_attr(dat, "tested.time",   tested.time)
  dat <- set_attr(dat, "diag.status",   diag.status)
  dat <- set_attr(dat, "diag.time",     diag.time)

  ## Summary statistics ##
  dat <- set_epi(dat, "nTest", at, nTest)               # total tests (any result)
  dat <- set_epi(dat, "nDiagNew", at, nDiagNew)         # newly diagnosed (positive)
  dat <- set_epi(dat, "nRest", at, length(idsTestedReset))
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
