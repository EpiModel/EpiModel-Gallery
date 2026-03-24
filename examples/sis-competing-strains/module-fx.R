
##
## Modeling Two Competing Strains in an SIS Epidemic
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##
## Author: Steven M. Goodreau (University of Washington)
## Date: September 2018
##


# Strain Initialization Module ---------------------------------------------

init_strain <- function(dat, at) {
  # Initialize the strain attribute for all nodes. Runs once at the first

  # module call (timestep 2). Infected nodes are randomly assigned strain 1
  # or strain 2, with strain 2 probability controlled by pct.st2.
  # Susceptible nodes receive NA (no strain).

  # Only run once: check if the strain attribute already exists
  if (is.null(get_attr(dat, "strain", override.null.error = TRUE))) {

    ## Attributes ##
    active <- get_attr(dat, "active")
    status <- get_attr(dat, "status")

    nActive <- sum(active == 1)
    idsInf <- which(active == 1 & status == "i")

    ## Parameters ##
    pct.st2 <- get_param(dat, "pct.st2")

    ## Assign strains ##
    # Default is NA (susceptible nodes have no strain)
    strain <- rep(NA, nActive)

    # Infected nodes: Bernoulli draw for strain 2, then +1 to get labels 1 or 2
    strain[idsInf] <- rbinom(length(idsInf), 1, pct.st2) + 1

    ## Write out attribute ##
    dat <- set_attr(dat, "strain", strain)
  }

  return(dat)
}


# Two-Strain Infection Module ----------------------------------------------

infection.2strains <- function(dat, at) {
  # Simulate transmission of two competing strains along discordant edges.
  # Each strain has its own per-act transmission probability (inf.prob for
  # strain 1, inf.prob.st2 for strain 2). Both can be time-varying vectors
  # indexed by infection duration. Newly infected individuals inherit the
  # strain of their infecting partner.
  #
  # When a susceptible node has discordant edges with multiple infected

  # partners and multiple transmissions occur, only the first (by row order
  # in the discordant edgelist) determines the strain assignment.

  ## Attributes ##
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  infTime <- get_attr(dat, "infTime")
  strain <- get_attr(dat, "strain")

  ## Parameters ##
  inf.prob <- get_param(dat, "inf.prob")
  act.rate <- get_param(dat, "act.rate")
  inf.prob.st2 <- get_param(dat, "inf.prob.st2")

  ## Eligible nodes ##
  idsInf <- which(active == 1 & status == "i")
  nActive <- sum(active == 1)
  nElig <- length(idsInf)

  ## Initialize flow counters ##
  nInf <- nInfST2 <- totInf <- 0

  ## Transmission process ##
  # Requires at least one infected AND at least one susceptible node
  if (nElig > 0 && nElig < nActive) {

    # Get discordant edgelist (edges between S and I nodes)
    del <- discord_edgelist(dat, at)

    if (!(is.null(del))) {

      # Infection duration for each infected partner
      del$infDur <- at - infTime[del$inf]
      del$infDur[del$infDur == 0] <- 1

      # Strain-specific transmission probabilities.
      # Both inf.prob and inf.prob.st2 can be vectors indexed by infection
      # duration (allowing time-varying infectiousness). If infDur exceeds
      # the vector length, the last element is used.
      linf.prob <- length(inf.prob)
      linf.prob.st2 <- length(inf.prob.st2)

      del$transProb <- ifelse(strain[del$inf] == 1,
                              ifelse(del$infDur <= linf.prob,
                                     inf.prob[del$infDur],
                                     inf.prob[linf.prob]),
                              ifelse(del$infDur <= linf.prob.st2,
                                     inf.prob.st2[del$infDur],
                                     inf.prob.st2[linf.prob.st2]))

      # Optional intervention: reduces transmission by (1 - inter.eff)
      # starting at timestep inter.start
      if (!is.null(dat$param$inter.eff) && at >= dat$param$inter.start) {
        del$transProb <- del$transProb * (1 - dat$param$inter.eff)
      }

      # Act rate (can also be a time-varying vector indexed by infDur)
      lact.rate <- length(act.rate)
      del$actRate <- ifelse(del$infDur <= lact.rate,
                            act.rate[del$infDur],
                            act.rate[lact.rate])

      # Per-timestep transmission probability:
      # P(transmit) = 1 - (1 - transProb)^actRate
      del$finalProb <- 1 - (1 - del$transProb) ^ del$actRate

      # Stochastic transmission (Bernoulli trial per discordant edge)
      transmit <- rbinom(nrow(del), 1, del$finalProb)
      del <- del[which(transmit == 1), ]

      # Identify newly infected nodes (unique susceptibles)
      idsNewInf <- unique(del$sus)

      # Resolve double infections: if a susceptible is infected by multiple
      # partners in the same timestep, take the strain of the first infector
      # (by row order in the discordant edgelist)
      infpairs <- sapply(idsNewInf, function(x) min(which(del$sus == x)))

      if (length(infpairs) > 0) {
        infectors <- del$inf[infpairs]
        strain <- get_attr(dat, "strain")
        infectors_strain <- strain[infectors]

        nInf <- sum(infectors_strain == 1)
        nInfST2 <- sum(infectors_strain == 2)
        totInf <- nInf + nInfST2
      } else {
        nInf <- nInfST2 <- totInf <- 0
      }

      ## Update attributes for newly infected ##
      if (totInf > 0) {
        status[idsNewInf] <- "i"
        dat <- set_attr(dat, "status", status)
        infTime[idsNewInf] <- at
        dat <- set_attr(dat, "infTime", infTime)
        strain[idsNewInf] <- strain[infectors]
        dat <- set_attr(dat, "strain", strain)
      }

    }
  }

  ## Summary statistics ##

  # Incidence by strain
  dat <- set_epi(dat, "si.flow", at, nInf)
  dat <- set_epi(dat, "si.flow.st2", at, nInfST2)

  # Prevalence by strain
  dat <- set_epi(dat, "i.num.st1", at, sum(status == "i" & strain == 1))
  dat <- set_epi(dat, "i.num.st2", at, sum(status == "i" & strain == 2))

  return(dat)
}


# Two-Strain Recovery Module -----------------------------------------------

recov.2strains <- function(dat, at) {
  # Simulate recovery (I -> S) with strain-dependent recovery rates.
  # Strain 1 recovers at rec.rate, strain 2 at rec.rate.st2.
  # On recovery, both disease status and strain are cleared.

  ## Attributes ##
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  strain <- get_attr(dat, "strain")

  ## Parameters ##
  rec.rate <- get_param(dat, "rec.rate")
  rec.rate.st2 <- get_param(dat, "rec.rate.st2")

  ## Determine eligible to recover ##
  idsElig <- which(active == 1 & status == "i")
  nElig <- length(idsElig)

  ## Strain-dependent recovery rates ##
  strain.elig <- strain[idsElig]
  ratesElig <- ifelse(strain.elig == 1, rec.rate, rec.rate.st2)

  ## Stochastic recovery (Bernoulli trial per infected node) ##
  vecRecov <- which(rbinom(nElig, 1, ratesElig) == 1)
  idsRecov <- idsElig[vecRecov]
  nRecov <- length(idsRecov)
  nRecov.st1 <- sum(strain[idsRecov] == 1)
  nRecov.st2 <- sum(strain[idsRecov] == 2)

  ## Update attributes for recovered individuals ##
  # Recovery clears both disease status and strain assignment
  status[idsRecov] <- "s"
  strain[idsRecov] <- NA

  ## Write out updated attributes ##
  dat <- set_attr(dat, "status", status)
  dat <- set_attr(dat, "strain", strain)

  ## Summary statistics ##
  dat <- set_epi(dat, "is.flow", at, nRecov)
  dat <- set_epi(dat, "is.flow.st1", at, nRecov.st1)
  dat <- set_epi(dat, "is.flow.st2", at, nRecov.st2)

  return(dat)
}
