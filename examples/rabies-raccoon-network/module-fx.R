
##
## Rabies on an Observed Raccoon Contact Network
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##
## Author: Samuel M. Jenness (Emory University)
## Date: September 2026
##

# Disease states, stored in the `status` attribute:
#   "s"  susceptible
#   "e"  exposed: infected and incubating, not yet infectious
#   "i"  infectious: clinical rabies, transmits for one weekly step, then dies
#   "d"  dead of rabies
#   "r"  vaccinated (immune for the year)
#
# The network is an observed census with a fixed node set, so a raccoon that
# dies is not removed from the network: it stays as a node with status "d",
# which is neither susceptible nor infectious, so it takes no further part in
# transmission. A raccoon that was not collared in a given week has no
# contacts that week and so cannot transmit or be infected then.


# Vaccination at baseline ----------------------------------------------------

# Sets a random fraction `vax.cov` of the raccoons to "r" at the first time
# step, standing in for oral rabies vaccination reaching that share of the
# population before the introduction. Runs once.
vaccinate <- function(dat, at) {
  if (at > 2) {
    return(dat)
  }
  vax.cov <- get_param(dat, "vax.cov")
  status <- get_attr(dat, "status")
  n <- length(status)
  n.vax <- round(vax.cov * n)
  if (n.vax > 0) {
    idsVax <- sample(which(status == "s"), n.vax)
    status[idsVax] <- "r"
    dat <- set_attr(dat, "status", status)
  }
  dat <- set_epi(dat, "vax.num", at, n.vax)
  return(dat)
}


# Introduction of rabies ------------------------------------------------------

# At week `intro.week`, one susceptible raccoon that has at least one contact
# that week becomes infectious. Choosing among raccoons with a contact makes
# the introduction an animal that was actually in the collared population that
# week, rather than one whose collar had already failed or been removed.
introduce <- function(dat, at) {
  intro.week <- get_param(dat, "intro.week")
  if (at != intro.week) {
    return(dat)
  }
  status <- get_attr(dat, "status")
  infTime <- get_attr(dat, "infTime")
  el <- get_edgelist(dat, network = 1)
  present <- unique(c(el[, 1], el[, 2]))
  idsElig <- present[status[present] == "s"]
  if (length(idsElig) == 0) {
    idsElig <- which(status == "s")
  }
  idsIntro <- idsElig[sample.int(length(idsElig), 1)]
  status[idsIntro] <- "i"
  infTime[idsIntro] <- at
  dat <- set_attr(dat, "status", status)
  dat <- set_attr(dat, "infTime", infTime)
  # infectious since the previous step, so that it transmits during this
  # week's infection step and dies at the end of it, like every other case
  dat <- set_attr(dat, "iTime", at - 1, posit_ids = idsIntro)
  dat <- set_attr(dat, "introduced", TRUE, posit_ids = idsIntro)
  return(dat)
}


# Infection --------------------------------------------------------------------

# Transmission over the contacts active this week. The census layer gives the
# pairs in contact; the contact time of each pair this week comes from the
# `contact.hours` parameter, a list with one named vector per week (names are
# "raccoon1_raccoon2" with the smaller id first). The per-pair probability that
# an infectious raccoon infects a susceptible one is
#
#   1 - exp(-inf.hazard * hours)
#
# a constant hazard per hour of proximity, so a pair that dens together for
# the whole week is almost certain to transmit and a pair with a few minutes
# of contact rarely does. Every discordant pair is a separate Bernoulli trial;
# a susceptible exposed by more than one infectious raccoon in the same week
# is infected once, with the infector drawn at random among the successes.
infect <- function(dat, at) {
  status <- get_attr(dat, "status")
  infTime <- get_attr(dat, "infTime")
  inf.hazard <- get_param(dat, "inf.hazard")
  contact.hours <- get_param(dat, "contact.hours")

  nInf <- 0
  del <- discord_edgelist(dat, at, network = 1, infstat = "i")

  if (!is.null(del) && nrow(del) > 0) {
    key <- paste(pmin(del$sus, del$inf), pmax(del$sus, del$inf), sep = "_")
    hours <- contact.hours[[at]][key]
    hours[is.na(hours)] <- 0
    del$hours <- unname(hours)
    del$transProb <- 1 - exp(-inf.hazard * del$hours)
    transmit <- rbinom(nrow(del), 1, del$transProb)
    del <- del[which(transmit == 1), , drop = FALSE]

    if (nrow(del) > 0) {
      del <- del[sample.int(nrow(del)), , drop = FALSE]
      idsNewInf <- unique(del$sus)
      status[idsNewInf] <- "e"
      infTime[idsNewInf] <- at
      dat <- set_attr(dat, "status", status)
      dat <- set_attr(dat, "infTime", infTime)
      nInf <- length(idsNewInf)
      dat <- set_transmat(dat, del, at)
    }
  }

  dat <- set_epi(dat, "se.flow", at, nInf)
  return(dat)
}


# Progression ------------------------------------------------------------------

# Exposed raccoons become infectious at rate ei.rate per week (a geometric
# incubation period with mean 1 / ei.rate weeks). Infectious raccoons transmit
# for the week in which they are infectious and die at the end of it, which
# matches a clinical course of a few days to a week and a case fatality of
# one. Within a step, infection runs before progression, so a raccoon that
# becomes infectious in this step's progression transmits during next week's
# infection step and dies at the end of that week. `iTime` records the step
# in which a raccoon became infectious.
progress <- function(dat, at) {
  status <- get_attr(dat, "status")
  infTime <- get_attr(dat, "infTime")
  iTime <- get_attr(dat, "iTime", override.null.error = TRUE)
  if (is.null(iTime)) {
    iTime <- rep(NA, length(status))
  }
  ei.rate <- get_param(dat, "ei.rate")

  ## infectious raccoons die at the end of their infectious week
  idsDie <- which(status == "i" & !is.na(iTime) & iTime < at)
  status[idsDie] <- "d"

  ## exposed raccoons become infectious; they transmit from next week
  idsElig <- which(status == "e" & infTime < at)
  nEI <- 0
  if (length(idsElig) > 0) {
    vecEI <- which(rbinom(length(idsElig), 1, ei.rate) == 1)
    if (length(vecEI) > 0) {
      idsEI <- idsElig[vecEI]
      status[idsEI] <- "i"
      iTime[idsEI] <- at
      nEI <- length(idsEI)
    }
  }

  dat <- set_attr(dat, "status", status)
  dat <- set_attr(dat, "iTime", iTime)
  dat <- set_epi(dat, "ei.flow", at, nEI)
  dat <- set_epi(dat, "id.flow", at, length(idsDie))
  return(dat)
}


# Prevalence -------------------------------------------------------------------

prevalence <- function(dat, at) {
  status <- get_attr(dat, "status")
  dat <- set_epi(dat, "s.num", at, sum(status == "s"))
  dat <- set_epi(dat, "e.num", at, sum(status == "e"))
  dat <- set_epi(dat, "i.num", at, sum(status == "i"))
  dat <- set_epi(dat, "d.num", at, sum(status == "d"))
  dat <- set_epi(dat, "r.num", at, sum(status == "r"))
  dat <- set_epi(dat, "num", at, length(status))
  return(dat)
}
