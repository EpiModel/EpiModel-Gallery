##
## Rabies on an Observed Raccoon Contact Network
## EpiModel Gallery (https://github.com/EpiModel/EpiModel-Gallery)
##
## Author: Samuel M. Jenness (Emory University)
## Date: September 2026
##

# Disease states, stored in the `status` attribute:
#   "s"  susceptible
#   "e"  exposed: infected and incubating, not yet infectious
#   "i"  infectious: clinical rabies, transmits for one weekly step, then dies
#   "d"  dead of rabies
#   "r"  vaccinated (immune for the rest of the run)
#
# The node set of an observed network is fixed, so a raccoon that dies is not
# removed from the network: it stays as a node with status "d", which is
# neither susceptible nor infectious, so it takes no further part in
# transmission. A raccoon with no recorded contact in a given week (most
# often because its collar was not working) cannot transmit or be infected
# that week.


# Vaccination at baseline ----------------------------------------------------

# Sets round(vax.cov * n) randomly chosen raccoons to "r" at step 2, the
# first step at which modules run, standing in for oral rabies vaccination
# that left that share of the population immune before any introduction.
# Runs once.
vaccinate <- function(dat, at) {
  if (at > 2) {
    return(dat)
  }
  vax.cov <- get_param(dat, "vax.cov")
  status <- get_attr(dat, "status")
  n.vax <- round(vax.cov * length(status))
  if (n.vax > 0) {
    idsVax <- sample(which(status == "s"), n.vax)
    status[idsVax] <- "r"
    dat <- set_attr(dat, "status", status)
  }
  return(dat)
}


# Introduction of rabies ------------------------------------------------------

# One infectious raccoon is introduced in week `intro.week`. When `intro.week`
# is a vector of candidate weeks, each simulation draws one of them at its
# first step and stores it with set_param(), so a single netsim() call covers
# every introduction week. The week is recorded in the epi output.
#
# The introduced raccoon is drawn from the susceptible raccoons with at least
# one contact that week, so it is an animal that was in the recorded network
# at the time. If every raccoon in contact that week is immune, it is drawn
# from the rest, and the introduction goes nowhere.
introduce <- function(dat, at) {
  intro.week <- get_param(dat, "intro.week")
  if (at == 2 && length(intro.week) > 1) {
    intro.week <- intro.week[sample.int(length(intro.week), 1)]
    dat <- set_param(dat, "intro.week", intro.week)
  }
  dat <- set_epi(dat, "intro.week", intro.week)
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
  return(dat)
}


# Infection --------------------------------------------------------------------

# Transmission over the contacts active this week. The census layer gives the
# pairs in contact; the contact time of each pair this week comes from the
# `contact.hours` parameter, a list with one named vector per week (names are
# "id1_id2" with the smaller id first). The per-pair probability that an
# infectious raccoon infects a susceptible one is
#
#   1 - exp(-inf.hazard * hours)
#
# a constant hazard per hour of proximity. This is an assumption: rabies is
# transmitted by bites, and proximity time stands in for the opportunity to
# bite. Every discordant pair is a separate Bernoulli trial; a susceptible
# exposed by more than one infectious raccoon in the same week is infected
# once, with the infector drawn at random among the successes.
#
# Each new infection draws its incubation period, the number of weeks from
# infection to the week in which it is infectious, from a gamma distribution
# with mean inc.mean and shape inc.shape, rounded to whole weeks (at least
# one).
infect <- function(dat, at) {
  status <- get_attr(dat, "status")
  infTime <- get_attr(dat, "infTime")
  incub <- get_attr(dat, "incub", override.null.error = TRUE)
  if (is.null(incub)) {
    incub <- rep(NA_real_, length(status))
  }
  inf.hazard <- get_param(dat, "inf.hazard")
  inc.mean <- get_param(dat, "inc.mean")
  inc.shape <- get_param(dat, "inc.shape")
  contact.hours <- get_param(dat, "contact.hours")

  nInf <- 0
  del <- discord_edgelist(dat, at, network = 1, infstat = "i")

  if (!is.null(del) && nrow(del) > 0) {
    key <- paste(pmin(del$sus, del$inf), pmax(del$sus, del$inf), sep = "_")
    hours <- contact.hours[[at]][key]
    if (anyNA(hours)) {
      stop("census edges without contact hours at step ", at,
           ": the network and contact.hours disagree")
    }
    del$hours <- unname(hours)
    del$transProb <- 1 - exp(-inf.hazard * del$hours)
    transmit <- rbinom(nrow(del), 1, del$transProb)
    del <- del[which(transmit == 1), , drop = FALSE]

    if (nrow(del) > 0) {
      # shuffle, then keep one successful exposure per newly infected raccoon
      del <- del[sample.int(nrow(del)), , drop = FALSE]
      del <- del[!duplicated(del$sus), , drop = FALSE]
      idsNewInf <- del$sus
      status[idsNewInf] <- "e"
      infTime[idsNewInf] <- at
      incub[idsNewInf] <- pmax(1, round(rgamma(length(idsNewInf), shape = inc.shape,
                                               scale = inc.mean / inc.shape)))
      dat <- set_attr(dat, "status", status)
      dat <- set_attr(dat, "infTime", infTime)
      dat <- set_attr(dat, "incub", incub)
      nInf <- length(idsNewInf)
      dat <- set_transmat(dat, del, at)
    }
  }

  dat <- set_epi(dat, "se.flow", nInf)
  return(dat)
}


# Progression ------------------------------------------------------------------

# A raccoon infected in week t with incubation period d is infectious in week
# t + d: it becomes "i" in the progression step of week t + d - 1, after that
# week's infection step, transmits in the infection step of week t + d, and
# dies at the end of that week. One week of clinical disease followed by
# death matches a clinical course of a few days to a week and a case
# fatality of one. `iTime` records the step in which a raccoon became
# infectious.
progress <- function(dat, at) {
  status <- get_attr(dat, "status")
  infTime <- get_attr(dat, "infTime")
  incub <- get_attr(dat, "incub", override.null.error = TRUE)
  iTime <- get_attr(dat, "iTime", override.null.error = TRUE)
  if (is.null(iTime)) {
    iTime <- rep(NA, length(status))
  }

  ## infectious raccoons die at the end of their infectious week
  idsDie <- which(status == "i" & !is.na(iTime) & iTime < at)
  status[idsDie] <- "d"

  ## exposed raccoons whose infectious week is next week become infectious
  idsEI <- which(status == "e" & at >= infTime + incub - 1)
  status[idsEI] <- "i"
  iTime[idsEI] <- at

  dat <- set_attr(dat, "status", status)
  dat <- set_attr(dat, "iTime", iTime)
  dat <- set_epi(dat, "ei.flow", length(idsEI))
  dat <- set_epi(dat, "id.flow", length(idsDie))
  return(dat)
}


# Prevalence -------------------------------------------------------------------

prevalence <- function(dat, at) {
  status <- get_attr(dat, "status")
  dat <- set_epi(dat, "s.num", sum(status == "s"))
  dat <- set_epi(dat, "e.num", sum(status == "e"))
  dat <- set_epi(dat, "i.num", sum(status == "i"))
  dat <- set_epi(dat, "d.num", sum(status == "d"))
  dat <- set_epi(dat, "r.num", sum(status == "r"))
  dat <- set_epi(dat, "num", length(status))
  return(dat)
}
