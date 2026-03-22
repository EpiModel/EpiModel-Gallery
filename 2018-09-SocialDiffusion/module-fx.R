
##
## Social Diffusion Model
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##
## Authors: Samuel M. Jenness
## Date: November 2018
##


# Threshold diffusion module ------------------------------------------------
#
# Models "complex contagion" where adoption requires social reinforcement.
# A susceptible individual only adopts when they have at least `min.degree`
# contacts who have already adopted. Below this threshold, the adoption
# probability is 0 regardless of other parameters.
#
# This captures phenomena like technology adoption (need to see several
# friends using it), behavior change (need multiple role models), or protest
# participation (need a critical mass of committed peers).

diffuse_threshold <- function(dat, at) {

  ## Attributes ##
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  infTime <- get_attr(dat, "infTime")

  ## Parameters ##
  inf.prob <- get_param(dat, "inf.prob")
  act.rate <- get_param(dat, "act.rate")
  min.degree <- get_param(dat, "min.degree")

  ## Find adopted ("infected") nodes ##
  idsInf <- which(active == 1 & status == "i")
  nActive <- sum(active == 1)
  nElig <- length(idsInf)

  ## Initialize adoption count ##
  nAdopt <- 0

  ## Diffusion process ##
  if (nElig > 0 && nElig < nActive) {

    # Get discordant edgelist: edges between adopters ("i") and
    # non-adopters ("s"). Each row represents one such edge.
    del <- discord_edgelist(dat, at)

    if (!is.null(del)) {

      # Count how many adopter contacts each non-adopter has.
      # This is the "exposure count": the number of current partners who
      # have already adopted. A susceptible node connected to 3 adopters
      # has 3 rows in the DEL, so aggregate gives exposure = 3.
      exposure <- aggregate(list(exposure = rep(1, nrow(del))),
                            by = list(sus = del$sus), FUN = sum)
      del <- merge(del, exposure, by = "sus")

      # Threshold rule: adoption probability is nonzero only when the
      # exposure count meets or exceeds the minimum threshold.
      # Below threshold: no chance of adoption, no matter how many acts.
      # At/above threshold: adoption probability = inf.prob per act.
      del$transProb <- ifelse(del$exposure >= min.degree, inf.prob, 0)

      del$actRate <- act.rate
      del$finalProb <- 1 - (1 - del$transProb)^del$actRate

      # Stochastic adoption process
      transmit <- rbinom(nrow(del), 1, del$finalProb)
      del <- del[which(transmit == 1), ]

      idsNewAdopt <- unique(del$sus)
      nAdopt <- length(idsNewAdopt)

      if (nAdopt > 0) {
        status[idsNewAdopt] <- "i"
        infTime[idsNewAdopt] <- at
        dat <- set_attr(dat, "status", status)
        dat <- set_attr(dat, "infTime", infTime)
      }
    }
  }

  ## Summary statistics ##
  dat <- set_epi(dat, "si.flow", at, nAdopt)

  return(dat)
}


# Dose-response diffusion module --------------------------------------------
#
# Models adoption probability as a continuous logistic function of the
# number of adopter contacts ("exposure count"):
#
#   P(adopt per act) = plogis(beta0 + beta1 * exposure)
#                    = 1 / (1 + exp(-(beta0 + beta1 * exposure)))
#
# Parameters:
#   beta0: intercept (log-odds of adoption with 1 adopter contact, minus beta1).
#          Typically negative so baseline probability is low.
#   beta1: slope (increase in log-odds per additional adopter contact).
#          Positive values mean more contacts -> higher adoption probability.
#
# This is a smooth generalization of the threshold model. Instead of a hard
# cutoff, adoption probability rises gradually with exposure. It produces
# intermediate dynamics between simple contagion (constant probability) and
# threshold contagion (all-or-nothing).

diffuse_dose_response <- function(dat, at) {

  ## Attributes ##
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  infTime <- get_attr(dat, "infTime")

  ## Parameters ##
  beta0 <- get_param(dat, "beta0")
  beta1 <- get_param(dat, "beta1")
  act.rate <- get_param(dat, "act.rate")

  ## Find adopted ("infected") nodes ##
  idsInf <- which(active == 1 & status == "i")
  nActive <- sum(active == 1)
  nElig <- length(idsInf)

  ## Initialize adoption count ##
  nAdopt <- 0

  ## Diffusion process ##
  if (nElig > 0 && nElig < nActive) {

    del <- discord_edgelist(dat, at)

    if (!is.null(del)) {

      # Count adopter contacts per non-adopter (same as threshold module)
      exposure <- aggregate(list(exposure = rep(1, nrow(del))),
                            by = list(sus = del$sus), FUN = sum)
      del <- merge(del, exposure, by = "sus")

      # Logistic dose-response: adoption probability is a smooth function
      # of exposure count, parameterized by beta0 (intercept) and beta1 (slope)
      del$transProb <- plogis(beta0 + beta1 * del$exposure)

      del$actRate <- act.rate
      del$finalProb <- 1 - (1 - del$transProb)^del$actRate

      # Stochastic adoption process
      transmit <- rbinom(nrow(del), 1, del$finalProb)
      del <- del[which(transmit == 1), ]

      idsNewAdopt <- unique(del$sus)
      nAdopt <- length(idsNewAdopt)

      if (nAdopt > 0) {
        status[idsNewAdopt] <- "i"
        infTime[idsNewAdopt] <- at
        dat <- set_attr(dat, "status", status)
        dat <- set_attr(dat, "infTime", infTime)
      }
    }
  }

  ## Summary statistics ##
  dat <- set_epi(dat, "si.flow", at, nAdopt)

  return(dat)
}
