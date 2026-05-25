
##
## Social Diffusion Model
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##
## Authors: Samuel M. Jenness
## Date: September 2018
##

# Both modules implement **complex contagion** as an individual-level
# decision process: at each timestep, each non-adopter draws once for
# adoption based on the number of adopter contacts they currently have
# (their "exposure"). This contrasts with EpiModel's built-in SI module,
# which fires a Bernoulli trial per adopter-non-adopter EDGE per timestep
# (a per-contact / "simple contagion" hazard).
#
# Diffusion requires exposure: a non-adopter with zero adopter contacts
# has adoption probability zero in both modules. People cannot adopt a
# behavior they have never been exposed to through the network.


# Threshold diffusion module ------------------------------------------------
#
# Complex contagion with a hard threshold. A non-adopter only adopts when
# at least `min.degree` of their current contacts have already adopted.
# Above threshold the adoption probability is `adopt.prob`; below it the
# probability is 0. This captures phenomena like technology adoption
# ("I'll switch only when enough friends use it") or collective action
# ("I'll join only when enough peers are committed").

diffuse_threshold <- function(dat, at) {

  ## Attributes ##
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  infTime <- get_attr(dat, "infTime")

  ## Parameters ##
  adopt.prob <- get_param(dat, "adopt.prob")
  min.degree <- get_param(dat, "min.degree")

  ## Candidate set: all active non-adopters ##
  idsSus <- which(active == 1 & status == "s")
  nAdopt <- 0

  if (length(idsSus) > 0) {
    # Per-individual exposure count (0 for non-adopters with no adopter
    # contact). `tabulate` returns a length-N vector counting how often
    # each node id appears as the susceptible end of a discordant edge.
    del <- discord_edgelist(dat, at)
    if (is.null(del) || nrow(del) == 0) {
      exposure_all <- integer(length(active))
    } else {
      exposure_all <- tabulate(del$sus, nbins = length(active))
    }
    sus_exposure <- exposure_all[idsSus]

    # Per-individual adoption probability: adopt.prob if at/above threshold,
    # 0 otherwise (which includes exposure = 0).
    p <- ifelse(sus_exposure >= min.degree, adopt.prob, 0)

    # Single Bernoulli draw per non-adopter per timestep
    idsNewAdopt <- idsSus[rbinom(length(idsSus), 1, p) == 1]
    nAdopt <- length(idsNewAdopt)

    if (nAdopt > 0) {
      status[idsNewAdopt] <- "i"
      infTime[idsNewAdopt] <- at
      dat <- set_attr(dat, "status", status)
      dat <- set_attr(dat, "infTime", infTime)
    }
  }

  ## Summary statistics ##
  dat <- set_epi(dat, "si.flow", at, nAdopt)

  return(dat)
}


# Dose-response diffusion module --------------------------------------------
#
# Complex contagion with a smooth logistic dose-response. Adoption
# probability rises gradually with the number of adopter contacts:
#
#   P(adopt) = plogis(beta0 + beta1 * exposure)  for exposure >= 1
#   P(adopt) = 0                                  for exposure = 0
#
# Parameters:
#   beta0 -- intercept of the logistic curve (log-odds at exposure 1 minus
#            beta1). Typically negative so the probability at low exposure
#            is small.
#   beta1 -- slope (increase in log-odds per additional adopter contact).
#            Positive: more contacts produce higher adoption probability.
#
# This is a smooth generalization of the threshold model: rather than a
# hard cutoff, the curve smoothly rises with exposure. The exposure-zero
# case is still excluded because diffusion requires exposure.

diffuse_dose_response <- function(dat, at) {

  ## Attributes ##
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  infTime <- get_attr(dat, "infTime")

  ## Parameters ##
  beta0 <- get_param(dat, "beta0")
  beta1 <- get_param(dat, "beta1")

  ## Candidate set: all active non-adopters ##
  idsSus <- which(active == 1 & status == "s")
  nAdopt <- 0

  if (length(idsSus) > 0) {
    del <- discord_edgelist(dat, at)
    if (is.null(del) || nrow(del) == 0) {
      exposure_all <- integer(length(active))
    } else {
      exposure_all <- tabulate(del$sus, nbins = length(active))
    }
    sus_exposure <- exposure_all[idsSus]

    # Per-individual adoption probability via logistic dose-response.
    # Exposure-zero individuals are excluded: diffusion requires exposure.
    p <- ifelse(sus_exposure >= 1,
                plogis(beta0 + beta1 * sus_exposure),
                0)

    # Single Bernoulli draw per non-adopter per timestep
    idsNewAdopt <- idsSus[rbinom(length(idsSus), 1, p) == 1]
    nAdopt <- length(idsNewAdopt)

    if (nAdopt > 0) {
      status[idsNewAdopt] <- "i"
      infTime[idsNewAdopt] <- at
      dat <- set_attr(dat, "status", status)
      dat <- set_attr(dat, "infTime", infTime)
    }
  }

  ## Summary statistics ##
  dat <- set_epi(dat, "si.flow", at, nAdopt)

  return(dat)
}
