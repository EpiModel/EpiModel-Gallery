
##
## Modeling Epidemics over Observed Networks
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##
## Author: Samuel M. Jenness (Emory University)
## Date: August 2018
##


# Initialization Module ---------------------------------------------------

# Custom initialization for observed (census) network data. The standard
# initialize.net expects a netest object (fitted ERGM), but observed networks
# have no model fit -- the network is placed directly into the dat object.
#
# When track.nw.attr = TRUE (set in param.net), time-varying disease status
# is stored directly on the networkDynamic object as a temporally extended
# attribute (TEA) named "testatus", enabling network visualization with
# plot(sim, type = "network", col.status = TRUE). This requires the low-level
# networkDynamic::activate.vertex.attribute() API because EpiModel's set_attr()
# does not support time-varying network attributes.

init_obsnw <- function(x, param, init, control, s) {

  # Master data list
  dat <- create_dat_object(param, init, control)

  # Place observed network directly (no ERGM fit to simulate from)
  dat$num.nw <- 1L
  dat$run$nw[[1]] <- x
  dat <- set_param(dat, "groups", 1)

  # Core attributes and infection status
  i.num <- get_init(dat, "i.num")
  n <- network.size(dat$run$nw[[1]])
  dat <- append_core_attr(dat, 1, n)

  status <- rep("s", n)
  status[sample(1:n, i.num)] <- "i"
  dat <- set_attr(dat, "status", status)

  infTime <- rep(NA, n)
  infTime[which(status == "i")] <- 1
  dat <- set_attr(dat, "infTime", infTime)

  # Optional: set time-varying status on network for visualization
  track.nw.attr <- get_param(dat, "track.nw.attr", override.null.error = TRUE)
  if (!is.null(track.nw.attr) && track.nw.attr) {
    dat$run$nw[[1]] <- networkDynamic::activate.vertex.attribute(
      dat$run$nw[[1]], prefix = "testatus",
      value = get_attr(dat, "status"), onset = 1, terminus = Inf
    )
  }

  dat <- prevalence.net(dat, 1)
  return(dat)
}


# Infection Module ---------------------------------------------------------

# Transmission over the observed network. Supports two modes:
#
#   Simple mode (inf.prob set): constant per-act transmission probability
#
#   Time-varying mode (inf.prob.stage1 set): transmission probability depends
#     on infection duration. Primary stage (duration <= dur.stage1) uses
#     inf.prob.stage1; secondary stage uses inf.prob.stage2.

infect_obsnw <- function(dat, at) {

  ## Attributes ##
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  infTime <- get_attr(dat, "infTime")

  ## Parameters ##
  act.rate <- get_param(dat, "act.rate")
  inf.prob.stage1 <- get_param(dat, "inf.prob.stage1",
                                override.null.error = TRUE)

  # Determine transmission mode
  time_varying <- !is.null(inf.prob.stage1)
  if (time_varying) {
    inf.prob.stage2 <- get_param(dat, "inf.prob.stage2")
    dur.stage1 <- get_param(dat, "dur.stage1")
  } else {
    inf.prob <- get_param(dat, "inf.prob")
  }

  # Vector of infected and susceptible IDs
  idsInf <- which(active == 1 & status == "i")
  nActive <- sum(active == 1)
  nElig <- length(idsInf)

  # Initialize vectors at 0
  totInf <- 0

  ## Processes ##
  # If some infected AND some susceptible, then proceed
  if (nElig > 0 && nElig < nActive) {

    # Get discordant edgelist
    del <- discord_edgelist(dat, at)

    # If some discordant edges, then proceed
    if (!is.null(del)) {

      # Infection probabilities
      if (time_varying) {
        infDur.del <- at - infTime[del$inf]
        del$transProb <- ifelse(infDur.del <= dur.stage1,
                                inf.prob.stage1, inf.prob.stage2)
      } else {
        del$transProb <- inf.prob
      }

      # Act rates
      del$actRate <- act.rate

      # Calculate final transmission probability per timestep
      del$finalProb <- 1 - (1 - del$transProb) ^ del$actRate

      # Randomize transmissions and subset df
      transmit <- rbinom(nrow(del), 1, del$finalProb)
      del <- del[which(transmit == 1), ]

      # Set new infections vector
      idsNewInf <- unique(del$sus)
      totInf <- length(idsNewInf)

      # Update attributes for newly infected
      if (totInf > 0) {
        status[idsNewInf] <- "i"
        infTime[idsNewInf] <- at
        dat <- set_attr(dat, "status", status)
        dat <- set_attr(dat, "infTime", infTime)

        # Update time-varying status on network for visualization
        track.nw.attr <- get_param(dat, "track.nw.attr",
                                   override.null.error = TRUE)
        if (!is.null(track.nw.attr) && track.nw.attr) {
          dat$run$nw[[1]] <- networkDynamic::activate.vertex.attribute(
            dat$run$nw[[1]], prefix = "testatus",
            value = "i", onset = at, terminus = Inf, v = idsNewInf
          )
        }
      }

    }
  }

  ## Summary statistics ##
  dat <- set_epi(dat, "si.flow", at, totInf)

  return(dat)
}
