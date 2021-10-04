
##
## Simple Cost-effectiveness Model (Simple SI model with cost and utility tracking)
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##
## Author: Gregory Knowlton (University of Minnesota)
## Date: October 2021
##

# Cost/Utility Tracking Module ------------------------------------------------

costeffect <- function(dat, at) {
  
  # If the start of the analytic time horizon has not been reached, skip module
  cea.start <- get_param(dat, "cea.start")
  if (at < cea.start) {
    return(dat)
  }
  
  # Import attributes
  active <- get_attr(dat, "active")
  age <- get_attr(dat, "age")
  status <- get_attr(dat, "status")
  
  # Import parameters
  sus.cost <- get_param(dat, "sus.cost")
  inf.cost <- get_param(dat, "inf.cost")
  sus.qaly <- get_param(dat, "sus.qaly")
  inf.qaly <- get_param(dat, "inf.qaly")
  age.decrement <- get_param(dat, "age.decrement")
  disc.rate <- get_param(dat, "disc.rate")
  inter.cost <- get_param(dat, "inter.cost", override.null.error = TRUE)
  inter.start <- get_param(dat, "inter.start", override.null.error = TRUE)
  end.horizon <- get_param(dat, "end.horizon", override.null.error = TRUE)
  
  # Identify relevant sub-populations
  idsSus <- which(active == 1 & status == "s")
  idsInf <- which(active == 1 & status == "i")
  
  # Account for costs of each sub-population
  pop.sus.cost <- length(idsSus) * sus.cost
  pop.inf.cost <- length(idsInf) * inf.cost
  
  # Account for intervention costs
  # Costs accrue uniformly from intervention start to end horizon
  if(!is.null(inter.start) & at < end.horizon) {
    inter.cost <- inter.cost / (end.horizon - inter.start)
  } else {
    inter.cost <- 0
  }
  
  # Account for effects (QALYs) of each sub-population
  # QALY parameters by health state and age decrements must be converted to weekly time-step
  pop.sus.qaly <- sum((age[idsSus] * age.decrement + sus.qaly) / 52, na.rm = TRUE)
  pop.inf.qaly <- sum((age[idsInf] * age.decrement + inf.qaly) / 52, na.rm = TRUE)
  
  # Aggregate costs and effects
  pop.cost <- pop.sus.cost + pop.inf.cost + inter.cost
  pop.qaly <- pop.sus.qaly + pop.inf.qaly
  
  # Discount aggregated costs and effects
  t <- at - cea.start
  pop.cost.disc <- pop.cost * (1 - disc.rate) ^ (t / 52)
  pop.qaly.disc <- pop.qaly * (1 - disc.rate) ^ (t / 52)
  
  ## Summary statistics ##
  dat <- set_epi(dat, "cost", at, pop.cost)
  dat <- set_epi(dat, "qaly", at, pop.qaly)
  dat <- set_epi(dat, "cost.disc", at, pop.cost.disc)
  dat <- set_epi(dat, "qaly.disc", at, pop.qaly.disc)
  
  return(dat)
}

# Updated Aging Module --------------------------------------------------------

aging <- function(dat, at) {

  # Update age on attr and also the network
  active <- get_attr(dat, "active")
  idsActive <- which(active == 1)
  age <- get_attr(dat, "age")
  age[idsActive] <- age[idsActive] + 1/52
  dat <- set_attr(dat, "age", age)

  ## Summary statistics ##
  dat <- set_epi(dat, "meanAge", at, mean(age[idsActive], na.rm = TRUE))

  return(dat)
}


# Updated Departure Module -----------------------------------------------------

dfunc <- function(dat, at) {
  
  ## Attributes
  active <- get_attr(dat, "active")
  active.s <- get_attr(dat, "active.s")
  exitTime <- get_attr(dat, "exitTime")
  age <- get_attr(dat, "age")
  status <- get_attr(dat, "status")

  ## Parameters
  death.rates <- get_param(dat, "death.rates")
  end.horizon <- get_param(dat, "end.horizon")
  
  ## Query active
  idsElig <- which(active == 1)
  nElig <- length(idsElig)
  nDeaths <- 0

  if (nElig > 0) {

    ## Calculate age-specific departure rates for each eligible node
    ## Everyone older than 85 gets the final mortality
    whole_ages_of_elig <- pmin(ceiling(age[idsElig]), 101)
    death_rates_of_elig <- death.rates[whole_ages_of_elig]
    
    ## Simulate departure process
    vecDeaths <- which(rbinom(nElig, 1, 1 - exp(-death_rates_of_elig)) == 1)
    idsDeaths <- idsElig[vecDeaths]
    nDeaths <- length(idsDeaths)

    ## Update nodal attributes
    if (nDeaths > 0) {
      active.s[idsDeaths] <- 0
      active[idsDeaths] <- 0
      exitTime[idsDeaths] <- at
    }
    
    ## 65+ who did not die this time-step
    idsRetire <- setdiff(which(age >= 65 & active.s == 1), idsDeaths)
    active.s[idsRetire] <- 0
    
  }
  
  # All individuals become sexually inactive at end horizon
  if (at == end.horizon) {
    active.s <- rep(0, length(active.s))
  }

  ## Reset attr
  dat <- set_attr(dat, "active.s", active.s)
  dat <- set_attr(dat, "active", active)
  dat <- set_attr(dat, "exitTime", exitTime)

  ## Summary statistics
  dat <- set_epi(dat, "d.flow", at, nDeaths)
  
  return(dat)
}


# Updated Arrivals Module ----------------------------------------------------

afunc <- function(dat, at) {
  
  # If the end horizon has been reached, skip arrival module
  end.horizon <- get_param(dat, "end.horizon")
  if (at >= end.horizon) {
    return(dat)
  }

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
    dat <- append_attr(dat, "age", 16, nArrivals)
    dat <- append_attr(dat, "active.s", 1, nArrivals)
  }

  ## Summary statistics
  dat <- set_epi(dat, "a.flow", at, nArrivals)

  return(dat)
}

# Updated Infection Module ---------------------------------------------------

ifunc <- function (dat, at) {
  
  # If the end horizon has been reached, skip infection module
  end.horizon <- get_param(dat, "end.horizon", override.null.error = TRUE)
  if (at >= end.horizon) {
    return(dat)
  }
  
  active.s <- get_attr(dat, "active.s")
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  infTime <- get_attr(dat, "infTime")
  inf.prob <- get_param(dat, "inf.prob")
  act.rate <- get_param(dat, "act.rate")
  inter.eff <- get_param(dat, "inter.eff", override.null.error = TRUE)
  inter.start <- get_param(dat, "inter.start", override.null.error = TRUE)
  
  # Identify which individuals are still sexually active
  idsActive.s <- which(active.s == 1)
  idsInf <- which(active == 1 & status == "i")
  nActive <- sum(active == 1)
  nElig <- length(idsInf)
  nInf <- 0
  if (nElig > 0 && nElig < nActive) {
    del <- discord_edgelist(dat, at, network = 1)
    if (!(is.null(del))) {
      # Remove rows of discordant edge list with an inactive partner
      del <- del[which(del$sus %in% idsActive.s & del$inf %in% idsActive.s),]
      del$infDur <- at - infTime[del$inf]
      del$infDur[del$infDur == 0] <- 1
      linf.prob <- length(inf.prob)
      del$transProb <- ifelse(del$infDur <= linf.prob, 
                              inf.prob[del$infDur], inf.prob[linf.prob])
      if (!is.null(inter.eff) && at >= inter.start) {
        # Apply reduction in transmission probability due to prophylaxis
        del$transProb <- del$transProb * (1 - inter.eff)
      }
      lact.rate <- length(act.rate)
      del$actRate <- ifelse(del$infDur <= lact.rate, act.rate[del$infDur], 
                            act.rate[lact.rate])
      del$finalProb <- 1 - (1 - del$transProb)^del$actRate
      transmit <- rbinom(nrow(del), 1, del$finalProb)
      del <- del[which(transmit == 1), ]
      idsNewInf <- unique(del$sus)
      status <- get_attr(dat, "status")
      status[idsNewInf] <- "i"
      dat <- set_attr(dat, "status", status)
      infTime[idsNewInf] <- at
      dat <- set_attr(dat, "infTime", infTime)
      nInf <- length(idsNewInf)
    }
  }
  if (nInf > 0) {
    dat <- set_transmat(dat, del, at)
  }
  dat <- set_epi(dat, "si.flow", at, nInf)
  return(dat)
}

# Updated Network Resimulation Module ----------------------------------------

resimfunc <- function (dat, at) {
  
  # This function is identical to the default of EpiModel::resim_nets()
  # except for the following line.
  # Skipping over this function after the start of the end horizon greatly
  # improves computational speed.
  end.horizon <- get_param(dat, "end.horizon")
  if (at >= end.horizon) {
    return(dat)
  }
  
  tergmLite <- get_control(dat, "tergmLite")
  isTERGM <- get_control(dat, "isTERGM")
  save.nwstats <- get_control(dat, "save.nwstats")
  resimulate.network <- get_control(dat, "resimulate.network")
  nwstats.formula <- get_control(dat, "nwstats.formula")
  set.control.stergm <- get_control(dat, "set.control.stergm")
  tergmLite.track.duration <- get_control(dat, "tergmLite.track.duration")
  dat <- edges_correct(dat, at)
  active <- get_attr(dat, "active")
  idsActive <- which(active == 1)
  anyActive <- ifelse(length(idsActive) > 0, TRUE, FALSE)
  if (dat$param$groups == 2) {
    group <- get_attr(dat, "group")
    groupids.1 <- which(group == 1)
    groupids.2 <- which(group == 2)
    nActiveG1 <- length(intersect(groupids.1, idsActive))
    nActiveG2 <- length(intersect(groupids.2, idsActive))
    anyActive <- ifelse(nActiveG1 > 0 & nActiveG2 > 0, TRUE, 
                        FALSE)
  }
  nwparam <- get_nwparam(dat)
  if (anyActive == TRUE & resimulate.network == TRUE) {
    if (tergmLite == FALSE) {
      if (isTERGM == TRUE) {
        suppressWarnings(dat$nw[[1]] <- simulate(dat$nw[[1]], 
                                                 formation = nwparam$formation, dissolution = nwparam$coef.diss$dissolution, 
                                                 coef.form = nwparam$coef.form, coef.diss = nwparam$coef.diss$coef.adj, 
                                                 constraints = nwparam$constraints, time.start = at, 
                                                 time.slices = 1, time.offset = 0, monitor = nwstats.formula, 
                                                 control = set.control.stergm))
      }
      else {
        dat$nw[[1]] <- simulate(object = nwparam$formation, 
                                basis = dat$nw[[1]], coef = nwparam$coef.form, 
                                constraints = nwparam$constraints, dynamic = FALSE, 
                                monitor = nwstats.formula, nsim = 1)
      }
      if (save.nwstats == TRUE) {
        new.nwstats <- tail(attributes(dat$nw[[1]])$stats, 
                            1)
        keep.cols <- which(!duplicated(colnames(new.nwstats)))
        new.nwstats <- new.nwstats[, keep.cols, drop = FALSE]
        dat$stats$nwstats[[1]] <- rbind(dat$stats$nwstats[[1]], 
                                        new.nwstats)
      }
    }
    if (tergmLite == TRUE) {
      dat <- tergmLite::updateModelTermInputs(dat)
      if (isTERGM == TRUE) {
        rv <- tergmLite::simulate_network(state = dat$p[[1]]$state, 
                                          coef = c(nwparam$coef.form, nwparam$coef.diss$coef.adj), 
                                          control = dat$control$mcmc.control[[1]], save.changes = TRUE)
        dat$el[[1]] <- rv$el
        if (tergmLite.track.duration == TRUE) {
          dat$p[[1]]$state$nw0 %n% "time" <- rv$state$nw0 %n% 
            "time"
          dat$p[[1]]$state$nw0 %n% "lasttoggle" <- rv$state$nw0 %n% 
            "lasttoggle"
        }
      }
      else {
        rv <- tergmLite::simulate_ergm(state = dat$p[[1]]$state, 
                                       coef = nwparam$coef.form, control = dat$control$mcmc.control[[1]])
        dat$el[[1]] <- rv$el
      }
      if (save.nwstats == TRUE) {
        nwL <- tergmLite::networkLite(dat$el[[1]], dat$attr)
        if (tergmLite.track.duration == TRUE) {
          nwL %n% "time" <- dat$p[[1]]$state$nw0 %n% 
            "time"
          nwL %n% "lasttoggle" <- dat$p[[1]]$state$nw0 %n% 
            "lasttoggle"
        }
        nwstats <- summary(dat$control$nwstats.formulas[[1]], 
                           basis = nwL, term.options = dat$control$mcmc.control[[1]]$term.options, 
                           dynamic = isTERGM)
        keep.cols <- which(!duplicated(names(nwstats)))
        dat$stats$nwstats[[1]] <- rbind(dat$stats$nwstats[[1]], 
                                        nwstats[keep.cols])
      }
    }
  }
  return(dat)
}
