##
## Tuberculosis: Household, Regular, and Casual Contacts Over a Multilayer Network
## EpiModel Gallery (https://github.com/EpiModel/EpiModel-Gallery)
##
## Author: Samuel M. Jenness (Emory University)
## Date: September 2026
##
## The time step is one month. Rates in the parameter list are annual and are
## converted to monthly probabilities with p = 1 - exp(-rate / 12). Modules
## are self-contained (they call no helper functions defined outside them)
## because netsim runs simulations on parallel workers that do not see the
## global environment.
##


# Household placement of newborns ----------------------------------------------

newborn_household <- function(dat, at, new_ids, network) {
  # Passed to netclique() as arrivals.FUN. Each newborn joins the household
  # of a randomly chosen active adult aged 15 to 49, so every child is born
  # into a household with an adult of childbearing age, and households with
  # more such adults receive proportionally more births.
  active <- get_attr(dat, "active")
  age <- get_attr(dat, "age")
  hh_id <- get_attr(dat, "hh_id")
  parents <- setdiff(which(active == 1 & !is.na(hh_id) & age >= 15 & age < 50),
                     new_ids)
  if (length(parents) == 0) {
    return(rep(NA, length(new_ids)))
  }
  hh_id[parents[sample.int(length(parents), length(new_ids), replace = TRUE)]]
}


# Initial TB states --------------------------------------------------------------

init_tb <- function(dat, at) {
  # One-shot setup of the TB attributes. Runs on the first call only (it does
  # nothing when a simulation is resumed from a saved burn-in).
  #
  # status: "s" never infected, "e" early latent infection,
  #         "l" late latent infection, "c" infection cleared (spontaneously
  #         or by preventive therapy), "i" active infectious TB
  #         (undiagnosed), "t" on treatment, "r" recovered from TB (treated
  #         or self-cured)
  # infTime:  step of the most recent infection
  # infType:  "first" or "reinf" for the most recent infection
  # inf_mult: individual infectiousness multiplier for the current or most
  #           recent episode of active TB (gamma distributed, mean 1)
  # txTime:   step treatment started; dxTime: step of diagnosis
  # hhcTime:  step at which the person was last a household contact of a
  #           newly diagnosed case (recorded in every scenario)
  #
  # Initial states are drawn by age from a constant annual risk of infection
  # (init.arti), so older people are more often infected; the burn-in then
  # brings the model to its own equilibrium.
  if (is.null(get_attr(dat, "inf_mult", override.null.error = TRUE))) {
    active <- get_attr(dat, "active")
    age <- get_attr(dat, "age")
    n <- length(active)

    arti <- get_param(dat, "init.arti")
    early.frac <- get_param(dat, "init.early.frac")
    prev.adult <- get_param(dat, "init.prev.adult")
    rec.adult <- get_param(dat, "init.rec.adult")
    inf.k <- get_param(dat, "inf.k")

    status <- rep("s", n)
    infected <- runif(n) < 1 - exp(-arti * age)
    early <- infected & runif(n) < early.frac
    status[infected] <- "l"
    status[early] <- "e"
    adult <- age >= 15
    status[adult & runif(n) < rec.adult] <- "r"
    status[adult & runif(n) < prev.adult] <- "i"

    inf_mult <- rep(NA_real_, n)
    ids_i <- which(status == "i")
    inf_mult[ids_i] <- rgamma(length(ids_i), shape = inf.k, scale = 1 / inf.k)

    infTime <- ifelse(status == "s", NA_integer_, 1L)
    dat <- set_attr(dat, "status", status)
    dat <- set_attr(dat, "infTime", infTime)
    dat <- set_attr(dat, "infType", ifelse(status == "s", NA, "first"))
    dat <- set_attr(dat, "inf_mult", inf_mult)
    dat <- set_attr(dat, "txTime", rep(NA_integer_, n))
    dat <- set_attr(dat, "dxTime", rep(NA_integer_, n))
    dat <- set_attr(dat, "hhcTime", rep(NA_integer_, n))
  }
  return(dat)
}


# Aging --------------------------------------------------------------------------

aging <- function(dat, at) {
  # Age advances by one month per step; agegrp, which the community layers'
  # mixing terms read, follows age.
  age <- get_attr(dat, "age") + 1 / 12
  dat <- set_attr(dat, "age", age)
  dat <- set_attr(dat, "agegrp", ifelse(age < 15, "child", "adult"))
  return(dat)
}


# Infection ------------------------------------------------------------------------

infect <- function(dat, at) {
  # Transmission over every contact layer. On layer k, an edge between an
  # infectious node i and a node j that can be infected transmits in a month
  # with probability
  #
  #   1 - exp(-beta * layer.hours[k] * infness_i * sus_j)
  #
  # where beta is the probability of transmission per contact-hour with an
  # average adult case, layer.hours[k] is the contact time per edge per month
  # on that layer, infness_i is the case's individual multiplier (inf_mult)
  # times the multiplier for its age (young children rarely transmit), and
  # sus_j is 1 for the never infected and for people who have recovered from
  # TB (sus.recovered), and sus.latent for people with late latent infection
  # or a cleared infection: the partial protection that prior infection
  # gives against reinfection. People in early latency, with active TB, or
  # on treatment are not infected again.
  #
  # A node with successful exposures on more than one edge in the same step
  # is infected once, and its infector (and hence the layer) is drawn at
  # random among the successes. Every transmission is recorded with
  # set_transmat().

  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  age <- get_attr(dat, "age")
  inf_mult <- get_attr(dat, "inf_mult")
  infTime <- get_attr(dat, "infTime")
  infType <- get_attr(dat, "infType")

  beta <- get_param(dat, "beta")
  hours <- get_param(dat, "layer.hours")
  layers <- get_param(dat, "layer.names")
  age.breaks <- get_param(dat, "inf.age.breaks")
  age.mult <- get_param(dat, "inf.age.mult")
  sus.latent <- get_param(dat, "sus.latent")
  sus.recovered <- get_param(dat, "sus.recovered")

  n <- length(active)
  infness <- numeric(n)
  ids_i <- which(active == 1 & status == "i")
  infness[ids_i] <- inf_mult[ids_i] *
    age.mult[findInterval(age[ids_i], age.breaks)]
  sus <- numeric(n)
  sus[status == "s"] <- 1
  sus[status %in% c("l", "c")] <- sus.latent
  sus[status == "r"] <- sus.recovered
  sus[active != 1] <- 0

  del <- NULL
  if (any(infness > 0)) {
    for (k in seq_along(hours)) {
      if (hours[k] == 0) next
      el <- get_edgelist(dat, network = k)
      if (NROW(el) == 0) next
      a <- el[, 1]
      b <- el[, 2]
      fwd <- infness[a] > 0 & sus[b] > 0
      rev <- infness[b] > 0 & sus[a] > 0
      inf <- c(a[fwd], b[rev])
      rec <- c(b[fwd], a[rev])
      if (length(inf) == 0) next
      p <- 1 - exp(-beta * hours[k] * infness[inf] * sus[rec])
      hit <- which(runif(length(p)) < p)
      if (length(hit) > 0) {
        del <- rbind(del, data.frame(sus = rec[hit], inf = inf[hit], layer = k))
      }
    }
  }

  n_new <- matrix(0, 2, length(hours), dimnames = list(c("child", "adult"), layers))
  n_type <- c(first = 0, reinf = 0)
  if (!is.null(del)) {
    del <- del[sample.int(nrow(del)), , drop = FALSE]
    del <- del[!duplicated(del$sus), , drop = FALSE]

    prior <- status[del$sus]
    status[del$sus] <- "e"
    infTime[del$sus] <- at
    infType[del$sus] <- ifelse(prior == "s", "first", "reinf")
    dat <- set_attr(dat, "status", status)
    dat <- set_attr(dat, "infTime", infTime)
    dat <- set_attr(dat, "infType", infType)

    del$at <- at
    del$layer <- layers[del$layer]
    del$susAge <- age[del$sus]
    del$infAge <- age[del$inf]
    del$susPrior <- prior
    del$infUid <- get_unique_ids(dat, del$inf)
    dat <- set_transmat(dat, del, at)

    grp <- factor(ifelse(del$susAge < 15, "child", "adult"), c("child", "adult"))
    n_new[] <- table(grp, factor(del$layer, layers))
    n_type[] <- c(sum(prior == "s"), sum(prior != "s"))
  }

  for (k in seq_along(layers)) {
    dat <- set_epi(dat, paste0("inf.", layers[k], ".child"), at, n_new["child", k])
    dat <- set_epi(dat, paste0("inf.", layers[k], ".adult"), at, n_new["adult", k])
  }
  dat <- set_epi(dat, "inf.first", at, n_type[["first"]])
  dat <- set_epi(dat, "inf.reinf", at, n_type[["reinf"]])
  return(dat)
}


# Natural history and passive diagnosis --------------------------------------------

progress <- function(dat, at) {
  # Transitions out of every TB state, with age-specific rates for the bands
  # 0-4, 5-14, and 15+ years where the evidence differs. A state with more
  # than one exit uses competing risks: the probability of leaving in a month
  # is 1 - exp(-(sum of rates) / 12) and the destination is drawn in
  # proportion to the rates. Candidates are taken from a snapshot of status
  # at the start of the module, so no one moves through two states in one
  # step, and infections from this step (infTime == at) stay in early latency
  # for at least one step.
  #
  #   e -> i  early progression          e -> l  stabilization
  #   l -> i  reactivation               l -> c  clearance of infection
  #   i -> t  diagnosis (passive)        i -> death  TB mortality
  #   i -> r  self-cure
  #   t -> r  cure at the end of treatment, t -> i  failure or loss
  #   t -> death  mortality on treatment
  #   r -> i  relapse                    r -> l  stabilization after recovery
  #
  # New episodes of active TB draw the case's infectiousness multiplier from
  # a gamma distribution with mean 1 and shape inf.k. A small shape gives the
  # overdispersion seen in TB: most cases transmit little and a few transmit
  # a great deal.

  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  age <- get_attr(dat, "age")
  infTime <- get_attr(dat, "infTime")
  infType <- get_attr(dat, "infType")
  inf_mult <- get_attr(dat, "inf_mult")
  txTime <- get_attr(dat, "txTime")
  dxTime <- get_attr(dat, "dxTime")
  exitTime <- get_attr(dat, "exitTime")
  hhcTime <- get_attr(dat, "hhcTime")

  prog.early <- get_param(dat, "prog.early")
  stab.early <- get_param(dat, "stab.early")
  react.late <- get_param(dat, "react.late")
  clear.late <- get_param(dat, "clear.late")
  diag.rate <- get_param(dat, "diag.rate")
  tb.mort <- get_param(dat, "tb.mort")
  self.cure <- get_param(dat, "self.cure")
  tx.duration <- get_param(dat, "tx.duration")
  tx.mort <- get_param(dat, "tx.mort")
  tx.success <- get_param(dat, "tx.success")
  relapse.rate <- get_param(dat, "relapse.rate")
  rec.stab <- get_param(dat, "rec.stab")
  inf.k <- get_param(dat, "inf.k")

  band <- findInterval(age, c(0, 5, 15))
  status0 <- status
  n <- length(active)
  u1 <- runif(n)
  u2 <- runif(n)

  # Competing exits for the nodes ids, with one column of annual rates per
  # exit. Returns 0 for no exit this step, else the column of the exit taken.
  exits <- function(ids, rates) {
    rates <- matrix(rates, nrow = length(ids))
    tot <- rowSums(rates)
    out <- integer(length(ids))
    go <- which(u1[ids] < 1 - exp(-tot / 12))
    if (length(go) > 0) {
      out[go] <- 1L
      cum <- 0
      for (j in seq_len(ncol(rates) - 1)) {
        cum <- cum + rates[go, j] / tot[go]
        out[go] <- out[go] + (u2[ids[go]] > cum)
      }
    }
    out
  }

  new_i <- integer(0)
  src_i <- character(0)
  died <- integer(0)

  ## early latent
  ids <- which(active == 1 & status0 == "e" & infTime < at)
  if (length(ids) > 0) {
    x <- exits(ids, cbind(prog.early[band[ids]], stab.early[band[ids]]))
    new_i <- c(new_i, ids[x == 1])
    src_i <- c(src_i, ifelse(infType[ids[x == 1]] == "reinf", "reinf", "recent"))
    status[ids[x == 2]] <- "l"
  }

  ## late latent: reactivation, or clearance of the infection
  ids <- which(active == 1 & status0 == "l")
  if (length(ids) > 0) {
    x <- exits(ids, cbind(react.late[band[ids]], clear.late))
    new_i <- c(new_i, ids[x == 1])
    src_i <- c(src_i, rep("react", sum(x == 1)))
    status[ids[x == 2]] <- "c"
  }

  ## recovered
  ids <- which(active == 1 & status0 == "r")
  if (length(ids) > 0) {
    x <- exits(ids, cbind(rep(relapse.rate, length(ids)), rec.stab))
    new_i <- c(new_i, ids[x == 1])
    src_i <- c(src_i, rep("relapse", sum(x == 1)))
    status[ids[x == 2]] <- "l"
  }

  ## active TB: diagnosis, TB death, self-cure
  ids <- which(active == 1 & status0 == "i")
  n_dx <- c(child = 0, adult = 0)
  if (length(ids) > 0) {
    x <- exits(ids, cbind(diag.rate, tb.mort[band[ids]], self.cure))
    dx <- ids[x == 1]
    status[dx] <- "t"
    txTime[dx] <- at
    dxTime[dx] <- at
    n_dx[] <- c(sum(age[dx] < 15), sum(age[dx] >= 15))
    died <- c(died, ids[x == 2])
    status[ids[x == 3]] <- "r"
  }

  ## on treatment: death during treatment, then cure or return to active TB
  ids <- which(active == 1 & status0 == "t")
  n_fail <- 0
  if (length(ids) > 0) {
    d <- u1[ids] < tx.mort
    died <- c(died, ids[d])
    done <- ids[!d & at - txTime[ids] >= tx.duration]
    cured <- done[u2[done] < tx.success]
    failed <- setdiff(done, cured)
    status[cured] <- "r"
    status[failed] <- "i"
    n_fail <- length(failed)
  }

  ## new episodes of active TB
  if (length(new_i) > 0) {
    status[new_i] <- "i"
    inf_mult[new_i] <- rgamma(length(new_i), shape = inf.k, scale = 1 / inf.k)
  }

  ## TB deaths leave the population at the end of the step
  if (length(died) > 0) {
    active[died] <- 0
    exitTime[died] <- at
    dat <- set_attr(dat, "active", active)
    dat <- set_attr(dat, "exitTime", exitTime)
  }

  dat <- set_attr(dat, "status", status)
  dat <- set_attr(dat, "inf_mult", inf_mult)
  dat <- set_attr(dat, "txTime", txTime)
  dat <- set_attr(dat, "dxTime", dxTime)

  ## incident TB by source and age group, notifications, deaths
  child_i <- age[new_i] < 15
  for (s in c("recent", "reinf", "react", "relapse")) {
    dat <- set_epi(dat, paste0("inc.", s, ".child"), at, sum(src_i == s & child_i))
    dat <- set_epi(dat, paste0("inc.", s, ".adult"), at, sum(src_i == s & !child_i))
  }
  # Onsets among people who were household contacts of a newly diagnosed
  # case in the previous 24 months, by age band at onset
  recent_hhc <- new_i[!is.na(hhcTime[new_i]) & at - hhcTime[new_i] <= 24]
  hhc_band <- cut(age[recent_hhc], c(0, 5, 15, Inf), right = FALSE,
                  labels = c("u5", "5to14", "adult"))
  for (b in levels(hhc_band)) {
    dat <- set_epi(dat, paste0("hhc.inc.", b), at, sum(hhc_band == b))
  }
  dat <- set_epi(dat, "dx.passive.child", at, n_dx[["child"]])
  dat <- set_epi(dat, "dx.passive.adult", at, n_dx[["adult"]])
  dat <- set_epi(dat, "tbdeath.child", at, sum(age[died] < 15))
  dat <- set_epi(dat, "tbdeath.adult", at, sum(age[died] >= 15))
  dat <- set_epi(dat, "tx.fail", at, n_fail)
  return(dat)
}


# Screening, contact investigation, and preventive therapy ---------------------------

screen <- function(dat, at) {
  # Runs after progress(), so the index cases of this step are the people
  # diagnosed this step (dxTime == at). Three optional programs, each off by
  # default:
  #
  #  1. Community-wide screening (acf = 1): every acf.interval months, each
  #     active person aged acf.min.age or older who is not on treatment is
  #     reached with probability acf.coverage, and reached people with active
  #     TB are detected with probability acf.sens and start treatment.
  #  2. Random screening (random.rate > 0): each month, each active person is
  #     screened with probability random.rate, with the same package as a
  #     household contact (below). This is the equal-effort comparator.
  #  3. Household contact investigation (hhci = 1): the household of each
  #     index case diagnosed this step (by any route) is investigated with
  #     probability hhci.reach; each co-resident is evaluated with probability
  #     hhci.eval. Evaluated contacts with active TB are detected with the
  #     age-specific sensitivity screen.sens and start treatment. Evaluated
  #     contacts younger than tpt.age.max who are not found to have TB are
  #     offered preventive therapy (TPT): they start with probability
  #     tpt.start, complete with probability tpt.complete, and a completed
  #     course clears early or late latent infection (status "c") with
  #     probability tpt.eff. A cleared infection keeps the partial
  #     protection of prior infection, but TPT adds none of its own.
  #
  # Whatever the programs, the module records what an investigation of every
  # index household would find (the status of every co-resident at the time
  # of the index diagnosis), which is the model's counterpart of the
  # household contact studies used to check it.

  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  age <- get_attr(dat, "age")
  hh_id <- get_attr(dat, "hh_id")
  txTime <- get_attr(dat, "txTime")
  dxTime <- get_attr(dat, "dxTime")
  hhcTime <- get_attr(dat, "hhcTime")

  acf <- get_param(dat, "acf")
  acf.interval <- get_param(dat, "acf.interval")
  acf.coverage <- get_param(dat, "acf.coverage")
  acf.sens <- get_param(dat, "acf.sens")
  acf.min.age <- get_param(dat, "acf.min.age")
  random.rate <- get_param(dat, "random.rate")
  hhci <- get_param(dat, "hhci")
  hhci.reach <- get_param(dat, "hhci.reach")
  hhci.eval <- get_param(dat, "hhci.eval")
  screen.sens <- get_param(dat, "screen.sens")
  tpt.age.max <- get_param(dat, "tpt.age.max")
  tpt.start <- get_param(dat, "tpt.start")
  tpt.complete <- get_param(dat, "tpt.complete")
  tpt.eff <- get_param(dat, "tpt.eff")

  n <- length(active)
  is_active <- active == 1
  cnt <- c(acf.screened = 0, acf.found = 0, rnd.screened = 0, rnd.found = 0,
           hh.screened = 0, hh.found = 0, tpt.started = 0, tpt.completed = 0,
           tpt.cleared = 0)

  # Screening package for a set of people: detect active TB (starting
  # treatment) and give TPT to eligible people not found to have TB. Returns
  # the updated vectors and counts.
  package <- function(ids, status, txTime, dxTime) {
    sens <- ifelse(age[ids] < 15, screen.sens[1], screen.sens[2])
    found <- ids[status[ids] == "i" & runif(length(ids)) < sens]
    status[found] <- "t"
    txTime[found] <- at
    dxTime[found] <- at
    elig <- setdiff(ids[age[ids] < tpt.age.max & status[ids] != "t"], found)
    started <- elig[runif(length(elig)) < tpt.start]
    completed <- started[runif(length(started)) < tpt.complete]
    latent <- completed[status[completed] %in% c("e", "l")]
    cleared <- latent[runif(length(latent)) < tpt.eff]
    status[cleared] <- "c"
    list(status = status, txTime = txTime, dxTime = dxTime,
         found = length(found), started = length(started),
         completed = length(completed), cleared = length(cleared))
  }

  ## 1. community-wide screening rounds
  if (acf == 1 && (at - 1) %% acf.interval == 0) {
    reached <- which(is_active & age >= acf.min.age & status != "t" &
                       runif(n) < acf.coverage)
    found <- reached[status[reached] == "i" & runif(length(reached)) < acf.sens]
    status[found] <- "t"
    txTime[found] <- at
    dxTime[found] <- at
    cnt[c("acf.screened", "acf.found")] <- c(length(reached), length(found))
  }

  ## 2. random screening with the household package
  if (random.rate > 0) {
    ids <- which(is_active & status != "t" & runif(n) < random.rate)
    out <- package(ids, status, txTime, dxTime)
    status <- out$status; txTime <- out$txTime; dxTime <- out$dxTime
    cnt[c("rnd.screened", "rnd.found", "tpt.started", "tpt.completed",
          "tpt.cleared")] <- c(length(ids), out$found, out$started,
                               out$completed, out$cleared)
  }

  ## 3. index households: what they hold, and the investigation
  index <- which(is_active & !is.na(dxTime) & dxTime == at)
  yield <- matrix(0, 3, 3, dimnames = list(c("u5", "5to14", "adult"),
                                           c("n", "active", "infected")))
  if (length(index) > 0) {
    idx_hh <- unique(hh_id[index])
    contacts <- which(is_active & hh_id %in% idx_hh & !is.na(hh_id))
    contacts <- setdiff(contacts, index)
    band <- cut(age[contacts], c(0, 5, 15, Inf), right = FALSE,
                labels = rownames(yield))
    yield[, "n"] <- table(band)
    yield[, "active"] <- table(band[status[contacts] == "i"])
    yield[, "infected"] <- table(band[status[contacts] %in% c("e", "l", "c", "r")])
    hhcTime[contacts] <- at

    if (hhci == 1 && length(contacts) > 0) {
      reached_hh <- idx_hh[runif(length(idx_hh)) < hhci.reach]
      ev <- contacts[hh_id[contacts] %in% reached_hh & status[contacts] != "t"]
      ev <- ev[runif(length(ev)) < hhci.eval]
      out <- package(ev, status, txTime, dxTime)
      status <- out$status; txTime <- out$txTime; dxTime <- out$dxTime
      cnt["hh.screened"] <- length(ev)
      cnt["hh.found"] <- out$found
      cnt[c("tpt.started", "tpt.completed", "tpt.cleared")] <-
        cnt[c("tpt.started", "tpt.completed", "tpt.cleared")] +
        c(out$started, out$completed, out$cleared)
    }
  }

  dat <- set_attr(dat, "status", status)
  dat <- set_attr(dat, "txTime", txTime)
  dat <- set_attr(dat, "dxTime", dxTime)
  dat <- set_attr(dat, "hhcTime", hhcTime)
  for (nm in names(cnt)) {
    dat <- set_epi(dat, nm, at, cnt[[nm]])
  }
  dat <- set_epi(dat, "index", at, length(index))
  for (b in rownames(yield)) {
    for (v in colnames(yield)) {
      dat <- set_epi(dat, paste0("hhc.", v, ".", b), at, yield[b, v])
    }
  }
  return(dat)
}


# Background mortality -------------------------------------------------------------

deaths <- function(dat, at) {
  # Age-specific background mortality; annual rates by the age bands that
  # start at mort.breaks. Everyone who reaches age 100 dies.
  active <- get_attr(dat, "active")
  age <- get_attr(dat, "age")
  exitTime <- get_attr(dat, "exitTime")
  mort.rates <- get_param(dat, "mort.rates")
  mort.breaks <- get_param(dat, "mort.breaks")

  ids <- which(active == 1)
  p <- 1 - exp(-mort.rates[findInterval(age[ids], mort.breaks)] / 12)
  dead <- ids[runif(length(ids)) < p | age[ids] >= 100]
  if (length(dead) > 0) {
    active[dead] <- 0
    exitTime[dead] <- at
    dat <- set_attr(dat, "active", active)
    dat <- set_attr(dat, "exitTime", exitTime)
  }
  dat <- set_epi(dat, "deaths", at, length(dead))
  return(dat)
}


# Births ----------------------------------------------------------------------------

births <- function(dat, at) {
  # Births at a constant per-capita rate. Every attribute a module reads is
  # appended here; an attribute left out would be filled by EpiModel with
  # values sampled from the current population. The household id is left
  # missing: the clique layer's arrival rule (newborn_household) places each
  # newborn in a household.
  active <- get_attr(dat, "active")
  birth.rate <- get_param(dat, "birth.rate")
  n_new <- rpois(1, birth.rate / 12 * sum(active == 1))
  if (n_new > 0) {
    dat <- append_core_attr(dat, at, n_new)
    dat <- append_attr(dat, "status", "s", n_new)
    dat <- append_attr(dat, "age", 0, n_new)
    dat <- append_attr(dat, "agegrp", "child", n_new)
    dat <- append_attr(dat, "hh_id", NA, n_new)
    dat <- append_attr(dat, "infTime", NA_integer_, n_new)
    dat <- append_attr(dat, "infType", NA_character_, n_new)
    dat <- append_attr(dat, "inf_mult", NA_real_, n_new)
    dat <- append_attr(dat, "txTime", NA_integer_, n_new)
    dat <- append_attr(dat, "dxTime", NA_integer_, n_new)
    dat <- append_attr(dat, "hhcTime", NA_integer_, n_new)
  }
  dat <- set_epi(dat, "births", at, n_new)
  return(dat)
}


# Household formation ----------------------------------------------------------------

households <- function(dat, at) {
  # Runs after nwupdate, when births have been placed and the dead removed.
  # Two kinds of moves keep the household structure stable over decades:
  #
  #  1. Leaving home: each month, people aged leave.ages[1] to leave.ages[2]
  #     who live in a household of three or more, with at least one other
  #     adult, leave with probability leave.rate. Leavers are paired at
  #     random into new two-person households (an odd leaver lives alone).
  #  2. Children left without an adult: children whose household no longer
  #     has anyone aged 15 or older move, together, into the household of a
  #     randomly chosen adult.
  #
  # move_to_group() changes the household id and rewires the clique layer in
  # one step; changing hh_id with set_attr() alone would leave the old
  # household edges in place.
  active <- get_attr(dat, "active")
  age <- get_attr(dat, "age")
  hh_id <- get_attr(dat, "hh_id")
  leave.rate <- get_param(dat, "leave.rate")
  leave.ages <- get_param(dat, "leave.ages")

  on <- which(active == 1 & !is.na(hh_id))
  hh <- match(hh_id, unique(hh_id[on]))
  size <- tabulate(hh[on])[hh]
  adults <- tabulate(hh[on][age[on] >= 15], nbins = max(hh[on]))[hh]
  elig <- which(active == 1 & !is.na(hh_id) & age >= leave.ages[1] &
                  age < leave.ages[2] & size >= 3 & adults >= 2)
  leavers <- elig[runif(length(elig)) < leave.rate]
  if (length(leavers) > 0) {
    leavers <- leavers[sample.int(length(leavers))]
    new_ids <- max(hh_id, na.rm = TRUE) + ceiling(seq_along(leavers) / 2)
    dat <- move_to_group(dat, ids = leavers, group = new_ids)
    hh_id <- get_attr(dat, "hh_id")
  }

  has_adult <- tapply(age[on] >= 15, hh_id[on], any)
  orphans <- which(active == 1 & age < 15 & !is.na(hh_id) &
                     !has_adult[as.character(hh_id)])
  if (length(orphans) > 0) {
    adults <- which(active == 1 & age >= 15 & !is.na(hh_id))
    old_hh <- unique(hh_id[orphans])
    new_hh <- hh_id[adults[sample.int(length(adults), length(old_hh),
                                      replace = TRUE)]]
    dat <- move_to_group(dat, ids = orphans,
                         group = new_hh[match(hh_id[orphans], old_hh)])
  }

  dat <- set_epi(dat, "hh.leave", at, length(leavers))
  dat <- set_epi(dat, "hh.orphan", at, length(orphans))
  return(dat)
}


# Summary statistics -------------------------------------------------------------------

tally <- function(dat, at) {
  # Population counts by TB state and age group, and the household
  # structure, recorded at the end of each step.
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  age <- get_attr(dat, "age")
  hh_id <- get_attr(dat, "hh_id")

  on <- active == 1
  child <- age < 15
  for (s in c("s", "e", "l", "c", "i", "t", "r")) {
    dat <- set_epi(dat, paste0(s, ".num.child"), at, sum(on & child & status == s))
    dat <- set_epi(dat, paste0(s, ".num.adult"), at, sum(on & !child & status == s))
  }
  dat <- set_epi(dat, "num.child", at, sum(on & child))
  dat <- set_epi(dat, "num.adult", at, sum(on & !child))
  dat <- set_epi(dat, "num.u5", at, sum(on & age < 5))
  infected <- status %in% c("e", "l", "c", "r")
  dat <- set_epi(dat, "infected.u5", at, sum(on & age < 5 & infected))
  dat <- set_epi(dat, "num.15to34", at, sum(on & age >= 15 & age < 35))
  dat <- set_epi(dat, "infected.15to34", at,
                 sum(on & age >= 15 & age < 35 & infected))

  sizes <- tabulate(match(hh_id[on], unique(hh_id[on])))
  dat <- set_epi(dat, "hh.num", at, length(sizes))
  dat <- set_epi(dat, "hh.alone", at, sum(sizes == 1))
  return(dat)
}
