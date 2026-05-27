
##
## HIV Transmission Model with Care Cascade and PrEP
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##


# Disease progression module -----------------------------------------------
#
# Custom nodal attributes used by these modules:
#
#   stage         "acute" / "chronic" / "aids" / NA (susceptible)
#   stage.time    timesteps since current stage entry
#   diag.status   0 (HIV-negative or undiagnosed) / 1 (diagnosed)
#   art.status    0 / 1 (currently on ART)
#   vl.supp       0 / 1 (virally suppressed; meaningful only on ART)
#   art.time      timestep at which the current ART episode began
#   prep.status   0 / 1 (on PrEP; meaningful only for susceptibles)
#
# All of these are initialized inside progress() on its first call.

progress <- function(dat, at) {
  # Acute -> Chronic -> AIDS progression.
  #
  # Uses a snapshot of stage at function entry. All eligibility checks read
  # from stage0, so a node cannot cascade through two stages in one step,
  # regardless of progression rates. Same pattern as the RSV example.
  #
  # Acute -> Chronic is gated on infTime < at so brand-new infections from
  # this step's infect() module always spend at least one step in acute.

  # ---- Attribute initializer (runs once on first call) ----
  # Done inline rather than via a separate init module so the helper is
  # serialized into each parallel worker along with progress() itself.
  if (is.null(get_attr(dat, "stage", override.null.error = TRUE))) {
    n <- length(get_attr(dat, "active"))
    .active <- get_attr(dat, "active")
    .status <- get_attr(dat, "status")

    # Seed-infected stage distribution proportional to mean stage duration
    # (acute 12 wk, chronic 520 wk, AIDS 104 wk) -- approximates the
    # stage-equilibrium distribution and avoids a synchronous cohort wave.
    stage0_init <- rep(NA_character_, n)
    stage.time0_init <- rep(NA_integer_, n)
    inf_ids <- which(.active == 1 & .status == "i")
    if (length(inf_ids) > 0) {
      durs <- c(acute = 12, chronic = 520, aids = 104)
      probs <- durs / sum(durs)
      stage0_init[inf_ids] <- sample(names(probs), length(inf_ids),
                                     replace = TRUE, prob = probs)
      stage.time0_init[inf_ids] <- 0L
    }

    dat <- set_attr(dat, "stage", stage0_init)
    dat <- set_attr(dat, "stage.time", stage.time0_init)
    dat <- set_attr(dat, "diag.status", rep(0L, n))
    dat <- set_attr(dat, "art.status", rep(0L, n))
    dat <- set_attr(dat, "vl.supp", rep(0L, n))
    dat <- set_attr(dat, "art.time", rep(NA_integer_, n))
    # prep.status starts NA; the prep() module replaces NA values on
    # its first call with an indication-based initial coverage, since
    # we need the simulated network to determine who is "high-risk"
    # and the network isn't accessible here.
    dat <- set_attr(dat, "prep.status", rep(NA_integer_, n))
  }

  active <- get_attr(dat, "active")
  stage <- get_attr(dat, "stage")
  stage.time <- get_attr(dat, "stage.time")
  art.status <- get_attr(dat, "art.status")
  vl.supp <- get_attr(dat, "vl.supp")
  infTime <- get_attr(dat, "infTime")

  acute.to.chronic.rate <- get_param(dat, "acute.to.chronic.rate")
  chronic.to.aids.rate <- get_param(dat, "chronic.to.aids.rate")
  art.prog.mult <- get_param(dat, "art.prog.mult")

  stage0 <- stage

  # Increment stage clock for everyone currently in a stage.
  in_stage <- !is.na(stage.time) & active == 1
  stage.time[in_stage] <- stage.time[in_stage] + 1L

  ## Acute -> Chronic ##
  ids_ac <- which(active == 1 & !is.na(stage0) & stage0 == "acute" &
                    !is.na(infTime) & infTime < at)
  n_ac <- 0
  if (length(ids_ac) > 0) {
    # Suppressive ART slows progression (a coarse pedagogical stand-in for
    # the real effect of HIV control on CD4 dynamics).
    supp <- art.status[ids_ac] == 1 & vl.supp[ids_ac] == 1
    rates <- ifelse(supp, acute.to.chronic.rate * art.prog.mult,
                          acute.to.chronic.rate)
    hits <- rbinom(length(ids_ac), 1, rates) == 1
    new_chronic <- ids_ac[hits]
    if (length(new_chronic) > 0) {
      stage[new_chronic] <- "chronic"
      stage.time[new_chronic] <- 0L
      n_ac <- length(new_chronic)
    }
  }

  ## Chronic -> AIDS ##
  ids_ca <- which(active == 1 & !is.na(stage0) & stage0 == "chronic")
  n_ca <- 0
  if (length(ids_ca) > 0) {
    supp <- art.status[ids_ca] == 1 & vl.supp[ids_ca] == 1
    rates <- ifelse(supp, chronic.to.aids.rate * art.prog.mult,
                          chronic.to.aids.rate)
    hits <- rbinom(length(ids_ca), 1, rates) == 1
    new_aids <- ids_ca[hits]
    if (length(new_aids) > 0) {
      stage[new_aids] <- "aids"
      stage.time[new_aids] <- 0L
      n_ca <- length(new_aids)
    }
  }

  dat <- set_attr(dat, "stage", stage)
  dat <- set_attr(dat, "stage.time", stage.time)

  is_act <- active == 1
  dat <- set_epi(dat, "ac.flow", at, n_ac)
  dat <- set_epi(dat, "ca.flow", at, n_ca)
  dat <- set_epi(dat, "acute.num", at,
                 sum(is_act & stage %in% "acute"))
  dat <- set_epi(dat, "chronic.num", at,
                 sum(is_act & stage %in% "chronic"))
  dat <- set_epi(dat, "aids.num", at,
                 sum(is_act & stage %in% "aids"))

  return(dat)
}


# Care cascade module ------------------------------------------------------

cascade <- function(dat, at) {
  # Tracks HIV+ individuals through the care continuum:
  #
  #   undiagnosed --(test.rate, or aids.dx.rate for AIDS)--> diagnosed
  #   diagnosed   --(linkage.rate)--> on ART (not yet suppressed)
  #   on ART      --(suppression.rate)--> virally suppressed
  #   on ART      --(art.disc.rate)--> off ART (returns to diagnosed)
  #
  # Snapshot pattern: all eligibility checks read from diag0/art0/supp0,
  # so no node moves more than one cascade step per timestep.

  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  stage <- get_attr(dat, "stage")
  diag.status <- get_attr(dat, "diag.status")
  art.status <- get_attr(dat, "art.status")
  vl.supp <- get_attr(dat, "vl.supp")
  art.time <- get_attr(dat, "art.time")

  test.rate <- get_param(dat, "test.rate")
  aids.dx.rate <- get_param(dat, "aids.dx.rate")
  linkage.rate <- get_param(dat, "linkage.rate")
  art.reinit.rate <- get_param(dat, "art.reinit.rate")
  suppression.rate <- get_param(dat, "suppression.rate")
  art.disc.rate <- get_param(dat, "art.disc.rate")

  diag0 <- diag.status
  art0 <- art.status
  supp0 <- vl.supp

  is_inf <- active == 1 & status == "i"

  ## 1. Testing: undiagnosed -> diagnosed.
  ## AIDS-stage infections are diagnosed at a higher rate to represent
  ## symptom-driven presentation to care; acute/chronic test at the
  ## scenario's background screening rate.
  ids_undx <- which(is_inf & diag0 == 0)
  if (length(ids_undx) > 0) {
    rates <- ifelse(!is.na(stage[ids_undx]) & stage[ids_undx] == "aids",
                    aids.dx.rate, test.rate)
    hits <- rbinom(length(ids_undx), 1, rates) == 1
    new_dx <- ids_undx[hits]
    if (length(new_dx) > 0) diag.status[new_dx] <- 1L
  }

  ## 2. Linkage / re-engagement: diagnosed, not on ART -> on ART.
  ## Newly-diagnosed (no prior ART history, art.time == NA) link at
  ## linkage.rate. Previously-on-ART individuals (have an art.time
  ## from a prior ART episode) re-engage at art.reinit.rate. This
  ## reflects the empirical pattern that initial linkage tends to be
  ## faster than re-engagement after disengagement from care.
  ids_link_new <- which(is_inf & diag0 == 1 & art0 == 0 &
                          is.na(art.time))
  if (length(ids_link_new) > 0) {
    hits <- rbinom(length(ids_link_new), 1, linkage.rate) == 1
    new_art <- ids_link_new[hits]
    if (length(new_art) > 0) {
      art.status[new_art] <- 1L
      art.time[new_art] <- at
    }
  }

  ids_link_re <- which(is_inf & diag0 == 1 & art0 == 0 &
                         !is.na(art.time))
  if (length(ids_link_re) > 0) {
    hits <- rbinom(length(ids_link_re), 1, art.reinit.rate) == 1
    reinit <- ids_link_re[hits]
    if (length(reinit) > 0) {
      art.status[reinit] <- 1L
      art.time[reinit] <- at
    }
  }

  ## 3. Suppression: on ART, not yet suppressed -> virally suppressed.
  ## Reflects the 2-3 month delay between ART initiation and undetectable
  ## viral load; until suppressed, the partial-suppression infectiousness
  ## multiplier applies in infect().
  ids_sup <- which(is_inf & art0 == 1 & supp0 == 0)
  if (length(ids_sup) > 0) {
    hits <- rbinom(length(ids_sup), 1, suppression.rate) == 1
    new_supp <- ids_sup[hits]
    if (length(new_supp) > 0) vl.supp[new_supp] <- 1L
  }

  ## 4. Discontinuation: on ART -> off ART (back to diagnosed only).
  ## Disengagement clears art.status and vl.supp but leaves art.time
  ## intact so that ids_link_re picks the node up next step at the
  ## re-engagement rate.
  ids_disc <- which(is_inf & art0 == 1)
  if (length(ids_disc) > 0) {
    hits <- rbinom(length(ids_disc), 1, art.disc.rate) == 1
    disc <- ids_disc[hits]
    if (length(disc) > 0) {
      art.status[disc] <- 0L
      vl.supp[disc] <- 0L
    }
  }

  dat <- set_attr(dat, "diag.status", diag.status)
  dat <- set_attr(dat, "art.status", art.status)
  dat <- set_attr(dat, "vl.supp", vl.supp)
  dat <- set_attr(dat, "art.time", art.time)

  # Compartment counts and per-step flows
  dat <- set_epi(dat, "dx.flow", at, sum(diag.status == 1 & diag0 == 0))
  dat <- set_epi(dat, "link.flow", at, sum(art.status == 1 & art0 == 0))
  dat <- set_epi(dat, "supp.flow", at, sum(vl.supp == 1 & supp0 == 0))
  dat <- set_epi(dat, "disc.flow", at, sum(art.status == 0 & art0 == 1))
  dat <- set_epi(dat, "dx.num", at, sum(is_inf & diag.status == 1))
  dat <- set_epi(dat, "art.num", at, sum(is_inf & art.status == 1))
  dat <- set_epi(dat, "supp.num", at, sum(is_inf & vl.supp == 1))

  return(dat)
}


# PrEP module --------------------------------------------------------------

prep <- function(dat, at) {
  # PrEP among indicated HIV-negative individuals. Indication mirrors
  # the CDC behavioral-risk criteria for PrEP eligibility in 2025:
  #
  #   - High partnership activity: total current degree across the main
  #     and casual layers >= prep.indic.deg (default 2, meaning at
  #     least one concurrent partner).
  #   - Serodifferent partnership: at least one current partner is HIV+.
  #
  # The indication is recomputed each step. Initiation applies only to
  # currently-indicated susceptibles; discontinuation applies to all
  # current PrEP users regardless of current indication, since real-
  # world drop-off is driven by adherence and access more than by
  # changes in indication.
  #
  # On the very first call, prep.status is NA (initialized in progress);
  # we replace those NAs with an indication-based initial coverage now
  # that the simulated network is available.

  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  prep.status <- get_attr(dat, "prep.status")

  prep.init.cov <- get_param(dat, "prep.init.cov")
  prep.start.rate <- get_param(dat, "prep.start.rate")
  prep.stop.rate <- get_param(dat, "prep.stop.rate")
  prep.indic.deg <- get_param(dat, "prep.indic.deg")

  ## ---- Compute PrEP indication ---- ##
  n <- length(active)
  total_deg <- integer(n)
  has_pos_partner <- logical(n)

  for (k in 1:2) {
    el <- get_edgelist(dat, network = k)
    if (is.null(el) || nrow(el) == 0) next
    total_deg <- total_deg + get_degree(el)
    has_pos_partner[el[, 2][status[el[, 1]] == "i"]] <- TRUE
    has_pos_partner[el[, 1][status[el[, 2]] == "i"]] <- TRUE
  }
  indicated <- active == 1 & status == "s" &
               (total_deg >= prep.indic.deg | has_pos_partner)

  ## ---- First-call initialization (NA -> 0 / 1) ---- ##
  needs_init <- is.na(prep.status) & active == 1
  if (any(needs_init)) {
    prep.status[needs_init] <- 0L
    if (prep.init.cov > 0) {
      ids_init <- which(needs_init & indicated)
      if (length(ids_init) > 0) {
        prep.status[ids_init] <- rbinom(length(ids_init), 1, prep.init.cov)
      }
    }
  }

  ## ---- Initiation among indicated susceptibles ---- ##
  if (prep.start.rate > 0) {
    ids_off <- which(indicated & prep.status == 0)
    if (length(ids_off) > 0) {
      hits <- rbinom(length(ids_off), 1, prep.start.rate) == 1
      prep.status[ids_off[hits]] <- 1L
    }
  }

  ## ---- Discontinuation among current users ---- ##
  if (prep.stop.rate > 0) {
    ids_on <- which(active == 1 & status == "s" & prep.status == 1)
    if (length(ids_on) > 0) {
      hits <- rbinom(length(ids_on), 1, prep.stop.rate) == 1
      prep.status[ids_on[hits]] <- 0L
    }
  }

  dat <- set_attr(dat, "prep.status", prep.status)

  is_sus <- active == 1 & status == "s"
  dat <- set_epi(dat, "prep.num", at, sum(is_sus & prep.status == 1))
  dat <- set_epi(dat, "prep.indic.num", at, sum(indicated))
  dat <- set_epi(dat, "prep.num.indic", at,
                 sum(indicated & prep.status == 1))

  return(dat)
}


# Infection module ---------------------------------------------------------

infect <- function(dat, at) {
  # Multilayer HIV transmission on the main and casual partnership
  # networks. Per-edge probability of transmission per timestep is:
  #
  #   p_act = inf.prob.act * stage_mult * art_mult
  #   p_edge = 1 - (1 - p_act)^acts_per_partnership
  #
  # Susceptible-side PrEP reduces p_edge by (1 - prep.efficacy).
  #
  # stage_mult is acute/aids/chronic relative infectiousness; art_mult
  # captures U=U: suppressed individuals are essentially non-infectious,
  # and unsuppressed-on-ART are partially reduced.

  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  stage <- get_attr(dat, "stage")
  art.status <- get_attr(dat, "art.status")
  vl.supp <- get_attr(dat, "vl.supp")
  prep.status <- get_attr(dat, "prep.status")
  infTime <- get_attr(dat, "infTime")
  stage.time <- get_attr(dat, "stage.time")

  inf.prob.act <- get_param(dat, "inf.prob.act")
  rel.inf.acute <- get_param(dat, "rel.inf.acute")
  rel.inf.aids <- get_param(dat, "rel.inf.aids")
  rel.inf.art.unsupp <- get_param(dat, "rel.inf.art.unsupp")
  rel.inf.art.supp <- get_param(dat, "rel.inf.art.supp")
  prep.efficacy <- get_param(dat, "prep.efficacy")
  acts.main <- get_param(dat, "acts.main")
  acts.casual <- get_param(dat, "acts.casual")

  layer_acts <- c(acts.main, acts.casual)
  all_new <- integer(0)

  for (k in 1:2) {
    # discord_edgelist returns one row per S-I pair across active nodes
    # on the requested network layer (columns: at, sus, inf).
    del <- discord_edgelist(dat, at, network = k)
    if (is.null(del) || nrow(del) == 0) next

    stg <- stage[del$inf]
    art <- art.status[del$inf]
    sup <- vl.supp[del$inf]

    stage_mult <- ifelse(!is.na(stg) & stg == "acute", rel.inf.acute,
                  ifelse(!is.na(stg) & stg == "aids",  rel.inf.aids, 1))
    art_mult <- ifelse(art == 1 & sup == 1, rel.inf.art.supp,
                ifelse(art == 1 & sup == 0, rel.inf.art.unsupp, 1))

    p_act <- inf.prob.act * stage_mult * art_mult

    # Susceptible-side PrEP reduces acquisition probability.
    p_act <- p_act * ifelse(prep.status[del$sus] == 1, 1 - prep.efficacy, 1)
    p_act <- pmax(pmin(p_act, 1), 0)

    p_edge <- 1 - (1 - p_act)^layer_acts[k]

    transmit <- rbinom(length(p_edge), 1, p_edge) == 1
    if (any(transmit)) all_new <- c(all_new, del$sus[transmit])
  }

  all_new <- unique(all_new)
  n_new <- length(all_new)

  if (n_new > 0) {
    status[all_new] <- "i"
    stage[all_new] <- "acute"
    stage.time[all_new] <- 0L
    infTime[all_new] <- at
    # PrEP becomes meaningless after seroconversion; clear so prep.num
    # tracks only HIV-negative users.
    prep.status[all_new] <- 0L

    dat <- set_attr(dat, "status", status)
    dat <- set_attr(dat, "stage", stage)
    dat <- set_attr(dat, "stage.time", stage.time)
    dat <- set_attr(dat, "infTime", infTime)
    dat <- set_attr(dat, "prep.status", prep.status)
  }

  dat <- set_epi(dat, "si.flow", at, n_new)

  return(dat)
}


# Departure module ---------------------------------------------------------

dfunc <- function(dat, at) {
  # Two competing departure processes:
  #   - Background: applies to all non-AIDS individuals at departure.rate.
  #   - AIDS-related: applies to AIDS-stage individuals at the elevated
  #     aids.depart.rate, reduced by art.aids.surv.mult for those on
  #     suppressive ART (modeling extended survival on treatment).

  active <- get_attr(dat, "active")
  stage <- get_attr(dat, "stage")
  art.status <- get_attr(dat, "art.status")
  vl.supp <- get_attr(dat, "vl.supp")
  exitTime <- get_attr(dat, "exitTime")

  departure.rate <- get_param(dat, "departure.rate")
  aids.depart.rate <- get_param(dat, "aids.depart.rate")
  art.aids.surv.mult <- get_param(dat, "art.aids.surv.mult")

  ## Background departure (all active, non-AIDS) ##
  ids_bg <- which(active == 1 & (is.na(stage) | stage != "aids"))
  bg_dep <- ids_bg[rbinom(length(ids_bg), 1, departure.rate) == 1]

  ## AIDS-related departure ##
  ids_aids <- which(active == 1 & !is.na(stage) & stage == "aids")
  aids_dep <- integer(0)
  if (length(ids_aids) > 0) {
    on_supp <- art.status[ids_aids] == 1 & vl.supp[ids_aids] == 1
    rates <- ifelse(on_supp, aids.depart.rate * art.aids.surv.mult,
                             aids.depart.rate)
    aids_dep <- ids_aids[rbinom(length(ids_aids), 1, rates) == 1]
  }

  all_dep <- c(bg_dep, aids_dep)
  if (length(all_dep) > 0) {
    active[all_dep] <- 0L
    exitTime[all_dep] <- at
    dat <- set_attr(dat, "active", active)
    dat <- set_attr(dat, "exitTime", exitTime)
  }

  dat <- set_epi(dat, "dep.flow", at, length(all_dep))
  dat <- set_epi(dat, "aids.dep.flow", at, length(aids_dep))

  return(dat)
}


# Arrival module -----------------------------------------------------------

afunc <- function(dat, at) {
  # New susceptibles arrive at a per-capita Poisson rate. Each custom
  # attribute vector is extended in lockstep with the EpiModel core
  # attributes via append_core_attr() + append_attr().

  active <- get_attr(dat, "active")
  a.rate <- get_param(dat, "arrival.rate")
  nArr <- rpois(1, sum(active == 1) * a.rate)

  if (nArr > 0) {
    dat <- append_core_attr(dat, at, nArr)
    dat <- append_attr(dat, "status", "s", nArr)
    dat <- append_attr(dat, "stage", NA_character_, nArr)
    dat <- append_attr(dat, "stage.time", NA_integer_, nArr)
    dat <- append_attr(dat, "infTime", NA_integer_, nArr)
    dat <- append_attr(dat, "diag.status", 0L, nArr)
    dat <- append_attr(dat, "art.status", 0L, nArr)
    dat <- append_attr(dat, "vl.supp", 0L, nArr)
    dat <- append_attr(dat, "art.time", NA_integer_, nArr)
    dat <- append_attr(dat, "prep.status", 0L, nArr)
  }

  dat <- set_epi(dat, "arr.flow", at, nArr)
  return(dat)
}
