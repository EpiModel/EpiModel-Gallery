
##
## RSV: Age-Stratified SEIR Over a Multilayer Network
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##
## Authors: Samuel M. Jenness (Emory University)
## Date: May 2026
##


# Attribute initializer ----------------------------------------------------

init_attrs <- function(dat, at) {
  # Set up nodal attributes used by the disease and intervention modules.
  # Runs once on the first call. The simulation's `active` and `status`
  # attributes are managed by EpiModel itself.
  #
  # vax_status:
  #   NA            -- not vaccinated
  #   "elderly_vax" -- received the older-adult vaccine (one-shot, leaky)
  #   "infant_proph"-- received passive antibody prophylaxis (infants only)
  #
  # inf_stage:
  #   NA   -- not infectious (status %in% c("s","e","r"))
  #   "ip" -- presymptomatic infectious
  #   "is" -- symptomatic infectious
  #   "ia" -- asymptomatic infectious

  if (is.null(get_attr(dat, "vax_status", override.null.error = TRUE))) {
    active <- get_attr(dat, "active")
    status <- get_attr(dat, "status")
    age <- get_attr(dat, "age")

    elderly.cov <- get_param(dat, "elderly.vax.coverage")
    infant.cov <- get_param(dat, "infant.proph.coverage")

    vax_status <- rep(NA_character_, length(active))
    inf_stage <- rep(NA_character_, length(active))

    # Place seed infections into the presymptomatic infectious stage so
    # they participate in the progression module and eventually recover.
    # Without this, init.net's i.num seeds would remain status = "i" with
    # inf_stage = NA forever (never picked up by progress()'s stage-keyed
    # transitions) and act as a permanent infectious reservoir.
    inf_stage[active == 1 & status == "i"] <- "ip"

    elig_e <- which(active == 1 & age == "elderly")
    if (length(elig_e) > 0 && elderly.cov > 0) {
      hit <- which(rbinom(length(elig_e), 1, elderly.cov) == 1)
      vax_status[elig_e[hit]] <- "elderly_vax"
    }
    elig_i <- which(active == 1 & age == "infant")
    if (length(elig_i) > 0 && infant.cov > 0) {
      hit <- which(rbinom(length(elig_i), 1, infant.cov) == 1)
      vax_status[elig_i[hit]] <- "infant_proph"
    }

    dat <- set_attr(dat, "vax_status", vax_status)
    dat <- set_attr(dat, "inf_stage", inf_stage)
  }

  return(dat)
}


# Infection module ---------------------------------------------------------

infect <- function(dat, at) {
  # S -> E transmission. For each of the two network layers we walk the
  # edgelist directly (instead of using discord_edgelist) so we can apply
  # per-layer transmission probabilities and per-layer NPI effects. The
  # layers are: 1 = family/close contacts, 2 = community.
  #
  # Per-edge transmission probability is the layer's per-step `inf.prob`
  # multiplied by:
  #   - asymp.inf.mult (if the infectious partner is asymptomatic)
  #   - (1 - npi.mask.efficacy) on the community layer when NPIs are active
  #   - (1 - vax efficacy) on the susceptible side (elderly or infant)

  ## Attributes ##
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  inf_stage <- get_attr(dat, "inf_stage")
  vax_status <- get_attr(dat, "vax_status")
  infTime <- get_attr(dat, "infTime")

  ## Parameters ##
  ip.fam <- get_param(dat, "inf.prob.family")
  ip.com <- get_param(dat, "inf.prob.community")
  asymp.mult <- get_param(dat, "asymp.inf.mult")

  npi.start <- get_param(dat, "npi.start")
  npi.end <- get_param(dat, "npi.end")
  npi.mask.eff <- get_param(dat, "npi.mask.efficacy")
  npi.contact.mult <- get_param(dat, "npi.contact.mult")

  vax.eff.elderly <- get_param(dat, "elderly.vax.efficacy")
  vax.eff.infant <- get_param(dat, "infant.proph.efficacy")

  npi.on <- (at >= npi.start && at <= npi.end)
  com.prob.mult <- if (npi.on) (1 - npi.mask.eff) else 1
  # npi.contact.mult applies to community contact RATE; we approximate by
  # thinning community edges -- include each one with this probability.
  com.contact.mult <- if (npi.on) npi.contact.mult else 1

  layer_probs <- c(ip.fam, ip.com * com.prob.mult)
  all_new <- integer(0)

  for (k in 1:2) {
    el <- dat$run$el[[k]]
    if (is.null(el) || nrow(el) == 0) next

    head <- el[, 1]
    tail <- el[, 2]

    # Optional community-edge thinning under NPIs
    if (k == 2 && com.contact.mult < 1) {
      keep <- which(rbinom(nrow(el), 1, com.contact.mult) == 1)
      head <- head[keep]; tail <- tail[keep]
      if (length(head) == 0) next
    }

    h_status <- status[head]
    t_status <- status[tail]
    h_active <- active[head] == 1
    t_active <- active[tail] == 1

    # Discordant edges (one S, one I, both active)
    sus_in_h <- h_status == "s" & t_status == "i" & h_active & t_active
    sus_in_t <- t_status == "s" & h_status == "i" & h_active & t_active
    disc <- sus_in_h | sus_in_t
    if (!any(disc, na.rm = TRUE)) next

    sus <- ifelse(sus_in_h[disc], head[disc], tail[disc])
    inf <- ifelse(sus_in_h[disc], tail[disc], head[disc])

    base.p <- layer_probs[k]
    stages <- inf_stage[inf]
    trans.p <- ifelse(stages %in% "ia", base.p * asymp.mult, base.p)

    # Susceptible-side vaccine effect (leaky reduction)
    vax <- vax_status[sus]
    eff <- rep(0, length(sus))
    eff[!is.na(vax) & vax == "elderly_vax"] <- vax.eff.elderly
    eff[!is.na(vax) & vax == "infant_proph"] <- vax.eff.infant
    trans.p <- trans.p * (1 - eff)

    new_inf <- sus[rbinom(length(trans.p), 1, trans.p) == 1]
    if (length(new_inf) > 0) all_new <- c(all_new, new_inf)
  }

  all_new <- unique(all_new)
  if (length(all_new) > 0) {
    status[all_new] <- "e"
    infTime[all_new] <- at
    dat <- set_attr(dat, "status", status)
    dat <- set_attr(dat, "infTime", infTime)
  }

  dat <- set_epi(dat, "se.flow", at, length(all_new))
  return(dat)
}


# Progression module -------------------------------------------------------

progress <- function(dat, at) {
  # E -> I(p) -> I(s) or I(a) -> R progression with age-independent rates.
  # Severity is age-dependent (modeled post-hoc through hospitalization
  # rates), but the progression timeline itself does not vary by age.
  #
  # Two safeguards keep the modeled stage durations consistent with the
  # parameter table:
  #
  # 1. Transition candidates are taken from a SNAPSHOT of status and
  #    inf_stage at the start of progress(). Updates are written back to
  #    the working vectors but not visible to subsequent transition
  #    blocks at this timestep -- so a node cannot cascade through
  #    multiple stages (e.g. E -> Ip -> Is) within a single step.
  #
  # 2. The E -> Ip transition additionally requires `infTime < at`,
  #    excluding new infections written by infect() earlier in the same
  #    timestep. This enforces a minimum of one timestep in the latent
  #    compartment regardless of `ei.rate`.

  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  inf_stage <- get_attr(dat, "inf_stage")
  age <- get_attr(dat, "age")
  infTime <- get_attr(dat, "infTime")

  ei.rate <- get_param(dat, "ei.rate")
  ip.rate <- get_param(dat, "ip.rate")
  ir.rate <- get_param(dat, "ir.rate")
  asymp.prob <- get_param(dat, "asymp.prob")

  # Snapshot the entry state. All eligibility checks below read from
  # these snapshots so transitions are determined before any updates.
  status0 <- status
  inf_stage0 <- inf_stage

  ## E -> I(presymp) ##
  # Exclude brand-new infections (infTime == at) so the latent period is
  # always at least one timestep.
  ids_e <- which(active == 1 & status0 == "e" &
                 (!is.na(infTime) & infTime < at))
  n_ei <- 0
  if (length(ids_e) > 0) {
    hit <- which(rbinom(length(ids_e), 1, ei.rate) == 1)
    new_ip <- ids_e[hit]
    n_ei <- length(new_ip)
    if (n_ei > 0) {
      status[new_ip] <- "i"
      inf_stage[new_ip] <- "ip"
    }
  }

  ## I(p) -> I(s) or I(a) ##
  ids_ip <- which(active == 1 & status0 == "i" &
                  !is.na(inf_stage0) & inf_stage0 == "ip")
  n_ips <- 0; n_ipa <- 0
  if (length(ids_ip) > 0) {
    hit <- which(rbinom(length(ids_ip), 1, ip.rate) == 1)
    new_clin <- ids_ip[hit]
    if (length(new_clin) > 0) {
      asymp <- rbinom(length(new_clin), 1, asymp.prob) == 1
      inf_stage[new_clin[asymp]] <- "ia"
      inf_stage[new_clin[!asymp]] <- "is"
      n_ipa <- sum(asymp)
      n_ips <- sum(!asymp)
    }
  }

  ## I(s/a) -> R ##
  ids_inf <- which(active == 1 & status0 == "i" &
                   !is.na(inf_stage0) & inf_stage0 %in% c("is", "ia"))
  n_ir <- 0
  if (length(ids_inf) > 0) {
    hit <- which(rbinom(length(ids_inf), 1, ir.rate) == 1)
    new_r <- ids_inf[hit]
    n_ir <- length(new_r)
    if (n_ir > 0) {
      status[new_r] <- "r"
      inf_stage[new_r] <- NA_character_
    }
  }

  dat <- set_attr(dat, "status", status)
  dat <- set_attr(dat, "inf_stage", inf_stage)

  dat <- set_epi(dat, "ei.flow", at, n_ei)
  dat <- set_epi(dat, "ips.flow", at, n_ips)
  dat <- set_epi(dat, "ipa.flow", at, n_ipa)
  dat <- set_epi(dat, "ir.flow", at, n_ir)

  ## Compartment counts ##
  is_active <- active == 1
  dat <- set_epi(dat, "e.num", at, sum(is_active & status == "e"))
  dat <- set_epi(dat, "ip.num", at,
                 sum(is_active & status == "i" & inf_stage %in% "ip"))
  dat <- set_epi(dat, "is.num", at,
                 sum(is_active & status == "i" & inf_stage %in% "is"))
  dat <- set_epi(dat, "ia.num", at,
                 sum(is_active & status == "i" & inf_stage %in% "ia"))
  dat <- set_epi(dat, "r.num", at, sum(is_active & status == "r"))

  ## Age-stratified cumulative infection counts ##
  for (a in c("infant", "young", "school", "adult", "elderly")) {
    in_a <- is_active & age == a
    dat <- set_epi(dat, paste0("cuminf.", a), at,
                   sum(in_a & status %in% c("e", "i", "r")))
  }

  return(dat)
}
