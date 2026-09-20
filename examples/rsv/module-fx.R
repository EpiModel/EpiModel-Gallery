##
## RSV: Age-Stratified SEIR Over a Multilayer Network
## EpiModel Gallery (https://github.com/EpiModel/EpiModel-Gallery)
##
## Authors: Samuel M. Jenness (Emory University)
## Date: September 2026
##


# Attribute initializer ----------------------------------------------------

init_attrs <- function(dat, at) {
  # One-shot setup of the nodal attributes used by the disease and
  # intervention modules. Runs on the first call only; EpiModel manages
  # the `active`, `status`, and `infTime` attributes itself.
  #
  # vax_status:
  #   NA             -- not immunized
  #   "elderly_vax"  -- received the older-adult vaccine
  #   "infant_proph" -- received the infant monoclonal antibody
  #   "cocoon"       -- adult co-resident of an infant, given a hypothetical
  #                     transmission-blocking product (household targeting)
  #
  # inf_stage (I substage; NA when status is "s", "e", or "r"):
  #   "ip" -- presymptomatic infectious
  #   "is" -- symptomatic infectious
  #   "ia" -- asymptomatic infectious

  if (is.null(get_attr(dat, "vax_status", override.null.error = TRUE))) {
    active <- get_attr(dat, "active")
    status <- get_attr(dat, "status")
    age <- get_attr(dat, "age")
    hh_id <- get_attr(dat, "hh_id")

    elderly.cov <- get_param(dat, "elderly.vax.coverage")
    infant.cov <- get_param(dat, "infant.proph.coverage")
    cocoon.cov <- get_param(dat, "cocoon.coverage")

    vax_status <- rep(NA_character_, length(active))
    inf_stage <- rep(NA_character_, length(active))

    # Seed infections from init.net() start in the presymptomatic stage so
    # that progress() moves them through the I substages to recovery.
    # Without this they would stay status = "i" with inf_stage = NA and
    # act as a permanent infectious reservoir.
    inf_stage[active == 1 & status == "i"] <- "ip"

    # Immunization is delivered before the season to a random fraction of
    # each eligible age group, independent of infection status.
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
    # Household targeting: the explicit household ids make "adults who live
    # with an infant" a definable target group.
    infant_hh <- unique(hh_id[active == 1 & age == "infant"])
    elig_c <- which(active == 1 & age == "adult" & hh_id %in% infant_hh)
    if (length(elig_c) > 0 && cocoon.cov > 0) {
      hit <- which(rbinom(length(elig_c), 1, cocoon.cov) == 1)
      vax_status[elig_c[hit]] <- "cocoon"
    }

    dat <- set_attr(dat, "vax_status", vax_status)
    dat <- set_attr(dat, "inf_stage", inf_stage)
  }

  return(dat)
}


# Infection module ---------------------------------------------------------

infect <- function(dat, at) {
  # S -> E transmission over two layers:
  #   1 = household: a fixed edgelist of within-household pairs (every
  #       household is a clique), passed in as the parameter hh.pairs. It is
  #       not an ERGM and is never resimulated, so it lives outside netsim's
  #       network machinery and is walked here directly.
  #   2 = community: the TERGM layer simulated by netsim, read each step with
  #       get_edgelist().
  # Walking each layer separately (rather than calling discord_edgelist)
  # lets the layers carry different per-contact transmission probabilities
  # and lets the NPI act on the community layer only.
  #
  # The per-edge, per-day transmission probability is the layer's inf.prob
  # multiplied by:
  #   - asymp.inf.mult if the infectious partner is asymptomatic
  #   - (1 - npi.mask.efficacy) on the community layer while the NPI is on
  #   - sus.mult[age] of the susceptible partner (prior-immunity proxy)
  #   - (1 - eff.inf) if the susceptible partner is immunized

  ## Attributes ##
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  age <- get_attr(dat, "age")
  inf_stage <- get_attr(dat, "inf_stage")
  vax_status <- get_attr(dat, "vax_status")
  infTime <- get_attr(dat, "infTime")

  ## Parameters ##
  ip.fam <- get_param(dat, "inf.prob.household")
  ip.com <- get_param(dat, "inf.prob.community")
  asymp.mult <- get_param(dat, "asymp.inf.mult")
  sus.mult <- get_param(dat, "sus.mult")

  npi.start <- get_param(dat, "npi.start")
  npi.end <- get_param(dat, "npi.end")
  npi.mask.eff <- get_param(dat, "npi.mask.efficacy")
  npi.contact.mult <- get_param(dat, "npi.contact.mult")

  eff.inf.elderly <- get_param(dat, "elderly.vax.eff.inf")
  eff.inf.infant <- get_param(dat, "infant.proph.eff.inf")
  eff.inf.cocoon <- get_param(dat, "cocoon.eff.inf")
  hh_el <- get_param(dat, "hh.pairs")

  npi.on <- (at >= npi.start && at <= npi.end)
  com.prob.mult <- if (npi.on) (1 - npi.mask.eff) else 1
  # npi.contact.mult reduces the community contact RATE; we approximate
  # this by keeping each community edge with that probability.
  com.contact.mult <- if (npi.on) npi.contact.mult else 1

  layer_probs <- c(ip.fam, ip.com * com.prob.mult)
  all_new <- integer(0)
  n_new <- c(0, 0)          # new infections by layer (household, community)
  n_new_infant <- c(0, 0)   # the same, infants only

  for (k in 1:2) {
    el <- if (k == 1) hh_el else get_edgelist(dat, network = 1)
    if (is.null(el) || nrow(el) == 0) next

    head <- el[, 1]
    tail <- el[, 2]

    # Optional community-edge thinning under the NPI
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

    # Infectious-side modifiers
    base.p <- layer_probs[k]
    stages <- inf_stage[inf]
    trans.p <- ifelse(stages %in% "ia", base.p * asymp.mult, base.p)

    # Susceptible-side modifiers: age-specific susceptibility, then the
    # leaky infection-blocking component of immunization
    trans.p <- trans.p * as.numeric(sus.mult[age[sus]])
    vax <- vax_status[sus]
    eff <- rep(0, length(sus))
    eff[!is.na(vax) & vax == "elderly_vax"] <- eff.inf.elderly
    eff[!is.na(vax) & vax == "infant_proph"] <- eff.inf.infant
    eff[!is.na(vax) & vax == "cocoon"] <- eff.inf.cocoon
    trans.p <- trans.p * (1 - eff)

    new_inf <- sus[rbinom(length(trans.p), 1, trans.p) == 1]
    # A node exposed on both layers in the same step is attributed to the
    # layer walked first (household), so the per-layer counts sum exactly.
    new_inf <- setdiff(unique(new_inf), all_new)
    n_new[k] <- length(new_inf)
    n_new_infant[k] <- sum(age[new_inf] == "infant")
    if (length(new_inf) > 0) all_new <- c(all_new, new_inf)
  }

  if (length(all_new) > 0) {
    status[all_new] <- "e"
    infTime[all_new] <- at
    dat <- set_attr(dat, "status", status)
    dat <- set_attr(dat, "infTime", infTime)
  }

  dat <- set_epi(dat, "se.flow", at, length(all_new))
  dat <- set_epi(dat, "se.flow.hh", at, n_new[1])
  dat <- set_epi(dat, "se.flow.com", at, n_new[2])
  dat <- set_epi(dat, "se.flow.infant.hh", at, n_new_infant[1])
  dat <- set_epi(dat, "se.flow.infant.com", at, n_new_infant[2])
  return(dat)
}


# Progression module -------------------------------------------------------

progress <- function(dat, at) {
  # E -> I(p) -> I(s) or I(a) -> R with age-independent rates. Severity is
  # age-dependent but handled post hoc through per-infection hospitalization
  # risks, so the progression timeline itself does not vary by age.
  #
  # Two safeguards keep the modeled stage durations consistent with the
  # parameter table:
  #
  # 1. Transition candidates are taken from a snapshot of status and
  #    inf_stage at the start of progress(), so a node cannot cascade
  #    through several stages (e.g. E -> Ip -> Is) within one step.
  #
  # 2. E -> Ip additionally requires infTime < at. Because infection runs
  #    before progression within each step (see module.order in
  #    control.net), this excludes the infections written by infect() in
  #    the same step and guarantees at least one step in E.

  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  inf_stage <- get_attr(dat, "inf_stage")
  age <- get_attr(dat, "age")
  vax_status <- get_attr(dat, "vax_status")
  infTime <- get_attr(dat, "infTime")

  ei.rate <- get_param(dat, "ei.rate")
  ip.rate <- get_param(dat, "ip.rate")
  ir.rate <- get_param(dat, "ir.rate")
  asymp.prob <- get_param(dat, "asymp.prob")

  status0 <- status
  inf_stage0 <- inf_stage

  ## E -> I(presymp) ##
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

  ## Cumulative incident infections by age ##
  # Seed infections from init.net() carry infTime = 1 and are excluded, so
  # these counters track infections acquired during the simulated season.
  # The immunized subsets (".prot") are tracked separately so that the
  # analysis can apply the severity-reducing component of each product
  # only to infections among immunized people.
  incident <- is_active & !is.na(infTime) & infTime > 1
  for (a in c("infant", "young", "school", "adult", "elderly")) {
    dat <- set_epi(dat, paste0("cuminf.", a), at, sum(incident & age == a))
  }
  inf_prot <- is_active & vax_status %in% "infant_proph"
  eld_prot <- is_active & vax_status %in% "elderly_vax"
  dat <- set_epi(dat, "cuminf.infant.prot", at, sum(incident & inf_prot))
  dat <- set_epi(dat, "cuminf.elderly.prot", at, sum(incident & eld_prot))
  dat <- set_epi(dat, "n.infant.prot", at, sum(inf_prot))
  dat <- set_epi(dat, "n.elderly.prot", at, sum(eld_prot))
  dat <- set_epi(dat, "n.cocoon", at, sum(is_active & vax_status %in% "cocoon"))

  return(dat)
}
