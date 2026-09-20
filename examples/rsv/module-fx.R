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
  #   "cocoon"       -- given a hypothetical transmission-blocking product,
  #                     targeted either at the co-residents of infants
  #                     (household cocooning) or, as a comparator with the
  #                     same number of doses and the same age mix, at people
  #                     chosen at random (cocoon.random = 1)
  #
  # inf_stage (I substage; NA when status is "s", "e", or "r"):
  #   "ip" -- presymptomatic infectious (infections that will be symptomatic)
  #   "is" -- symptomatic infectious
  #   "ia" -- asymptomatic infectious (never symptomatic)

  if (is.null(get_attr(dat, "vax_status", override.null.error = TRUE))) {
    active <- get_attr(dat, "active")
    status <- get_attr(dat, "status")
    age <- get_attr(dat, "age")
    hh_id <- get_attr(dat, "hh_id")

    elderly.cov <- get_param(dat, "elderly.vax.coverage")
    infant.cov <- get_param(dat, "infant.proph.coverage")
    cocoon.cov <- get_param(dat, "cocoon.coverage")
    cocoon.random <- get_param(dat, "cocoon.random")
    asymp.prob <- get_param(dat, "asymp.prob")

    vax_status <- rep(NA_character_, length(active))
    inf_stage <- rep(NA_character_, length(active))

    # Seed infections from init.net() are placed in an infectious substage
    # so that progress() moves them through to recovery. Without this they
    # would stay status = "i" with inf_stage = NA and act as a permanent
    # infectious reservoir. Seeds are split into symptomatic-track and
    # asymptomatic infections with the same probability as later cases.
    seeds <- which(active == 1 & status == "i")
    if (length(seeds) > 0) {
      asymp <- rbinom(length(seeds), 1, asymp.prob) == 1
      inf_stage[seeds[asymp]] <- "ia"
      inf_stage[seeds[!asymp]] <- "ip"
    }

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
    # Household targeting: the explicit household ids make "people who live
    # with an infant" a definable target group; every co-resident of an
    # infant, of any age, is eligible. The random comparator gives the same
    # number of doses to people of the same ages drawn from the whole
    # population, which isolates what the household link contributes.
    infant_hh <- unique(hh_id[active == 1 & age == "infant"])
    elig_c <- which(active == 1 & age != "infant" & hh_id %in% infant_hh)
    if (length(elig_c) > 0 && cocoon.cov > 0) {
      hit <- elig_c[rbinom(length(elig_c), 1, cocoon.cov) == 1]
      if (cocoon.random == 1 && length(hit) > 0) {
        n_by_age <- table(age[hit])
        hit <- unlist(lapply(names(n_by_age), function(a) {
          pool <- which(active == 1 & age == a)
          pool[sample.int(length(pool), n_by_age[[a]])]
        }))
      }
      vax_status[hit] <- "cocoon"
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
  #
  # Every discordant edge on both layers is a separate Bernoulli trial. A
  # susceptible node with more than one successful exposure in the same
  # step is infected once, and its infector (and hence its layer) is chosen
  # at random among the successful exposures. This is the tie-breaking rule
  # EpiModel's built-in infection module uses, and it keeps the attribution
  # of infections to layers unbiased. Each transmission is recorded with
  # set_transmat() so that who-infected-whom can be analyzed afterwards.

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
  seas.amp <- get_param(dat, "seas.amp")
  seas.peak <- get_param(dat, "seas.peak")
  hh_el <- get_param(dat, "hh.pairs")

  npi.on <- (at >= npi.start && at <= npi.end)
  com.prob.mult <- if (npi.on) (1 - npi.mask.eff) else 1
  # npi.contact.mult reduces the community contact RATE; we approximate
  # this by keeping each community edge with that probability.
  com.contact.mult <- if (npi.on) npi.contact.mult else 1

  # Seasonal forcing: a cosine multiplier on both layers' transmission
  # probabilities with a period of one year, equal to 1 + seas.amp on day
  # seas.peak. RSV seasons end because transmissibility falls, not only
  # because susceptibles run out.
  seas.mult <- 1 + seas.amp * cos(2 * pi * (at - seas.peak) / 365)

  layer_probs <- c(ip.fam, ip.com * com.prob.mult) * seas.mult
  del <- NULL   # one row per successful exposure: sus, inf, layer

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

    hit <- which(rbinom(length(trans.p), 1, trans.p) == 1)
    if (length(hit) > 0) {
      del <- rbind(del, data.frame(sus = sus[hit], inf = inf[hit], layer = k))
    }
  }

  n_new <- c(0, 0)          # new infections by layer (household, community)
  n_new_infant <- c(0, 0)   # the same, infants only

  if (!is.null(del) && nrow(del) > 0) {
    # Competing exposures: shuffle the successful exposures, then keep one
    # per susceptible node. Shuffling first makes the kept row a uniform
    # random draw among that node's successes, whatever layer they came from.
    del <- del[sample.int(nrow(del)), , drop = FALSE]
    del <- del[!duplicated(del$sus), , drop = FALSE]

    new_inf <- del$sus
    status[new_inf] <- "e"
    infTime[new_inf] <- at
    dat <- set_attr(dat, "status", status)
    dat <- set_attr(dat, "infTime", infTime)

    # Transmission record: who infected whom, on which layer, with both
    # ages and the infector's own infection time (1 for the seeds).
    del$at <- at
    del$susAge <- age[del$sus]
    del$infAge <- age[del$inf]
    del$infTime <- infTime[del$inf]
    dat <- set_transmat(dat, del, at)

    n_new <- c(sum(del$layer == 1), sum(del$layer == 2))
    is_inf <- del$susAge == "infant"
    n_new_infant <- c(sum(is_inf & del$layer == 1), sum(is_inf & del$layer == 2))
  }

  dat <- set_epi(dat, "se.flow", at, sum(n_new))
  dat <- set_epi(dat, "se.flow.hh", at, n_new[1])
  dat <- set_epi(dat, "se.flow.com", at, n_new[2])
  dat <- set_epi(dat, "se.flow.infant.hh", at, n_new_infant[1])
  dat <- set_epi(dat, "se.flow.infant.com", at, n_new_infant[2])
  return(dat)
}


# Progression module -------------------------------------------------------

progress <- function(dat, at) {
  # E -> I(p) -> I(s) -> R for infections that will become symptomatic, and
  # E -> I(a) -> R for those that never will. Rates are age-independent and
  # each stage duration is geometric with the stated mean. Severity is
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
  # 2. E -> I additionally requires infTime < at. Because infection runs
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

  ## E -> I(presymp) or I(asymp) ##
  ids_e <- which(active == 1 & status0 == "e" &
                 (!is.na(infTime) & infTime < at))
  n_eip <- 0; n_eia <- 0
  if (length(ids_e) > 0) {
    hit <- which(rbinom(length(ids_e), 1, ei.rate) == 1)
    new_i <- ids_e[hit]
    if (length(new_i) > 0) {
      asymp <- rbinom(length(new_i), 1, asymp.prob) == 1
      status[new_i] <- "i"
      inf_stage[new_i[asymp]] <- "ia"
      inf_stage[new_i[!asymp]] <- "ip"
      n_eia <- sum(asymp)
      n_eip <- sum(!asymp)
    }
  }

  ## I(p) -> I(s) ##
  ids_ip <- which(active == 1 & status0 == "i" &
                  !is.na(inf_stage0) & inf_stage0 == "ip")
  n_ips <- 0
  if (length(ids_ip) > 0) {
    hit <- which(rbinom(length(ids_ip), 1, ip.rate) == 1)
    new_is <- ids_ip[hit]
    n_ips <- length(new_is)
    if (n_ips > 0) {
      inf_stage[new_is] <- "is"
    }
  }

  ## I(s) or I(a) -> R ##
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

  dat <- set_epi(dat, "ei.flow", at, n_eip + n_eia)
  dat <- set_epi(dat, "eip.flow", at, n_eip)
  dat <- set_epi(dat, "eia.flow", at, n_eia)
  dat <- set_epi(dat, "ips.flow", at, n_ips)
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
  # only to infections among immunized people, and so that the realized
  # effectiveness among recipients can be computed.
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
