##
## SEIR with Contact Tracing for an Acute, Immunizing Infection
## EpiModel Gallery (https://github.com/EpiModel/EpiModel-Gallery)
##
## Author: Samuel M. Jenness (Emory University)
## Date: September 2026
##


# Attribute initializer ----------------------------------------------------

init_attrs <- function(dat, at) {
  # One-shot setup on the first module pass. EpiModel manages the built-in
  # attributes `status` (s, e, i, r), `active`, and `infTime`. Six
  # attributes are added for the natural history and the interventions:
  #
  #   inf.stage   -- substage of I: "ip" presymptomatic, "is" symptomatic,
  #                  "ia" asymptomatic; NA outside I
  #   symp.time   -- step of symptom onset (NA before onset or if asymptomatic)
  #   dx.due      -- step on which a symptomatic case will be diagnosed
  #                  (NA if the case never seeks a test)
  #   dx.time     -- step of diagnosis (NA if never diagnosed)
  #   iso.until   -- last step of a diagnosed index's isolation (NA if none)
  #   quar.until  -- last step of a traced contact's quarantine (NA if none)
  if (is.null(get_attr(dat, "inf.stage", override.null.error = TRUE))) {
    status <- get_attr(dat, "status")
    n <- length(status)
    asymp.prob <- get_param(dat, "asymp.prob")

    # Seeds from init.net() enter as status "i" and are split into the
    # presymptomatic and asymptomatic substages like every later infection.
    inf.stage <- rep(NA_character_, n)
    seeds <- which(status == "i")
    asymp <- rbinom(length(seeds), 1, asymp.prob) == 1
    inf.stage[seeds] <- ifelse(asymp, "ia", "ip")

    dat <- set_attr(dat, "inf.stage", inf.stage)
    for (nm in c("symp.time", "dx.due", "dx.time", "iso.until", "quar.until")) {
      dat <- set_attr(dat, nm, rep(NA_integer_, n))
    }
  }
  return(dat)
}


# Infection module ---------------------------------------------------------

infect <- function(dat, at) {
  # S -> E across the discordant edges of the current network. The per-edge
  # daily transmission probability depends on the infector's substage, and
  # is multiplied by iso.mult when the infector is a diagnosed index in
  # isolation, or by quar.mult when either partner is a traced contact in
  # quarantine. Every transmission is recorded with set_transmat() so the
  # analysis can attribute infections to substages and restriction status.
  # The numbers of people isolated and quarantined on this day's contacts
  # are recorded here, where the restrictions are applied.
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  infTime <- get_attr(dat, "infTime")
  inf.stage <- get_attr(dat, "inf.stage")
  iso.until <- get_attr(dat, "iso.until")
  quar.until <- get_attr(dat, "quar.until")

  inf.prob <- get_param(dat, "inf.prob")
  is.inf.mult <- get_param(dat, "is.inf.mult")
  ia.inf.mult <- get_param(dat, "ia.inf.mult")
  iso.mult <- get_param(dat, "iso.mult")
  quar.mult <- get_param(dat, "quar.mult")

  isolated <- !is.na(iso.until) & at <= iso.until
  quarantined <- !is.na(quar.until) & at <= quar.until

  del <- discord_edgelist(dat, at, network = 1, infstat = "i")

  nInf <- 0
  nInf.stage <- c(ip = 0, is = 0, ia = 0)
  if (!is.null(del) && nrow(del) > 0) {
    stage <- inf.stage[del$inf]
    p <- inf.prob * ifelse(stage == "is", is.inf.mult,
                           ifelse(stage == "ia", ia.inf.mult, 1))
    # The stronger of the two restrictions applies when both are present.
    restrict <- ifelse(isolated[del$inf], iso.mult, 1)
    restrict <- pmin(restrict,
                     ifelse(quarantined[del$inf] | quarantined[del$sus],
                            quar.mult, 1))
    p <- p * restrict
    hit <- which(rbinom(length(p), 1, p) == 1)

    if (length(hit) > 0) {
      # A susceptible exposed by several partners in one step is infected
      # once, by a random one of them: discord_edgelist() returns the edges
      # in random order, so keeping the first row per susceptible does it.
      del <- del[hit, , drop = FALSE]
      del <- del[!duplicated(del$sus), , drop = FALSE]
      newInf <- del$sus
      nInf <- length(newInf)
      status[newInf] <- "e"
      infTime[newInf] <- at
      dat <- set_attr(dat, "status", status)
      dat <- set_attr(dat, "infTime", infTime)

      # Transmission record: the infector's substage and restriction
      # status ride along with the standard columns.
      del$infStage <- inf.stage[del$inf]
      del$infIsolated <- as.integer(isolated[del$inf])
      del$anyQuarantined <- as.integer(quarantined[del$inf] |
                                         quarantined[del$sus])
      del$infTime <- infTime[del$inf]
      dat <- set_transmat(dat, del, at)

      tab <- table(factor(del$infStage, levels = c("ip", "is", "ia")))
      nInf.stage[] <- as.numeric(tab)
    }
  }

  dat <- set_epi(dat, "se.flow", at, nInf)
  dat <- set_epi(dat, "se.flow.ip", at, nInf.stage[["ip"]])
  dat <- set_epi(dat, "se.flow.is", at, nInf.stage[["is"]])
  dat <- set_epi(dat, "se.flow.ia", at, nInf.stage[["ia"]])
  dat <- set_epi(dat, "iso.num", at, sum(active == 1 & isolated))
  dat <- set_epi(dat, "quar.num", at,
                 sum(active == 1 & quarantined & !isolated))
  return(dat)
}


# Progression and diagnosis module ----------------------------------------

progress <- function(dat, at) {
  # Stage transitions, all drawn from a snapshot of the state at entry so a
  # node moves at most one stage per step:
  #   E  -> Ip (prob 1 - asymp.prob) or Ia (prob asymp.prob) at ei.rate
  #   Ip -> Is at ips.rate (symptom onset)
  #   Is -> R  at isr.rate
  #   Ia -> R  at iar.rate
  # Diagnosis is scheduled at symptom onset: with probability dx.prob the
  # case will be diagnosed after a delay of 1 + Geometric(1 / dx.delay)
  # days, so the mean delay from onset to diagnosis is dx.delay days.
  # Diagnosis starts isolation, which applies to the next iso.duration
  # days of contacts. Asymptomatic infections are never diagnosed through
  # this pathway.
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  inf.stage <- get_attr(dat, "inf.stage")
  infTime <- get_attr(dat, "infTime")
  symp.time <- get_attr(dat, "symp.time")
  dx.due <- get_attr(dat, "dx.due")
  dx.time <- get_attr(dat, "dx.time")
  iso.until <- get_attr(dat, "iso.until")

  ei.rate <- get_param(dat, "ei.rate")
  asymp.prob <- get_param(dat, "asymp.prob")
  ips.rate <- get_param(dat, "ips.rate")
  isr.rate <- get_param(dat, "isr.rate")
  iar.rate <- get_param(dat, "iar.rate")
  dx.prob <- get_param(dat, "dx.prob")
  dx.delay <- get_param(dat, "dx.delay")
  iso.duration <- get_param(dat, "iso.duration")

  status0 <- status
  stage0 <- inf.stage

  ## E -> Ip or Ia (infTime < at guarantees at least one step in E)
  ids_e <- which(active == 1 & status0 == "e" & infTime < at)
  n_ei <- 0
  if (length(ids_e) > 0) {
    new_i <- ids_e[rbinom(length(ids_e), 1, ei.rate) == 1]
    n_ei <- length(new_i)
    if (n_ei > 0) {
      asymp <- rbinom(n_ei, 1, asymp.prob) == 1
      status[new_i] <- "i"
      inf.stage[new_i] <- ifelse(asymp, "ia", "ip")
    }
  }

  ## Ip -> Is: symptom onset, and the diagnosis draw
  ids_ip <- which(active == 1 & status0 == "i" & stage0 %in% "ip")
  n_ips <- 0
  if (length(ids_ip) > 0) {
    new_is <- ids_ip[rbinom(length(ids_ip), 1, ips.rate) == 1]
    n_ips <- length(new_is)
    if (n_ips > 0) {
      inf.stage[new_is] <- "is"
      symp.time[new_is] <- at
      seek <- rbinom(n_ips, 1, dx.prob) == 1
      dx.due[new_is[seek]] <- at + 1 + rgeom(sum(seek), 1 / dx.delay)
    }
  }

  ## Is -> R and Ia -> R
  ids_is <- which(active == 1 & status0 == "i" & stage0 %in% "is")
  ids_ia <- which(active == 1 & status0 == "i" & stage0 %in% "ia")
  new_r <- c(ids_is[rbinom(length(ids_is), 1, isr.rate) == 1],
             ids_ia[rbinom(length(ids_ia), 1, iar.rate) == 1])
  n_ir <- length(new_r)
  if (n_ir > 0) {
    status[new_r] <- "r"
    inf.stage[new_r] <- NA_character_
  }

  ## Diagnosis of the cases scheduled for today. A case that has already
  ## recovered is still diagnosed (a test detects infection, not
  ## infectiousness): its isolation has no effect on transmission, but it
  ## still triggers tracing of its contacts.
  new_dx <- which(active == 1 & !is.na(dx.due) & dx.due <= at &
                  is.na(dx.time))
  n_dx <- length(new_dx)
  if (n_dx > 0) {
    dx.time[new_dx] <- at
    iso.until[new_dx] <- at + iso.duration
  }

  dat <- set_attr(dat, "status", status)
  dat <- set_attr(dat, "inf.stage", inf.stage)
  dat <- set_attr(dat, "symp.time", symp.time)
  dat <- set_attr(dat, "dx.due", dx.due)
  dat <- set_attr(dat, "dx.time", dx.time)
  dat <- set_attr(dat, "iso.until", iso.until)

  is_active <- active == 1
  dat <- set_epi(dat, "ei.flow", at, n_ei)
  dat <- set_epi(dat, "ips.flow", at, n_ips)
  dat <- set_epi(dat, "ir.flow", at, n_ir)
  dat <- set_epi(dat, "dx.flow", at, n_dx)
  dat <- set_epi(dat, "e.num", at, sum(is_active & status == "e"))
  dat <- set_epi(dat, "ip.num", at, sum(is_active & inf.stage %in% "ip"))
  dat <- set_epi(dat, "is.num", at, sum(is_active & inf.stage %in% "is"))
  dat <- set_epi(dat, "ia.num", at, sum(is_active & inf.stage %in% "ia"))
  dat <- set_epi(dat, "r.num", at, sum(is_active & status == "r"))
  return(dat)
}


# Contact tracing module ---------------------------------------------------

trace <- function(dat, at) {
  # For each index diagnosed trace.delay steps ago, look up the partners
  # recorded in the cumulative edgelist, keep the partnerships that overlap
  # the index's contact elicitation window (from trace.window days before
  # symptom onset to the day of diagnosis), reach each partner with
  # probability trace.reach.prob, and place the reached contacts in
  # quarantine for quar.duration days.
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  inf.stage <- get_attr(dat, "inf.stage")
  symp.time <- get_attr(dat, "symp.time")
  dx.time <- get_attr(dat, "dx.time")
  quar.until <- get_attr(dat, "quar.until")

  trace.reach.prob <- get_param(dat, "trace.reach.prob")
  trace.delay <- get_param(dat, "trace.delay")
  trace.window <- get_param(dat, "trace.window")
  quar.duration <- get_param(dat, "quar.duration")

  n_index <- 0
  n_part <- 0
  n_part_ended <- 0
  n_reach <- 0
  n_quar_start <- 0
  reach_state <- c(s = 0, e = 0, ip = 0, is = 0, ia = 0, r = 0)

  if (trace.reach.prob > 0) {
    # 1. Indices whose trace is due today
    idsIndex <- which(active == 1 & !is.na(dx.time) &
                      (at - dx.time) == trace.delay)
    n_index <- length(idsIndex)

    if (n_index > 0) {
      # 2. Their partners from the cumulative edgelist. get_partners()
      #    takes positional ids and returns one row per partnership with
      #    the index and partner as unique ids and the partnership's start
      #    and stop steps (stop is NA while the partnership is active).
      part_df <- get_partners(dat, idsIndex, only.active.nodes = TRUE)

      if (!is.null(part_df) && nrow(part_df) > 0) {
        # 3. Keep the partnerships that overlap each index's elicitation
        #    window: still active or ended no earlier than trace.window
        #    days before the index's symptom onset, and begun no later
        #    than the index's diagnosis.
        index_pid <- get_posit_ids(dat, part_df$index)
        window_start <- symp.time[index_pid] - trace.window
        window_end <- dx.time[index_pid]
        in_window <- (is.na(part_df$stop) | part_df$stop >= window_start) &
                     part_df$start <= window_end
        part_df <- part_df[in_window, , drop = FALSE]

        # 4. Partner unique ids back to positional ids, one row per
        #    partner, dropping partners who are already diagnosed (their
        #    isolation is already in force).
        partner_pid <- get_posit_ids(dat, part_df$partner)
        ended <- !is.na(part_df$stop)
        keep <- !duplicated(partner_pid) & is.na(dx.time[partner_pid])
        partner_pid <- partner_pid[keep]
        ended <- ended[keep]
        n_part <- length(partner_pid)
        n_part_ended <- sum(ended)

        if (n_part > 0) {
          # 5. Reach each partner with probability trace.reach.prob and
          #    quarantine the reached contacts, extending any quarantine
          #    already in force rather than shortening it.
          reached <- partner_pid[rbinom(n_part, 1, trace.reach.prob) == 1]
          n_reach <- length(reached)

          if (n_reach > 0) {
            state <- ifelse(status[reached] == "i", inf.stage[reached],
                            status[reached])
            reach_state[] <- as.numeric(table(factor(state,
                                                     levels = names(reach_state))))
            in_quar <- !is.na(quar.until[reached]) & at <= quar.until[reached]
            n_quar_start <- sum(!in_quar)
            quar.until[reached] <- pmax(quar.until[reached],
                                        at + quar.duration, na.rm = TRUE)
          }
        }
      }
    }
  }

  dat <- set_attr(dat, "quar.until", quar.until)

  dat <- set_epi(dat, "trace.index.flow", at, n_index)
  dat <- set_epi(dat, "trace.part.flow", at, n_part)
  dat <- set_epi(dat, "trace.part.ended.flow", at, n_part_ended)
  dat <- set_epi(dat, "trace.reach.flow", at, n_reach)
  dat <- set_epi(dat, "quar.start.flow", at, n_quar_start)
  dat <- set_epi(dat, "reach.s.flow", at, reach_state[["s"]])
  dat <- set_epi(dat, "reach.e.flow", at, reach_state[["e"]])
  dat <- set_epi(dat, "reach.i.flow", at,
                 reach_state[["ip"]] + reach_state[["is"]] + reach_state[["ia"]])
  dat <- set_epi(dat, "reach.r.flow", at, reach_state[["r"]])
  return(dat)
}
