
##
## Norovirus in a Nursing Home: Person-to-Person and Environmental Transmission
## EpiModel Gallery Example
## Custom module functions
##
## Each time step is one 8-hour shift: step `at` is shift (at - 1) %% 3 + 1
## (1 day, 2 evening, 3 night) of day (at - 1) %/% 3 + 1. Every module is
## self-contained (no helper defined outside it), so the modules run unchanged
## on parallel workers.
##


# Shared calendar arithmetic -----------------------------------------------------
#
# Several modules need the shift, the day of the season, and the weekday. Each
# computes them inline:
#
#   shift <- (at - 1) %% 3 + 1          # 1 day, 2 evening, 3 night
#   day <- (at - 1) %/% 3 + 1
#   weekday <- (day - 1) %% 7 + 1


# Staffing: who is at work this shift, and who cares for whom --------------------

staffing <- function(dat, at) {
  # The shift clock. Each staff member belongs to a pool for one unit and one
  # shift (or to the float pool, on the day shift) and has two fixed days off
  # each week. A staff member is on duty when it is their shift and weekday and
  # they are not absent: absent means symptomatic and staying home (everyone
  # who does not work while ill), or excluded under an outbreak policy.
  #
  # Care is assigned by bed. Each bed belongs to one nursing assistant (CNA)
  # care group per shift: 4 groups of 8 beds on days, 3 of about 11 on
  # evenings, 2 of 16 on nights. Pool member k is the primary CNA for group k;
  # when a primary is off or absent, a relief member of the same pool covers
  # the group, and when no relief is left, the group is added to the load of
  # another CNA on duty on the unit. Licensed nurses split the unit's residents
  # between them. Assignments are stored as staff ids (an attribute), never as
  # positional ids, because positions shift when residents are discharged.

  shift <- (at - 1) %% 3 + 1
  weekday <- ((at - 1) %/% 3) %% 7 + 1

  active <- get_attr(dat, "active")
  type <- get_attr(dat, "type")
  role <- get_attr(dat, "role")
  unit <- get_attr(dat, "unit")
  bed <- get_attr(dat, "bed")
  status <- get_attr(dat, "status")
  staff_id <- get_attr(dat, "staff_id")
  staff_shift <- get_attr(dat, "staff_shift")
  pool_k <- get_attr(dat, "pool_k")
  off1 <- get_attr(dat, "off1")
  off2 <- get_attr(dat, "off2")
  presentee <- get_attr(dat, "presentee")
  excluded <- get_attr(dat, "excluded")

  groups <- get_param(dat, "cna.groups")[shift]
  beds.per.unit <- get_param(dat, "beds.per.unit")

  is_staff <- active == 1 & type == "staff"
  absent <- (status == "i" & presentee == 0) | excluded == 1
  scheduled <- is_staff & staff_shift == shift &
    weekday != off1 & weekday != off2
  on_duty <- as.numeric(scheduled & !absent)

  carer <- rep(NA_real_, length(active))
  nurse_of <- rep(NA_real_, length(active))
  res <- which(active == 1 & type == "res")

  for (u in unique(unit[res])) {
    res_u <- res[unit[res] == u]
    cnas <- which(on_duty == 1 & role == "cna" & unit == u)
    nurses <- which(on_duty == 1 & role == "nurse" & unit == u)
    cnas <- cnas[order(pool_k[cnas])]

    # Cover each care group: primary, then relief, then shared load.
    cover <- rep(NA_real_, groups)
    if (length(cnas) > 0) {
      prim <- match(seq_len(groups), pool_k[cnas])
      cover[!is.na(prim)] <- staff_id[cnas[prim[!is.na(prim)]]]
      relief <- cnas[pool_k[cnas] > groups]
      gaps <- which(is.na(cover))
      n_fill <- min(length(gaps), length(relief))
      if (n_fill > 0) cover[gaps[seq_len(n_fill)]] <- staff_id[relief[seq_len(n_fill)]]
      gaps <- which(is.na(cover))
      if (length(gaps) > 0) {
        cover[gaps] <- staff_id[cnas[(seq_along(gaps) - 1) %% length(cnas) + 1]]
      }
    } else if (length(nurses) > 0) {
      cover[] <- staff_id[nurses[1]]
    }
    grp <- ceiling(bed[res_u] / (beds.per.unit / groups))
    carer[res_u] <- cover[grp]

    if (length(nurses) > 0) {
      half <- ceiling(bed[res_u] / (beds.per.unit / length(nurses)))
      nurse_of[res_u] <- staff_id[nurses[half]]
    }
  }

  dat <- set_attr(dat, "on_duty", on_duty)
  dat <- set_attr(dat, "carer", carer)
  dat <- set_attr(dat, "nurse_of", nurse_of)
  dat <- set_epi(dat, "staff.absent", sum(scheduled & absent))
  return(dat)
}


# Infection: every route in one module ------------------------------------------

infect <- function(dat, at) {
  # Person-to-person transmission over six kinds of contact, plus pickup from
  # each unit's environmental reservoir, plus staff infected in the community.
  #
  # Each exposure k of a susceptible person j carries a hazard lambda_k for
  # this shift. For a contact with an infectious person i it is
  #
  #   lambda = beta * w * infness_i * mult
  #
  # where w is the contact weight for that kind of contact on this shift (in
  # hours of close contact, with hands-on care weighted up), infness_i is 1
  # for a symptomatic case (vomit.mult if the case vomits) and rel.inf.a for
  # someone shedding without symptoms (one value for residents, one for
  # staff), and mult carries any contact precaution or hand hygiene
  # multiplier. For the reservoir of unit u it is
  #
  #   lambda = w_env * E_u / N_u
  #
  # where E_u is the unit's contamination at the end of the previous shift,
  # N_u the number of people on the unit this shift, and w_env the person's
  # time on the unit's shared surfaces. Staff touch surfaces with their hands,
  # so soap and water during an outbreak also lowers staff pickup. Each
  # exposure draws an exponential
  # waiting time with its hazard; a susceptible is infected if the first
  # arrival falls inside the shift, and that exposure is the route and source.
  # This is the competing-risks form of 1 - exp(-sum(lambda)).

  shift <- (at - 1) %% 3 + 1
  day <- (at - 1) %/% 3 + 1

  active <- get_attr(dat, "active")
  type <- get_attr(dat, "type")
  role <- get_attr(dat, "role")
  unit <- get_attr(dat, "unit")
  status <- get_attr(dat, "status")
  vomit <- get_attr(dat, "vomit")
  isolated <- get_attr(dat, "isolated")
  on_duty <- get_attr(dat, "on_duty")
  staff_id <- get_attr(dat, "staff_id")
  carer <- get_attr(dat, "carer")
  nurse_of <- get_attr(dat, "nurse_of")
  float_unit <- get_attr(dat, "float_unit")

  beta <- get_param(dat, "beta")
  rel.inf.a <- get_param(dat, "rel.inf.a")
  vomit.mult <- get_param(dat, "vomit.mult")
  w.room <- get_param(dat, "w.room")[shift]
  w.common <- get_param(dat, "w.common")[shift]
  w.team <- get_param(dat, "w.team")[shift]
  w.care <- get_param(dat, "w.care")[shift]
  w.nurse <- get_param(dat, "w.nurse")[shift]
  w.float <- get_param(dat, "w.float")[shift]
  n.float.res <- get_param(dat, "n.float.res")
  env.w.res <- get_param(dat, "env.w.res")[shift]
  prec.mult <- get_param(dat, "prec.mult")
  hh.mult <- get_param(dat, "hh.mult")
  intro.staff.days <- get_param(dat, "intro.staff.days")
  n.units <- get_param(dat, "n.units")

  ob_active <- isTRUE(dat$run$ob_active)
  cohort <- ob_active && get_param(dat, "ctl.cohort") == 1
  handhyg <- ob_active && get_param(dat, "ctl.handhyg") == 1
  hands <- if (handhyg) hh.mult else 1

  n <- length(active)
  is_res <- active == 1 & type == "res"
  present <- is_res | (active == 1 & type == "staff" & on_duty == 1)

  infness <- numeric(n)
  infness[status == "i"] <- 1
  infness[status == "i" & vomit == 1] <- vomit.mult
  infness[status == "a" & type == "res"] <- rel.inf.a[1]
  infness[status == "a" & type == "staff"] <- rel.inf.a[2]
  infness[!present] <- 0
  sus <- present & status == "s"

  rows <- list()
  add_pairs <- function(rows, a, b, w, route) {
    # Both directions of each pair, keeping infectious-to-susceptible ones.
    if (length(a) == 0 || all(w == 0)) return(rows)
    w <- rep_len(w, length(a))
    inf <- c(a, b)
    rec <- c(b, a)
    ww <- c(w, w)
    keep <- infness[inf] > 0 & sus[rec] & ww > 0
    if (!any(keep)) return(rows)
    rows[[length(rows) + 1]] <- list(
      sus = rec[keep], inf = inf[keep], route = rep(route, sum(keep)),
      lambda = beta * ww[keep] * infness[inf[keep]])
    return(rows)
  }

  # Most shifts of a season have no one infectious and a clean environment;
  # the contact and pickup blocks below are skipped on those shifts.
  env <- dat$run$env
  if (is.null(env)) env <- numeric(n.units)

  if (any(infness > 0)) {
    # Layer 1: roommates.
    el <- get_edgelist(dat, network = 1)
    if (NROW(el) > 0) rows <- add_pairs(rows, el[, 1], el[, 2], w.room, "room")

    # Layer 2: residents of the same unit in the dining room and common areas;
    # isolated residents stay in their rooms.
    el <- get_edgelist(dat, network = 2)
    if (NROW(el) > 0 && w.common > 0) {
      ok <- isolated[el[, 1]] == 0 & isolated[el[, 2]] == 0
      rows <- add_pairs(rows, el[ok, 1], el[ok, 2], w.common, "common")
    }

    # Layer 3: staff on the same team (unit and shift, or the float team), when
    # both are at work.
    el <- get_edgelist(dat, network = 3)
    if (NROW(el) > 0) {
      ok <- on_duty[el[, 1]] == 1 & on_duty[el[, 2]] == 1
      rows <- add_pairs(rows, el[ok, 1], el[ok, 2], w.team, "team")
    }

    # Care contacts with the assigned CNA and nurse, built from attributes.
    res_ids <- which(is_res & !is.na(carer))
    cna_pos <- match(carer[res_ids], staff_id)
    mult <- ifelse(isolated[res_ids] == 1, prec.mult, 1) * hands
    rows <- add_pairs(rows, res_ids, cna_pos, w.care * mult, "care")
    res_ids <- which(is_res & !is.na(nurse_of))
    nurse_pos <- match(nurse_of[res_ids], staff_id)
    mult <- ifelse(isolated[res_ids] == 1, prec.mult, 1) * hands
    rows <- add_pairs(rows, res_ids, nurse_pos, w.nurse * mult, "nurse")

    # Float staff (activities, therapy, dietary, housekeeping) each see a few
    # residents who are not isolated, anywhere in the building, or on their
    # assigned unit while staff are cohorted during an outbreak.
    floats <- which(present & role == "float")
    if (length(floats) > 0 && w.float > 0) {
      fa <- fb <- integer(0)
      for (f in floats) {
        pool <- which(is_res & isolated == 0)
        if (cohort) pool <- pool[unit[pool] == float_unit[f]]
        if (length(pool) == 0) next
        pick <- pool[sample.int(length(pool), min(n.float.res, length(pool)))]
        fa <- c(fa, rep(f, length(pick)))
        fb <- c(fb, pick)
      }
      rows <- add_pairs(rows, fa, fb, w.float * hands, "float")
    }
  }

  # The environment: everyone present picks up from the unit they are on.
  # Residents are on their own unit; unit staff on theirs; float staff spread
  # their shift over all units, or stay on one while cohorted.
  if (any(env > 0)) {
    w_on <- matrix(0, n, n.units)
    for (u in seq_len(n.units)) {
      w_on[is_res & unit == u, u] <- env.w.res
      w_on[present & type == "staff" & role != "float" & unit == u, u] <- 1
      if (cohort) {
        w_on[present & role == "float" & float_unit == u, u] <- 1
      } else {
        w_on[present & role == "float", u] <- 1 / n.units
      }
    }
    n_on <- pmax(colSums(w_on > 0 & present), 1)
    pick <- ifelse(type == "staff", hands, 1)
    for (u in seq_len(n.units)) {
      if (env[u] <= 0) next
      j <- which(sus & w_on[, u] > 0)
      if (length(j) == 0) next
      rows[[length(rows) + 1]] <- list(
        sus = j, inf = rep(NA_integer_, length(j)), route = rep("env", length(j)),
        lambda = pick[j] * w_on[j, u] * env[u] / n_on[u])
    }
  }

  route_names <- c("room", "common", "team", "care", "nurse", "float", "env")
  n_route <- setNames(numeric(length(route_names)), route_names)
  n_env_res <- n_env_staff <- 0
  status_new <- status
  if (length(rows) > 0) {
    get_v <- function(v) unlist(lapply(rows, `[[`, v))
    lambda <- get_v("lambda")
    pos <- which(lambda > 0)
    t <- rexp(length(pos), rate = lambda[pos])
    hit <- pos[t < 1]
    ex <- data.frame(sus = get_v("sus")[hit], inf = get_v("inf")[hit],
                     route = get_v("route")[hit], t = t[t < 1])
    if (nrow(ex) > 0) {
      ex <- ex[order(ex$sus, ex$t), , drop = FALSE]
      ex <- ex[!duplicated(ex$sus), , drop = FALSE]
      status_new[ex$sus] <- "e"
      n_route[] <- table(factor(ex$route, levels = route_names))
      n_env_res <- sum(ex$route == "env" & type[ex$sus] == "res")
      n_env_staff <- sum(ex$route == "env" & type[ex$sus] == "staff")

      # Person-to-person transmissions go in the transmission matrix, which
      # needs an infector; environmental infections are counted by route.
      del <- ex[ex$route != "env", c("sus", "inf", "route")]
      if (nrow(del) > 0) {
        del$at <- at
        del$susType <- type[del$sus]
        del$infType <- type[del$inf]
        del$infStatus <- status[del$inf]
        dat <- set_transmat(dat, del, at)
      }
      route <- get_attr(dat, "infRoute")
      route[ex$sus] <- ex$route
      dat <- set_attr(dat, "infRoute", route)
      infTime <- get_attr(dat, "infTime")
      infTime[ex$sus] <- at
      dat <- set_attr(dat, "infTime", infTime)
    }
  }

  # Community introductions. The season's introductions are drawn before the
  # simulation (in model.R) and passed in as the days on which a staff member
  # is exposed outside the facility. On the night step of each such day, a
  # staff member is drawn at random and infected if susceptible.
  n_intro <- 0
  k <- sum(intro.staff.days == day)
  if (shift == 3 && k > 0) {
    st <- which(active == 1 & type == "staff")
    hit <- unique(st[sample.int(length(st), k, replace = TRUE)])
    hit <- hit[status_new[hit] == "s"]
    if (length(hit) > 0) {
      status_new[hit] <- "e"
      route <- get_attr(dat, "infRoute")
      route[hit] <- "community"
      dat <- set_attr(dat, "infRoute", route)
      infTime <- get_attr(dat, "infTime")
      infTime[hit] <- at
      dat <- set_attr(dat, "infTime", infTime)
      n_intro <- length(hit)
    }
  }

  # Draw the incubation period of everyone infected this step.
  new_e <- which(status_new == "e" & status != "e")
  if (length(new_e) > 0) {
    stage_end <- get_attr(dat, "stage_end")
    days <- rlnorm(length(new_e), get_param(dat, "inc.meanlog"),
                   get_param(dat, "inc.sdlog"))
    stage_end[new_e] <- at + pmax(1, round(days * 3))
    dat <- set_attr(dat, "stage_end", stage_end)
  }
  dat <- set_attr(dat, "status", status_new)

  for (r in route_names[route_names != "env"]) {
    dat <- set_epi(dat, paste0("inf.", r), n_route[[r]])
  }
  dat <- set_epi(dat, "inf.env.res", n_env_res)
  dat <- set_epi(dat, "inf.env.staff", n_env_staff)
  dat <- set_epi(dat, "intro.staff", n_intro)
  return(dat)
}


# The environmental reservoir ---------------------------------------------------

environ <- function(dat, at) {
  # One reservoir per unit, E_u, in arbitrary units of contamination. Each
  # shift, every infectious person on the unit adds env.dep times their
  # infectiousness times their time on the unit (for staff, times the hand
  # hygiene multiplier while soap and water is in use); the reservoir decays at
  # env.decay per day; routine cleaning at the start of each day shift removes
  # a fraction clean.routine; and during an outbreak with enhanced cleaning,
  # hypochlorite cleaning on the day and evening shifts removes clean.bleach.
  # Cleaning as discrete, scheduled pulses rather than a faster decay rate
  # keeps its timing, which matters for telling the routes apart.

  shift <- (at - 1) %% 3 + 1

  active <- get_attr(dat, "active")
  type <- get_attr(dat, "type")
  role <- get_attr(dat, "role")
  unit <- get_attr(dat, "unit")
  status <- get_attr(dat, "status")
  vomit <- get_attr(dat, "vomit")
  on_duty <- get_attr(dat, "on_duty")
  float_unit <- get_attr(dat, "float_unit")

  n.units <- get_param(dat, "n.units")
  env.dep <- get_param(dat, "env.dep")
  env.decay <- get_param(dat, "env.decay")
  env.w.res <- get_param(dat, "env.w.res")[shift]
  rel.inf.a <- get_param(dat, "rel.inf.a")
  vomit.mult <- get_param(dat, "vomit.mult")
  clean.routine <- get_param(dat, "clean.routine")
  clean.bleach <- get_param(dat, "clean.bleach")
  hh.mult <- get_param(dat, "hh.mult")

  ob_active <- isTRUE(dat$run$ob_active)
  bleach <- ob_active && get_param(dat, "ctl.clean") == 1
  cohort <- ob_active && get_param(dat, "ctl.cohort") == 1
  hands <- if (ob_active && get_param(dat, "ctl.handhyg") == 1) hh.mult else 1

  env <- dat$run$env
  if (is.null(env)) env <- numeric(n.units)

  is_res <- active == 1 & type == "res"
  present <- is_res | (active == 1 & type == "staff" & on_duty == 1)
  shed <- numeric(length(active))
  shed[status == "i"] <- 1
  shed[status == "i" & vomit == 1] <- vomit.mult
  shed[status == "a" & type == "res"] <- rel.inf.a[1]
  shed[status == "a" & type == "staff"] <- rel.inf.a[2]
  shed[type == "staff"] <- shed[type == "staff"] * hands
  shed[!present] <- 0

  for (u in seq_len(n.units)) {
    dep <- sum(shed[is_res & unit == u]) * env.w.res +
      sum(shed[present & type == "staff" & role != "float" & unit == u])
    if (cohort) {
      dep <- dep + sum(shed[present & role == "float" & float_unit == u])
    } else {
      dep <- dep + sum(shed[present & role == "float"]) / n.units
    }
    env[u] <- env[u] * exp(-env.decay / 3) + env.dep * dep
  }
  if (shift == 1) env <- env * (1 - clean.routine)
  if (bleach && shift %in% c(1, 2)) env <- env * (1 - clean.bleach)
  # Contamination that has decayed to a negligible level is set to zero, so
  # the infection module can skip the pickup step on quiet shifts.
  env[env < 1e-6] <- 0

  dat$run$env <- env
  for (u in seq_len(n.units)) dat <- set_epi(dat, paste0("env.u", u), env[u])
  return(dat)
}


# Natural history ----------------------------------------------------------------

progress <- function(dat, at) {
  # e -> i (symptomatic) or a (asymptomatic shedding) at the end of the
  # incubation period; i -> a (post-symptomatic shedding) when symptoms end;
  # a -> r when shedding ends. Durations are drawn on entry to each state and
  # stored as the step the state ends. At symptom onset each case draws
  # whether they vomit, and each staff case whether they keep working while
  # ill. Every onset is logged in dat$run$onsets, which keeps the record of
  # residents who are later discharged.

  active <- get_attr(dat, "active")
  type <- get_attr(dat, "type")
  role <- get_attr(dat, "role")
  unit <- get_attr(dat, "unit")
  status <- get_attr(dat, "status")
  stage_end <- get_attr(dat, "stage_end")
  vomit <- get_attr(dat, "vomit")
  presentee <- get_attr(dat, "presentee")
  onsetTime <- get_attr(dat, "onsetTime")
  symEnd <- get_attr(dat, "symEnd")

  p.asym <- get_param(dat, "p.asym")
  p.vomit <- get_param(dat, "p.vomit")
  p.present <- get_param(dat, "p.present")
  sym.mean <- get_param(dat, "sym.mean")
  sym.shape <- get_param(dat, "sym.shape")
  shed.mean <- get_param(dat, "shed.mean")
  shed.shape <- get_param(dat, "shed.shape")

  due <- active == 1 & !is.na(stage_end) & stage_end <= at

  # Shedding ends.
  ids <- which(due & status == "a")
  status[ids] <- "r"
  stage_end[ids] <- NA

  # Symptoms end; shedding continues at a low level.
  ids <- which(due & status == "i")
  if (length(ids) > 0) {
    status[ids] <- "a"
    symEnd[ids] <- at
    days <- rgamma(length(ids), shape = shed.shape, scale = shed.mean / shed.shape)
    stage_end[ids] <- at + pmax(1, round(days * 3))
  }

  # Incubation ends: symptoms or asymptomatic shedding.
  ids <- which(due & status == "e")
  n_onset_res <- n_onset_staff <- 0
  on_res <- integer(0)
  if (length(ids) > 0) {
    asym <- runif(length(ids)) < p.asym
    a_ids <- ids[asym]
    s_ids <- ids[!asym]
    status[a_ids] <- "a"
    days <- rgamma(length(a_ids), shape = shed.shape, scale = shed.mean / shed.shape)
    stage_end[a_ids] <- at + pmax(1, round(days * 3))
    if (length(s_ids) > 0) {
      status[s_ids] <- "i"
      onsetTime[s_ids] <- at
      vomit[s_ids] <- as.numeric(runif(length(s_ids)) < p.vomit)
      presentee[s_ids] <- ifelse(type[s_ids] == "staff",
                                 as.numeric(runif(length(s_ids)) < p.present), 0)
      days <- rgamma(length(s_ids), shape = sym.shape, scale = sym.mean / sym.shape)
      stage_end[s_ids] <- at + pmax(1, round(days * 3))
      log <- data.frame(at = at, type = type[s_ids], role = role[s_ids],
                        unit = unit[s_ids])
      dat$run$onsets <- rbind(dat$run$onsets, log)
      n_onset_res <- sum(type[s_ids] == "res")
      n_onset_staff <- sum(type[s_ids] == "staff")
      on_res <- s_ids[type[s_ids] == "res"]
    }
  }
  for (u in seq_len(get_param(dat, "n.units"))) {
    dat <- set_epi(dat, paste0("onset.res.u", u), sum(unit[on_res] == u))
  }

  dat <- set_attr(dat, "status", status)
  dat <- set_attr(dat, "stage_end", stage_end)
  dat <- set_attr(dat, "vomit", vomit)
  dat <- set_attr(dat, "presentee", presentee)
  dat <- set_attr(dat, "onsetTime", onsetTime)
  dat <- set_attr(dat, "symEnd", symEnd)
  dat <- set_epi(dat, "onset.staff", n_onset_staff)
  return(dat)
}


# Outbreak detection and control ------------------------------------------------

outbreak <- function(dat, at) {
  # The facility declares an outbreak when trig.n or more people (residents
  # or staff) have had symptom onset within the last trig.window steps, after
  # a reporting delay of trig.delay steps, and declares it over after
  # end.window steps with no new onset. If cases resume and the outbreak is
  # declared again within ob.gap steps of being declared over, it is the same
  # outbreak resurging rather than a new one. While an outbreak is declared,
  # the measures switched on in the parameters apply:
  #
  #   ctl.isolate  symptomatic residents stay in their rooms (no dining or
  #                common-area contact, contact precautions for care) until
  #                iso.hold steps after symptoms end
  #   ctl.exclude  symptomatic staff stay home until excl.hold steps after
  #                symptoms end
  #   ctl.cohort   float staff work on one unit each (read by infect, environ)
  #   ctl.close    no admissions to units with a case in this outbreak
  #   ctl.clean    hypochlorite cleaning twice a day (read by environ)
  #   ctl.handhyg  soap and water after care (read by infect)
  #
  # The current outbreak is kept in dat$run$cur_ob; each finished outbreak is
  # written as a row of dat$run$ob_log.

  nsteps <- get_control(dat, "nsteps")
  trig.n <- get_param(dat, "trig.n")
  trig.window <- get_param(dat, "trig.window")
  trig.delay <- get_param(dat, "trig.delay")
  end.window <- get_param(dat, "end.window")
  ob.gap <- get_param(dat, "ob.gap")

  if (is.null(dat$run$ob_log)) dat$run$ob_log <- data.frame()
  onsets <- dat$run$onsets
  on_at <- if (is.null(onsets)) numeric(0) else onsets$at
  ob <- dat$run$cur_ob

  # Declare, or declare again.
  if (is.null(ob) || !ob$active) {
    if (is.null(dat$run$pending) && sum(on_at > at - trig.window) >= trig.n) {
      dat$run$pending <- at + trig.delay
    }
    if (!is.null(dat$run$pending) && at >= dat$run$pending) {
      dat$run$pending <- NULL
      if (is.null(ob)) {
        ob <- list(first = min(on_at[on_at > at - trig.window - trig.delay]),
                   declared = at, declarations = 0)
      }
      ob$active <- TRUE
      ob$declarations <- ob$declarations + 1
    }
  } else if (sum(on_at > at - end.window) == 0) {
    # Declared over.
    ob$active <- FALSE
    ob$over <- at
  }

  # Write the outbreak to the log once it has stayed over for ob.gap steps,
  # or at the end of the season.
  if (!is.null(ob) && (at == nsteps || (!ob$active && at >= ob$over + ob.gap))) {
    o <- onsets[on_at >= ob$first, , drop = FALSE]
    if (!ob$active) o <- o[o$at <= ob$over, , drop = FALSE]
    rec <- data.frame(
      first = ob$first, declared = ob$declared, last = max(o$at),
      declarations = ob$declarations, ongoing = ob$active,
      cases.res = sum(o$type == "res"), cases.staff = sum(o$type == "staff"),
      units = length(unique(o$unit[o$role != "float"])))
    dat$run$ob_log <- rbind(dat$run$ob_log, rec)
    ob <- NULL
  }
  dat$run$cur_ob <- ob
  ob_active <- !is.null(ob) && ob$active
  dat$run$ob_active <- ob_active

  # Units touched by the current outbreak (for closure to admissions).
  if (ob_active) {
    o <- onsets[on_at >= ob$first, , drop = FALSE]
    dat$run$closed <- unique(o$unit[o$role != "float"])
  } else {
    dat$run$closed <- NULL
  }

  # Isolation and exclusion.
  active <- get_attr(dat, "active")
  type <- get_attr(dat, "type")
  status <- get_attr(dat, "status")
  isolated <- get_attr(dat, "isolated")
  excluded <- get_attr(dat, "excluded")
  hold_until <- get_attr(dat, "hold_until")
  symEnd <- get_attr(dat, "symEnd")

  if (ob_active && get_param(dat, "ctl.isolate") == 1) {
    isolated[active == 1 & type == "res" & status == "i"] <- 1
  }
  if (ob_active && get_param(dat, "ctl.exclude") == 1) {
    excluded[active == 1 & type == "staff" & status == "i"] <- 1
  }
  held <- (isolated == 1 | excluded == 1)
  # Start the hold clock when symptoms end, and release when it runs out.
  start <- held & status != "i" & is.na(hold_until)
  hold <- ifelse(type[start] == "res", get_param(dat, "iso.hold"),
                 get_param(dat, "excl.hold"))
  hold_until[start] <- symEnd[start] + hold
  release <- held & !is.na(hold_until) & at >= hold_until
  isolated[release] <- 0
  excluded[release] <- 0
  hold_until[release] <- NA

  dat <- set_attr(dat, "isolated", isolated)
  dat <- set_attr(dat, "excluded", excluded)
  dat <- set_attr(dat, "hold_until", hold_until)
  dat <- set_epi(dat, "ob.active", as.numeric(ob_active))
  dat <- set_epi(dat, "n.isolated", sum(isolated == 1 & active == 1))
  dat <- set_epi(dat, "n.excluded", sum(excluded == 1 & active == 1))
  return(dat)
}


# Discharges and admissions --------------------------------------------------------

discharge <- function(dat, at) {
  # On each day shift, residents leave (discharge home, transfer, or death)
  # with a daily probability set by the mean length of stay of their unit.
  # Their beds go on the free list, keyed by unit and bed, never by position.
  shift <- (at - 1) %% 3 + 1
  if (shift != 1) return(dat)

  active <- get_attr(dat, "active")
  type <- get_attr(dat, "type")
  unit <- get_attr(dat, "unit")
  bed <- get_attr(dat, "bed")
  room <- get_attr(dat, "room")
  los <- get_param(dat, "los")

  res <- which(active == 1 & type == "res")
  out <- res[runif(length(res)) < 1 - exp(-1 / los[unit[res]])]
  if (length(out) > 0) {
    active[out] <- 0
    exitTime <- get_attr(dat, "exitTime")
    exitTime[out] <- at
    dat <- set_attr(dat, "active", active)
    dat <- set_attr(dat, "exitTime", exitTime)
    freed <- data.frame(unit = unit[out], bed = bed[out], room = room[out])
    dat$run$free <- rbind(dat$run$free, freed)
  }
  dat <- set_epi(dat, "discharges", length(out))
  return(dat)
}

admit <- function(dat, at) {
  # On each day shift, each free bed is filled with probability fill.prob,
  # except on units closed to admissions during an outbreak. Introductions by
  # admission are drawn before the simulation (in model.R) as days on which
  # an infected person is admitted; the next new resident on or after each
  # such day is that admission. An infected admission was infected at a
  # random time before arrival, so they arrive incubating, symptomatic, or
  # shedding without symptoms in proportion to the mean time spent in each
  # state, partway through it. A resident who arrives symptomatic is placed
  # on contact precautions in their room on arrival, as is standard practice,
  # until 48 hours after symptoms end. Other new residents are immune with
  # probability p.immune and otherwise susceptible. Every attribute is
  # appended, and the room and unit clique layers place the new node with
  # the room and unit read back from its attributes.
  shift <- (at - 1) %% 3 + 1
  day <- (at - 1) %/% 3 + 1
  if (shift != 1) return(dat)
  free <- dat$run$free
  n_new <- 0
  owed <- dat$run$adm_owed
  if (is.null(owed)) owed <- 0
  dat$run$adm_owed <- owed + sum(get_param(dat, "intro.adm.days") == day)
  if (!is.null(free) && nrow(free) > 0) {
    closed <- if (get_param(dat, "ctl.close") == 1) dat$run$closed else NULL
    open <- !(free$unit %in% closed)
    fill <- open & runif(nrow(free)) < get_param(dat, "fill.prob")
    new <- free[fill, , drop = FALSE]
    dat$run$free <- free[!fill, , drop = FALSE]
    n_new <- nrow(new)
  }
  if (n_new == 0) {
    dat <- set_epi(dat, "intro.adm", 0)
    dat <- set_epi(dat, "admissions", 0)
    return(dat)
  }

  inc.meanlog <- get_param(dat, "inc.meanlog")
  inc.sdlog <- get_param(dat, "inc.sdlog")
  sym.mean <- get_param(dat, "sym.mean")
  shed.mean <- get_param(dat, "shed.mean")
  p.asym <- get_param(dat, "p.asym")

  status <- ifelse(runif(n_new) < get_param(dat, "p.immune"), "r", "s")
  n_inf <- min(dat$run$adm_owed, n_new)
  dat$run$adm_owed <- dat$run$adm_owed - n_inf
  stage_end <- rep(NA_real_, n_new)
  isolated <- vomit <- rep(0, n_new)
  if (n_inf > 0) {
    ids <- seq_len(n_inf)
    # Mean days in each state over an infection: incubation, symptoms (for
    # those who have them), and shedding without symptoms.
    t_state <- c(e = exp(inc.meanlog + inc.sdlog^2 / 2),
                 i = (1 - p.asym) * sym.mean, a = shed.mean)
    status[ids] <- sample(names(t_state), n_inf, replace = TRUE, prob = t_state)
    full <- ifelse(status[ids] == "e", rlnorm(n_inf, inc.meanlog, inc.sdlog),
                   ifelse(status[ids] == "i", sym.mean, shed.mean))
    stage_end[ids] <- at + pmax(1, round(runif(n_inf) * full * 3))
    isolated[ids] <- as.numeric(status[ids] == "i")
    vomit[ids] <- as.numeric(status[ids] == "i" &
                               runif(n_inf) < get_param(dat, "p.vomit"))
  }
  infected <- status %in% c("e", "i", "a")

  dat <- append_core_attr(dat, at, n_new)
  dat <- append_attr(dat, "status", status, n_new)
  dat <- append_attr(dat, "infTime", ifelse(infected, at, NA), n_new)
  dat <- append_attr(dat, "infRoute", ifelse(infected, "admission", NA), n_new)
  dat <- append_attr(dat, "stage_end", stage_end, n_new)
  dat <- append_attr(dat, "type", "res", n_new)
  dat <- append_attr(dat, "role", "resident", n_new)
  dat <- append_attr(dat, "unit", new$unit, n_new)
  dat <- append_attr(dat, "res_unit", new$unit, n_new)
  dat <- append_attr(dat, "room", new$room, n_new)
  dat <- append_attr(dat, "bed", new$bed, n_new)
  dat <- append_attr(dat, "team", NA, n_new)
  for (a in c("staff_id", "staff_shift", "pool_k", "off1", "off2",
              "float_unit", "carer", "nurse_of", "onsetTime", "symEnd",
              "hold_until")) {
    dat <- append_attr(dat, a, NA, n_new)
  }
  for (a in c("on_duty", "presentee", "excluded")) {
    dat <- append_attr(dat, a, 0, n_new)
  }
  dat <- append_attr(dat, "isolated", isolated, n_new)
  dat <- append_attr(dat, "vomit", vomit, n_new)
  dat <- set_epi(dat, "intro.adm", n_inf)
  dat <- set_epi(dat, "admissions", n_new)
  return(dat)
}
