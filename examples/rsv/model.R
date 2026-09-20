##
## RSV: Age-Stratified SEIR Over a Household + Community Multilayer Network
## EpiModel Gallery (https://github.com/EpiModel/EpiModel-Gallery)
##
## Authors: Samuel M. Jenness (Emory University)
## Date: September 2026
##

# Load EpiModel
suppressMessages(library(EpiModel))

# Standard Gallery unit test lines
rm(list = ls())
eval(parse(text = print(commandArgs(TRUE)[1])))

if (interactive()) {
  N <- 10000
  nsims <- 10
  ncores <- 5
  nsteps <- 150
} else {
  N <- 1000
  nsims <- 1
  ncores <- 1
  nsteps <- 50
}


# 1. Population and Households ----------------------------------------------

# The population is built from households rather than from individuals. Each
# entry below is a household type: the names list its members by age group
# (infant < 1 year, young 1-4, school 5-17, adult 18-64, elderly 65+) and the
# values are the probability that a sampled household is of that type. The
# mix was chosen so that the person-level age shares approximate the United
# States (about 1.1% infants, 5% young, 17% school-age, 58% adults, 19%
# older adults, mean household size 2.3), every infant lives with at least
# one adult, about 40% of infants have an older sibling, and about 28% of
# older adults live alone.
hh_types <- c(
  "adult"                          = 0.150,
  "adult adult"                    = 0.184,
  "adult adult adult"              = 0.050,
  "elderly"                        = 0.120,
  "elderly elderly"                = 0.120,
  "adult elderly"                  = 0.030,
  "adult adult infant"             = 0.010,
  "adult adult infant young"       = 0.005,
  "adult adult infant school"      = 0.006,
  "adult adult infant elderly"     = 0.003,
  "adult infant"                   = 0.002,
  "adult adult young"              = 0.030,
  "adult adult young young"        = 0.010,
  "adult adult young school"       = 0.035,
  "adult young"                    = 0.010,
  "adult adult school"             = 0.060,
  "adult adult school school"      = 0.070,
  "adult adult school school school" = 0.020,
  "adult school"                   = 0.035,
  "adult school school"            = 0.020,
  "adult adult school elderly"     = 0.020,
  "adult adult young elderly"      = 0.010
)
stopifnot(abs(sum(hh_types) - 1) < 1e-8)

# Samples households until the population reaches N, then truncates to N
# (at most the last household is cut short). Returns each person's age group
# and household id; households occupy consecutive node ids.
generate_households <- function(N, hh_types, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  members <- strsplit(names(hh_types), " ")
  sizes <- lengths(members)
  n_draw <- ceiling(1.5 * N / sum(hh_types * sizes)) + 10
  draw <- sample.int(length(hh_types), n_draw, replace = TRUE, prob = hh_types)
  draw <- draw[seq_len(which(cumsum(sizes[draw]) >= N)[1])]
  age <- unlist(members[draw])[1:N]
  hh_id <- rep(seq_along(draw), sizes[draw])[1:N]
  list(age = age, hh_id = hh_id)
}

# Every pair of co-residents is a household edge, so each household is a
# clique. This edgelist is the whole household layer: it is fixed for the
# season and is handed to the infection module as a parameter.
household_edgelist <- function(hh_id) {
  members <- split(seq_along(hh_id), hh_id)
  members <- members[lengths(members) > 1]
  el <- do.call(rbind, lapply(members, function(m) t(combn(m, 2))))
  unname(el)
}

pop <- generate_households(N, hh_types, seed = 123)
age <- pop$age
hh_id <- pop$hh_id
hh_el <- household_edgelist(hh_id)

counts <- table(factor(age, levels = c("adult", "elderly", "infant",
                                       "school", "young")))
cat(sprintf("N = %d in %d households. Age counts:\n", N, max(hh_id)))
print(counts)
cat("\nHousehold size distribution:\n")
print(table(tabulate(hh_id)))

# Household-layer diagnostics: mean degree by age is a property of the
# clique structure, read straight from the edgelist.
hh_deg <- tabulate(c(hh_el[, 1], hh_el[, 2]), nbins = N)
deg_hh_by_age <- round(tapply(hh_deg, age, mean), 2)
cat("\nHousehold-layer mean degree by age:\n"); print(deg_hh_by_age)


# 2. Community Network (TERGM layer) ----------------------------------------

# The community layer is an ERGM with an edges term plus nodemix("age"),
# which gives every cell of the age-by-age mixing matrix its own target.
# Targeting the full matrix (rather than a subset of cells) fixes the degree
# of every age group by design; cells left out of nodemix would otherwise
# absorb whatever edge count remains from the edges target, spread uniformly
# over the untargeted dyads.
#
# Both attributes are set on the network so that netsim carries them into
# the simulation: age drives contacts, susceptibility, severity, and
# eligibility; hh_id makes household-targeted strategies definable.

nw <- network_initialize(N)
nw <- set_vertex_attribute(nw, "age", age)
nw <- set_vertex_attribute(nw, "hh_id", hh_id)

# Canonical nodemix cell names in ergm's order: the upper triangle of the
# mixing matrix, column-major, with alphabetical levels.
cells <- sub("^mix\\.age\\.", "",
             names(summary(nw ~ nodemix("age", levels2 = TRUE))))

# Converts a contact profile into an edge-count target for every cell. Each
# profile entry "a.b" is the mean number of b-partners per a-node; the edge
# count is that value times the size of group a, halved for within-group
# cells because each edge is counted from both ends. Targets are rounded to
# whole edges (floor of one) so ergm can match them exactly.
mix_targets <- function(profile, counts, cells) {
  out <- setNames(numeric(length(cells)), cells)
  for (nm in names(profile)) {
    ab <- strsplit(nm, ".", fixed = TRUE)[[1]]
    n_a <- as.numeric(counts[ab[1]])
    e <- if (ab[1] == ab[2]) n_a * profile[[nm]] / 2 else n_a * profile[[nm]]
    cell <- paste(sort(ab), collapse = ".")
    stopifnot(cell %in% cells)
    out[cell] <- out[cell] + e
  }
  pmax(round(out), 1)
}

# Daily contacts outside the household, with the strong age assortativity
# seen in contact surveys. School-age children have the highest degree and
# mix mostly with each other (the school amplifier); infants have the
# fewest contacts. Cross-age entries are written from the perspective of
# the smaller group, e.g. school.adult = 1.5 means each school-age child
# has 1.5 adult community contacts on average.
com_profile <- c(
  school.school   = 5.0,   # classmates
  adult.adult     = 4.0,   # workplace, social
  young.young     = 2.5,   # daycare
  elderly.elderly = 1.5,
  school.adult    = 1.5,   # teachers, coaches, friends' parents
  young.adult     = 2.0,   # daycare staff
  elderly.adult   = 1.2,
  school.elderly  = 0.3,   # grandparents
  young.school    = 0.5,
  infant.adult    = 0.8,   # non-household caregivers
  infant.school   = 0.2,
  infant.young    = 0.2,
  infant.elderly  = 0.2,
  young.elderly   = 0.3,
  infant.infant   = 0.02
)
com_targets <- mix_targets(com_profile, counts, cells)

# levels2 = -1 drops the first cell (adult.adult) as the reference; its
# count is implied by the edges target minus the other 14 cells.
formation <- ~edges + nodemix("age", levels2 = -1)

# With 14 targeted cells, ergm's default simulated annealing (SAN) step can
# fail to hit every target exactly, in which case ergm falls back to slow
# MCMC estimation. A larger SAN budget lets this dyad-independent model be
# fit by maximum pseudolikelihood in a few seconds.
san_ctrl <- control.ergm(SAN = control.san(SAN.maxit = 20,
                                           SAN.nsteps = 2^21))

# Community ties last one day: each day brings a new set of contacts.
coef.diss_com <- dissolution_coefs(~offset(edges), duration = 1)

cat("\nFitting community layer ERGM...\n")
est_com <- netest(nw, formation,
                  target.stats = c(sum(com_targets),
                                   as.numeric(com_targets[-1])),
                  coef.diss = coef.diss_com,
                  set.control.ergm = san_ctrl, verbose = FALSE)

# Diagnostics: realized mean degree by age from static simulations of the
# formation model. This is the check that the mixing targets translated
# into the intended age-specific contact rates.
mean_degree_by_age <- function(est, counts) {
  dx <- netdx(est, nsims = 10, dynamic = FALSE, verbose = FALSE,
              nwstats.formula = ~edges + nodemix("age", levels2 = TRUE) +
                                 nodefactor("age", levels = TRUE))
  st <- colMeans(do.call(rbind, lapply(dx$stats, as.matrix)))
  st <- st[grep("^nodefactor", names(st))]
  names(st) <- sub("nodefactor.age.", "", names(st))
  round(st / as.numeric(counts[names(st)]), 2)
}
deg_by_age <- rbind(household = deg_hh_by_age[names(counts)],
                    community = mean_degree_by_age(est_com, counts))
cat("\nRealized mean degree by age:\n"); print(deg_by_age)


# 3. Disease + Intervention Parameters --------------------------------------

source("examples/rsv/module-fx.R")

# Natural history (daily time step):
#   ei.rate        -- 1/4  mean 4-day latent period
#   ip.rate        -- 1/2  mean 2-day presymptomatic infectious period
#   ir.rate        -- 1/7  mean 7-day symptomatic or asymptomatic period
#   asymp.prob     -- 0.3  share of infections that never become symptomatic
#   asymp.inf.mult -- 0.5  asymptomatic cases are half as infectious
#
# Prior immunity. Nearly everyone is infected with RSV by age 2 and
# reinfected throughout life, with each reinfection milder and less likely
# (Glezen et al. 1986; Pitzer et al. 2015). Age is used here as a proxy for
# exposure history: sus.mult scales the per-contact probability of
# infection for a susceptible node of each age relative to a never-infected
# infant. The per-contact probabilities (inf.prob.household,
# inf.prob.community) and sus.mult were chosen together so that the
# baseline season produces seasonal attack rates near the published age
# gradient (about 60% of infants, 40-60% of 1-4 year olds, roughly 20-30%
# of school-age children, 7% of adults, and 3-7% of older adults; Glezen
# et al. 1986, Hall et al. 2001, Falsey et al. 2005). They are
# illustrative, not fitted.
#
# Immunization products. Each has two leaky components, following how the
# products are evaluated and how other RSV scenario models represent them:
#   *.eff.inf  -- reduction in the per-contact probability of infection
#                 (this is the only component that produces indirect
#                 protection of others through reduced transmission)
#   *.eff.hosp -- reduction in the hospitalization risk given infection
# The two components combine to 1 - (1 - eff.inf) * (1 - eff.hosp) = 0.80
# against hospitalization for both products, matching the first-season
# effectiveness assumed by the RSV Scenario Modeling Hub. Trial estimates
# against medically attended RSV illness are 74-80% for nirsevimab in
# infants (Hammitt et al. 2022; Simoes et al. 2023) and 67-83% for the
# older-adult vaccines (Papi et al. 2023; Walsh et al. 2023). The split
# between the two components is an assumption; setting eff.inf = 0 gives a
# product that protects only the recipient.
#
# Household targeting ("cocooning"). Because households are explicit, the
# adult co-residents of infants are a definable target group. The cocoon
# scenario gives them a HYPOTHETICAL product with the older-adult vaccine's
# infection-blocking component and no severity component, to ask how much
# infant protection blocking household transmission can deliver compared
# with the direct infant product. No such adult product is currently
# recommended; maternal vaccination protects the infant through antibody
# transfer, not by blocking the parent's transmission.
#
# Eligibility is simplified to all infants and all adults 65+. Current CDC
# guidance is nirsevimab (or clesrovimab, or maternal vaccination) for
# infants entering their first RSV season, and a single vaccine dose for
# adults 75+ and adults 50-74 at increased risk of severe disease.
#
# NPI window: between npi.start and npi.end the community-layer per-contact
# probability is multiplied by (1 - npi.mask.efficacy) and each community
# edge is kept with probability npi.contact.mult.

init <- init.net(i.num = round(0.01 * N))

# Base parameter set. Intervention coverages default to zero and the NPI
# window defaults to inactive, so the "none" scenario is this set as is.
# hh.pairs carries the fixed household edgelist into the infection module.
param_base <- param.net(
  inf.prob.household = 0.45,
  inf.prob.community = 0.10,
  sus.mult = c(infant = 1.00, young = 0.60, school = 0.16,
               adult = 0.08, elderly = 0.13),
  asymp.inf.mult = 0.5,
  ei.rate = 1 / 4,
  ip.rate = 1 / 2,
  ir.rate = 1 / 7,
  asymp.prob = 0.3,
  elderly.vax.coverage = 0,
  elderly.vax.eff.inf = 0.5,
  elderly.vax.eff.hosp = 0.6,
  infant.proph.coverage = 0,
  infant.proph.eff.inf = 0.3,
  infant.proph.eff.hosp = 0.71,
  cocoon.coverage = 0,
  cocoon.eff.inf = 0.5,
  npi.start = -1,
  npi.end = -1,
  npi.mask.efficacy = 0.4,
  npi.contact.mult = 0.7,
  hh.pairs = hh_el
)

# module.order is set explicitly. EpiModel's default runs user-supplied
# modules before the built-in ones, which would put progress() ahead of
# infect() within each step; the order here matches the natural history
# (transmission, then progression) and the safeguards in progress().
control <- control.net(
  type = NULL,
  nsims = nsims,
  ncores = ncores,
  nsteps = nsteps,
  tergmLite = TRUE,
  resimulate.network = TRUE,
  initAttr.FUN = init_attrs,
  infection.FUN = infect,
  progress.FUN = progress,
  module.order = c("resim_nets.FUN", "summary_nets.FUN", "initAttr.FUN",
                   "infection.FUN", "progress.FUN", "nwupdate.FUN",
                   "prevalence.FUN"),
  verbose = FALSE
)


# 4. Scenarios --------------------------------------------------------------

# One row per scenario, columns matching parameter names in param_base.
# Product coverage follows the "typical" 2026-27 assumptions of the RSV
# Scenario Modeling Hub (about 55-60% of infants, 50% of adults 75+).
scenarios.df <- data.frame(
  .scenario.id          = c("none", "elderly_vax", "infant_proph",
                            "both", "cocoon", "npi"),
  .at                   = 0,
  elderly.vax.coverage  = c(0, 0.5, 0, 0.5, 0, 0),
  infant.proph.coverage = c(0, 0, 0.6, 0.6, 0, 0),
  cocoon.coverage       = c(0, 0, 0, 0, 0.6, 0),
  npi.start             = c(-1, -1, -1, -1, -1, 30),
  npi.end               = c(-1, -1, -1, -1, -1, 90)
)
scenarios.list <- create_scenario_list(scenarios.df)

labels <- c(none = "No intervention",
            elderly_vax = "Older-adult vaccine (50%)",
            infant_proph = "Infant antibody (60%)",
            both = "Both products",
            cocoon = "Cocooning (60% of infant-household adults)",
            npi = "NPI (days 30-90)")

sims <- list()
for (scn in scenarios.list) {
  cat(sprintf("Scenario: %s\n", scn$id))
  sims[[scn$id]] <- netsim(est_com, use_scenario(param_base, scn),
                           init, control)
}


# 5. Analysis ---------------------------------------------------------------

# Hospitalization risk per infection by age. Hospitalization is not a
# compartment: expected hospitalizations are computed after the fact as
# infections times these risks. The values are illustrative, chosen so that
# the baseline season lands near the age pattern of RSV-NET hospitalization
# rates (highest in infants, then older adults, lowest in school-age
# children) given the attack rates above. For immunized people the risk is
# further multiplied by (1 - eff.hosp).
hosp_rate <- c(infant = 0.030, young = 0.007, school = 0.001,
               adult = 0.004, elderly = 0.030)
eff_hosp <- c(infant = param_base$infant.proph.eff.hosp,
              elderly = param_base$elderly.vax.eff.hosp)

age_groups <- c("infant", "young", "school", "adult", "elderly")
age_pop <- as.numeric(counts[age_groups])
names(age_pop) <- age_groups

# Season-end summary for one scenario, averaged over simulations
summarize_scenario <- function(sim, scn_row) {
  df <- as.data.frame(sim)
  last <- df[df$time == max(df$time), ]           # one row per simulation
  # Expected hospitalizations per simulation (rows) and age group (columns):
  # infections times the per-infection risk, minus the severity reduction
  # for infections among immunized people.
  hosp_sim <- sapply(age_groups, function(a) {
    last[[paste0("cuminf.", a)]] * hosp_rate[a]
  })
  hosp_sim <- matrix(hosp_sim, ncol = length(age_groups),
                     dimnames = list(NULL, age_groups))
  for (a in c("infant", "elderly")) {
    hosp_sim[, a] <- hosp_sim[, a] -
      last[[paste0("cuminf.", a, ".prot")]] * hosp_rate[a] * eff_hosp[a]
  }
  inf <- sapply(age_groups, function(a) mean(last[[paste0("cuminf.", a)]]))
  inf_prot <- c(infant = mean(last$cuminf.infant.prot),
                elderly = mean(last$cuminf.elderly.prot))
  n_prot <- c(infant = mean(last$n.infant.prot),
              elderly = mean(last$n.elderly.prot))
  attack_unprot <- (inf[c("infant", "elderly")] - inf_prot) /
                   (age_pop[c("infant", "elderly")] - n_prot)
  # Doses: product coverage times the eligible population, plus the adults
  # actually reached by household targeting (a count from the simulation).
  doses <- scn_row$elderly.vax.coverage * age_pop["elderly"] +
           scn_row$infant.proph.coverage * age_pop["infant"] +
           mean(last$n.cocoon)
  hosp <- colMeans(hosp_sim)
  # Share of infant infections acquired from household contacts, pooled
  # over the season and the simulations
  hh_share_infant <- sum(df$se.flow.infant.hh, na.rm = TRUE) /
    sum(df$se.flow.infant.hh + df$se.flow.infant.com, na.rm = TRUE)
  list(inf = inf, attack = inf / age_pop, hosp = hosp,
       hosp_per100k = 1e5 * hosp / age_pop,
       hosp_total_sims = 1e5 * rowSums(hosp_sim) / N,
       attack_unprot = attack_unprot, doses = as.numeric(doses),
       hh_share_infant = hh_share_infant)
}

res <- lapply(names(sims), function(s) {
  summarize_scenario(sims[[s]], scenarios.df[scenarios.df$.scenario.id == s, ])
})
names(res) <- names(sims)

cat("\n=== Cumulative attack rate by age (%) ===\n")
attack_tbl <- sapply(res, function(r) round(100 * r$attack, 1))
print(attack_tbl)

cat("\n=== Expected hospitalizations per 100,000 by age (total row: all ages) ===\n")
hosp100k_tbl <- rbind(sapply(res, function(r) round(r$hosp_per100k)),
                      total = sapply(res, function(r) round(1e5 * sum(r$hosp) / N)))
print(hosp100k_tbl)

cat("\n=== Total hospitalizations per 100,000: mean and range across simulations ===\n")
print(sapply(res, function(r) round(c(mean = mean(r$hosp_total_sims),
                                      min = min(r$hosp_total_sims),
                                      max = max(r$hosp_total_sims)))))

# Age groups each strategy is designed to protect. Hospitalizations averted
# and NNI are computed within these groups. Neither product has a measurable
# indirect effect (see the unimmunized attack rates below), and the
# between-simulation noise in the large untargeted strata is larger than a
# product's whole effect, so an all-ages difference would mostly be noise.
target_groups <- list(none = age_groups, elderly_vax = "elderly",
                      infant_proph = "infant", both = c("infant", "elderly"),
                      cocoon = "infant", npi = age_groups)
averted_in_target <- function(s) {
  tg <- target_groups[[s]]
  sum(res$none$hosp[tg]) - sum(res[[s]]$hosp[tg])
}
cat("\n=== Intervention summary (hospitalizations averted in the target groups) ===\n")
int_tbl <- data.frame(
  scenario = labels[names(res)],
  target = sapply(names(res), function(s)
    if (length(target_groups[[s]]) == 5) "all ages" else
      paste(target_groups[[s]], collapse = " + ")),
  doses = sapply(res, function(r) round(r$doses)),
  hosp_averted_per100k = sapply(names(res), function(s)
    round(1e5 * averted_in_target(s) / N, 1)),
  pct_averted = sapply(names(res), function(s)
    round(100 * averted_in_target(s) / sum(res$none$hosp[target_groups[[s]]]), 1)),
  NNI = sapply(names(res), function(s)
    if (res[[s]]$doses > 0) round(res[[s]]$doses / averted_in_target(s)) else NA),
  row.names = NULL
)
print(int_tbl)

# Where infants get infected: the share of infant infections acquired from
# household contacts. This is what limits cocooning, which blocks only the
# adult-to-infant part of that share.
cat("\n=== Share of infant infections acquired from household contacts (%) ===\n")
print(sapply(res, function(r) round(100 * r$hh_share_infant, 1)))

# Indirect protection: attack rate among UNimmunized infants and older
# adults relative to the no-intervention baseline. Any reduction here is
# transmission blocked by the eff.inf component in other people. Under
# cocooning no infant is immunized, so the infant row is the whole effect.
cat("\n=== Attack rate among unimmunized (%), by scenario ===\n")
print(sapply(res, function(r) round(100 * r$attack_unprot, 1)))

# Between-simulation variability of the baseline: CV of season-end
# cumulative infections by age. The infant stratum is small (about 1% of
# N), so it is the noisiest; this is the reason for the large default N.
if (nsims > 1) {
  df_none <- as.data.frame(sims$none)
  last_none <- df_none[df_none$time == max(df_none$time), ]
  cv <- sapply(age_groups, function(a) {
    v <- last_none[[paste0("cuminf.", a)]]; sd(v) / mean(v)
  })
  cat("\n=== Baseline CV of cumulative infections across simulations ===\n")
  print(round(cv, 2))
}


# 6. Plots ------------------------------------------------------------------

cols_scn <- c(none = "gray40", elderly_vax = "seagreen",
              infant_proph = "purple", both = "darkblue",
              cocoon = "darkorange", npi = "firebrick")
cols_age <- c(infant = "#3498db", young = "#f39c12", school = "#e74c3c",
              adult = "#27ae60", elderly = "#8e44ad")

## --- Plot 1: Age-stratified cumulative attack rates ---
par(mfrow = c(2, 3), mar = c(3, 3.5, 2, 1), mgp = c(2.2, 0.7, 0))
for (a in age_groups) {
  col_name <- paste0("cuminf.", a)
  curves <- lapply(names(sims), function(s) {
    df <- as.data.frame(sims[[s]])
    tapply(df[[col_name]], df$time, mean, na.rm = TRUE) / age_pop[a]
  })
  names(curves) <- names(sims)
  max_y <- max(sapply(curves, max, na.rm = TRUE))
  plot(NA, xlim = c(1, nsteps), ylim = c(0, max(max_y, 0.01) * 1.05),
       xlab = "Day", ylab = "Cumulative attack rate",
       main = paste0(toupper(substr(a, 1, 1)), substr(a, 2, nchar(a)),
                    " (N=", age_pop[a], ")"))
  for (s in names(sims)) {
    lines(as.numeric(names(curves[[s]])), curves[[s]], lwd = 2,
          col = cols_scn[s])
  }
}
plot.new()
legend("center", legend = labels, col = cols_scn, lwd = 2,
       bty = "n", cex = 0.85)

## --- Plot 2: Hospitalizations per 100,000 population by age ---
hosp_mat <- sapply(res, function(r) 1e5 * r$hosp / N)
tot <- colSums(hosp_mat)
par(mfrow = c(1, 1), mar = c(7, 4.5, 5, 1), mgp = c(3, 1, 0))
bp <- barplot(hosp_mat, names.arg = names(sims), las = 2,
              col = cols_age[age_groups],
              ylab = "Hospitalizations per 100,000 population",
              main = "Season Hospitalization Burden by Age",
              ylim = c(0, max(tot) * 1.15))
text(bp, tot + max(tot) * 0.04, sprintf("%.0f", tot), cex = 0.9, font = 2)
legend("top", legend = age_groups, horiz = TRUE, fill = cols_age[age_groups],
       bty = "n", cex = 0.9, inset = c(0, -0.18), xpd = TRUE)

## --- Plot 3: Hospitalizations averted and number needed to immunize ---
averted <- int_tbl$hosp_averted_per100k
names(averted) <- names(res)
short <- c(elderly_vax = "Older-adult vaccine", infant_proph = "Infant antibody",
           both = "Both products", cocoon = "Cocooning", npi = "NPI")
par(mfrow = c(1, 2), mar = c(9, 5, 3, 1), mgp = c(3.5, 1, 0))
bp2 <- barplot(averted[-1], names.arg = short[names(averted)[-1]],
               col = cols_scn[names(averted)[-1]], las = 2, cex.names = 0.8,
               ylab = "Averted per 100,000 population",
               main = "Averted in Target Groups",
               ylim = c(0, max(averted) * 1.25))
text(bp2, averted[-1] + max(averted) * 0.05, sprintf("%.1f", averted[-1]),
     cex = 0.85, font = 2)
keep <- which(!is.na(int_tbl$NNI))
nni <- int_tbl$NNI[keep]; names(nni) <- names(res)[keep]
bp3 <- barplot(nni, names.arg = short[names(nni)], col = cols_scn[names(nni)],
               las = 2, cex.names = 0.8, ylab = "Doses per hospitalization averted",
               main = "Number Needed to Immunize",
               ylim = c(0, max(nni) * 1.25))
text(bp3, nni + max(nni) * 0.05, nni, cex = 0.85, font = 2)
