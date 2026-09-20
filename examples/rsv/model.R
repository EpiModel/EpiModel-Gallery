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

# Run settings. The full settings are used when the script is run
# interactively (for example sourced in RStudio). Rscript is not interactive,
# so it uses the small CI settings unless the first command-line argument
# defines run_full, which the unit test line above evaluates as R code:
#   Rscript examples/rsv/model.R "run_full <- TRUE"
# The full run takes several minutes; CI mode runs in well under a minute
# and its results are not meant to be interpreted.
if (interactive() || exists("run_full")) {
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
# States (about 1.2% infants, 5% young, 18% school-age, 58% adults, 18%
# older adults, mean household size 2.3), every infant lives with at least
# one adult, about 60% of infants have an older sibling (in the United
# States about 60% of births are second or later births), about 10% of
# infants live with a single adult, and about 28% of older adults live
# alone.
hh_types <- c(
  "adult"                          = 0.150,
  "adult adult"                    = 0.182,
  "adult adult adult"              = 0.050,
  "elderly"                        = 0.120,
  "elderly elderly"                = 0.120,
  "adult elderly"                  = 0.030,
  "adult adult infant"             = 0.007,
  "adult adult infant young"       = 0.007,
  "adult adult infant school"      = 0.006,
  "adult adult infant young school" = 0.003,
  "adult adult infant elderly"     = 0.002,
  "adult infant"                   = 0.002,
  "adult infant young"             = 0.001,
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
# mix mostly with each other (the school amplifier). Infants have the
# fewest contacts; theirs are non-household caregivers and, for the share
# in child care, other infants and young children. Cross-age entries are
# written from the perspective of the smaller group, e.g. school.adult =
# 1.5 means each school-age child has 1.5 adult community contacts on
# average.
#
# Every community tie lasts one day, so classmates and coworkers are
# redrawn daily: the layer reproduces the daily contact RATE by age but not
# the persistence of real school and workplace contacts. This is the
# simplest TERGM specification and a consequential one; a longer tie
# duration at the same mean degree would make repeated contacts with the
# same infectious person more likely and spread infection through fewer
# households.
com_profile <- c(
  school.school   = 5.0,   # classmates
  adult.adult     = 4.0,   # workplace, social
  young.young     = 2.5,   # child care
  elderly.elderly = 1.5,
  school.adult    = 1.5,   # teachers, coaches, friends' parents
  young.adult     = 2.0,   # child care staff
  elderly.adult   = 1.2,
  school.elderly  = 0.3,   # grandparents
  young.school    = 0.5,
  infant.adult    = 1.0,   # non-household caregivers
  infant.young    = 0.6,   # child care
  infant.infant   = 0.3,   # child care
  infant.school   = 0.2,
  infant.elderly  = 0.2,
  young.elderly   = 0.3
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

# The ERGM does not exclude co-resident pairs from the community layer. The
# expected number of community edges that fall on a household pair is the
# number of household pairs times the per-dyad tie probability, which is a
# few edges per day among the 25,000 or so in the layer and negligible
# next to the household transmission those pairs already have.


# 3. Disease + Intervention Parameters --------------------------------------

# model.R is written to be run from the repository root; the fallback lets
# the downloaded script run from a directory that holds both files.
source(if (file.exists("examples/rsv/module-fx.R")) "examples/rsv/module-fx.R"
       else "module-fx.R")

# Natural history (daily time step; every stage duration is geometric with
# the stated mean, because each transition is a daily Bernoulli draw):
#   ei.rate        -- 1/4  mean 4-day latent period
#   ip.rate        -- 1/2  mean 2-day presymptomatic infectious period
#   ir.rate        -- 1/7  mean 7-day symptomatic or asymptomatic period
#   asymp.prob     -- 0.3  share of infections that never become symptomatic;
#                          these skip the presymptomatic stage and are
#                          infectious for a mean of 7 days
#   asymp.inf.mult -- 0.5  asymptomatic cases are half as infectious
#
# Seasonal forcing. Both layers' transmission probabilities are multiplied
# by 1 + seas.amp * cos(2 * pi * (at - seas.peak) / 365): transmissibility
# is at its seasonal maximum on day seas.peak (here day 1, so the run
# starts at the top of the season) and declines through the run. RSV
# seasons end because transmissibility falls, not only because
# susceptibles run out, and a forced epidemic also has a much narrower
# final-size distribution across simulations than one that sits near the
# epidemic threshold for the whole run. Without forcing (seas.amp = 0) the
# same attack rates require an epidemic that grows only barely, and the
# between-simulation CV of season totals roughly doubles.
#
# Prior immunity. Nearly everyone is infected with RSV by age 2 and
# reinfected throughout life, with each reinfection milder and less likely
# (Glezen et al. 1986; Pitzer et al. 2015). Age is used here as a proxy for
# exposure history: sus.mult scales the per-contact probability of
# infection for a susceptible node of each age relative to a never-infected
# infant. The per-contact probabilities (inf.prob.household,
# inf.prob.community), the forcing amplitude, and sus.mult were chosen
# together, by iterating the baseline scenario, so that the season-long
# attack rates land near the
# published age gradient: 50-70% of infants (Glezen et al. 1986 report 69%
# in the first year of life), 40-60% of 1-4 year olds, 20-30% of school-age
# children, about 7-10% of adults (7% among healthy working adults, Hall
# et al. 2001; higher for parents of young children), and 3-7% of older
# adults (Falsey et al. 2005). The baseline attack rates are printed below
# so the reader can check that the run in hand still lands there. They are
# illustrative, not fitted to data.
#
# Immunization products. Each has two leaky components, following how the
# products are evaluated and how other RSV scenario models represent them:
#   *.eff.inf  -- reduction in the per-contact probability of infection
#                 (this is the only component that produces indirect
#                 protection of others through reduced transmission)
#   *.eff.hosp -- reduction in the hospitalization risk given infection
# Against a single exposure the two combine to
# 1 - (1 - eff.inf) * (1 - eff.hosp) = 0.80 for both products, the
# first-season effectiveness against hospitalization assumed by the RSV
# Scenario Modeling Hub (80% for infant monoclonals and for the older-adult
# vaccines in the year of vaccination; Round 4, 2026-27 season). Because
# eff.inf is leaky and acts per contact, a recipient exposed repeatedly over
# the season is protected against infection by less than eff.inf, so the
# realized effectiveness against hospitalization over the season is an
# OUTPUT of the model, computed in the analysis from the attack rates among
# immunized and unimmunized members of each group, rather than an input.
# Trial estimates against medically attended RSV illness are 74-80% for
# nirsevimab in infants (Hammitt et al. 2022; Simoes et al. 2023) and
# 67-83% for the older-adult vaccines (Papi et al. 2023; Walsh et al.
# 2023). The split between the two components is an assumption; setting
# eff.inf = 0 and eff.hosp = 0.80 gives a product with the same
# per-exposure effectiveness, no indirect effect, and a realized
# effectiveness of exactly 80% by construction.
#
# Household targeting ("cocooning"). Because households are explicit, the
# co-residents of infants are a definable target group. The cocoon
# scenario gives every co-resident of an infant (parents, siblings,
# grandparents) a HYPOTHETICAL product with the older-adult vaccine's
# infection-blocking component and no severity component, to ask how much
# infant protection blocking household transmission can deliver, at most,
# compared with the direct infant product. A comparator scenario gives the
# same number of doses of the same product to people of the same ages
# chosen at random, which isolates what the household link itself
# contributes. No such product is currently recommended for this purpose;
# maternal vaccination protects the infant through antibody transfer, not
# by blocking the parent's transmission.
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
  inf.prob.household = 0.35,
  inf.prob.community = 0.08,
  seas.amp = 0.5,
  seas.peak = 1,
  sus.mult = c(infant = 1.00, young = 0.55, school = 0.16,
               adult = 0.07, elderly = 0.13),
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
  cocoon.random = 0,
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
# Product coverage follows the "usual" assumptions of the RSV Scenario
# Modeling Hub for the 2026-27 season (56% of infants receiving a
# monoclonal antibody, 50% of adults 75+ vaccinated). The cocoon scenario
# covers every co-resident of an infant, an upper bound on household
# targeting; the random arm is its equal-dose, age-matched comparator.
scenarios.df <- data.frame(
  .scenario.id          = c("none", "elderly_vax", "infant_proph", "both",
                            "cocoon", "cocoon_random", "npi"),
  .at                   = 0,
  elderly.vax.coverage  = c(0, 0.5, 0, 0.5, 0, 0, 0),
  infant.proph.coverage = c(0, 0, 0.6, 0.6, 0, 0, 0),
  cocoon.coverage       = c(0, 0, 0, 0, 1, 1, 0),
  cocoon.random         = c(0, 0, 0, 0, 0, 1, 0),
  npi.start             = c(-1, -1, -1, -1, -1, -1, 30),
  npi.end               = c(-1, -1, -1, -1, -1, -1, 90)
)
scenarios.list <- create_scenario_list(scenarios.df)

labels <- c(none = "No intervention",
            elderly_vax = "Older-adult vaccine (50%)",
            infant_proph = "Infant antibody (60%)",
            both = "Both products",
            cocoon = "Household cocoon (all co-residents of infants)",
            cocoon_random = "Same doses, random people of the same ages",
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
hosp_rate <- c(infant = 0.030, young = 0.006, school = 0.001,
               adult = 0.004, elderly = 0.045)
eff_hosp <- c(infant = param_base$infant.proph.eff.hosp,
              elderly = param_base$elderly.vax.eff.hosp)

age_groups <- c("infant", "young", "school", "adult", "elderly")
age_pop <- as.numeric(counts[age_groups])
names(age_pop) <- age_groups

# Season-end summary for one scenario. Per-simulation quantities are kept
# (one row per simulation) so that Monte Carlo intervals can be computed.
summarize_scenario <- function(sim) {
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
  # Attack rates among immunized and unimmunized infants and older adults,
  # per simulation. The realized effectiveness of a product among its
  # recipients is one minus their ratio.
  attack_prot <- attack_unprot <- matrix(NA_real_, nrow(last), 2,
                                         dimnames = list(NULL, c("infant", "elderly")))
  for (a in c("infant", "elderly")) {
    n_p <- last[[paste0("n.", a, ".prot")]]
    i_p <- last[[paste0("cuminf.", a, ".prot")]]
    attack_prot[, a] <- ifelse(n_p > 0, i_p / n_p, NA)
    attack_unprot[, a] <- (last[[paste0("cuminf.", a)]] - i_p) / (age_pop[a] - n_p)
  }
  # Doses actually delivered in the simulation: immunized older adults and
  # infants, plus people given the cocooning product.
  doses <- mean(last$n.elderly.prot + last$n.infant.prot + last$n.cocoon)
  # Share of infant infections acquired from household contacts, pooled
  # over the season and the simulations
  hh_share_infant <- sum(df$se.flow.infant.hh, na.rm = TRUE) /
    sum(df$se.flow.infant.hh + df$se.flow.infant.com, na.rm = TRUE)
  # Infections still in progress at the end of the run, as a check on the
  # observation window
  active_end <- mean(last$e.num + last$ip.num + last$is.num + last$ia.num)
  list(inf = inf, attack = inf / age_pop,
       hosp = colMeans(hosp_sim), hosp_sim = hosp_sim,
       hosp_per100k = 1e5 * colMeans(hosp_sim) / age_pop,
       attack_prot = attack_prot, attack_unprot = attack_unprot,
       doses = doses, hh_share_infant = hh_share_infant,
       active_end = active_end, cuminf_end = sum(inf))
}

res <- lapply(sims, summarize_scenario)

cat("\n=== Cumulative attack rate through day", nsteps, "by age (%) ===\n")
attack_tbl <- sapply(res, function(r) round(100 * r$attack, 1))
print(attack_tbl)

cat("\n=== Expected hospitalizations per 100,000 by age (total row: all ages) ===\n")
hosp100k_tbl <- rbind(sapply(res, function(r) round(r$hosp_per100k)),
                      total = sapply(res, function(r) round(1e5 * sum(r$hosp) / N)))
print(hosp100k_tbl)

# Observation window: infections still in progress at the last step, as a
# share of the season's cumulative infections. The scenario outcomes are
# "through day nsteps"; a large value here would mean the comparison is
# partly about timing rather than final size.
cat("\n=== Infections in progress at day", nsteps, "as % of cumulative infections ===\n")
print(sapply(res, function(r) round(100 * r$active_end / r$cuminf_end, 1)))

# Monte Carlo interval for a difference in means between two sets of
# simulations, from the between-simulation variance (NA with one simulation)
mc_diff <- function(x0, x1) {
  d <- mean(x0) - mean(x1)
  if (length(x0) < 2 || length(x1) < 2) return(c(est = d, lo = NA, hi = NA))
  se <- sqrt(var(x0) / length(x0) + var(x1) / length(x1))
  c(est = d, lo = d - 1.96 * se, hi = d + 1.96 * se)
}
fmt_ci <- function(v, digits = 1) {
  if (is.na(v["lo"])) return(sprintf("%.*f", digits, v["est"]))
  sprintf("%.*f (%.*f, %.*f)", digits, v["est"], digits, v["lo"], digits, v["hi"])
}

# Age groups each strategy is designed to protect. Hospitalizations averted
# and NNI are computed within these groups; the all-ages difference is
# reported next to them with its own interval. The between-simulation
# noise in the large untargeted strata is comparable to a product's whole
# effect, which the all-ages intervals make visible.
target_groups <- list(none = age_groups, elderly_vax = "elderly",
                      infant_proph = "infant", both = c("infant", "elderly"),
                      cocoon = "infant", cocoon_random = "infant",
                      npi = age_groups)
per100k <- function(r, groups) 1e5 * rowSums(r$hosp_sim[, groups, drop = FALSE]) / N
averted <- lapply(names(res), function(s) {
  tg <- target_groups[[s]]
  list(target = mc_diff(per100k(res$none, tg), per100k(res[[s]], tg)),
       all = mc_diff(per100k(res$none, age_groups), per100k(res[[s]], age_groups)),
       base_target = mean(per100k(res$none, tg)))
})
names(averted) <- names(res)

# NNI is doses per hospitalization averted in the target group. It is
# undefined when the point estimate of averted hospitalizations is not
# positive, and its interval is reported only when the whole interval for
# averted hospitalizations is positive.
nni_of <- function(s) {
  a <- averted[[s]]$target
  doses <- res[[s]]$doses
  averted_n <- a * N / 1e5                  # per 100,000 -> count in N
  if (doses == 0 || is.na(averted_n["est"]) || averted_n["est"] <= 0) return(NA_character_)
  est <- doses / averted_n["est"]
  if (is.na(averted_n["lo"]) || averted_n["lo"] <= 0) {
    return(sprintf("%.0f", est))
  }
  sprintf("%.0f (%.0f, %.0f)", est, doses / averted_n["hi"], doses / averted_n["lo"])
}

cat("\n=== Intervention summary: hospitalizations averted per 100,000 population (95% Monte Carlo interval) ===\n")
int_tbl <- data.frame(
  scenario = labels[names(res)],
  target = sapply(names(res), function(s)
    if (length(target_groups[[s]]) == 5) "all ages" else
      paste(target_groups[[s]], collapse = " + ")),
  doses = sapply(res, function(r) round(r$doses)),
  averted_target = sapply(names(res), function(s) fmt_ci(averted[[s]]$target)),
  pct_averted_target = sapply(names(res), function(s) {
    b <- averted[[s]]$base_target
    if (b > 0) round(100 * averted[[s]]$target["est"] / b, 1) else NA
  }),
  averted_all_ages = sapply(names(res), function(s) fmt_ci(averted[[s]]$all)),
  NNI = sapply(names(res), nni_of),
  row.names = NULL
)
print(int_tbl, right = FALSE)

# Realized effectiveness among recipients. The infection-blocking component
# is leaky and per contact, so the reduction in a recipient's season-long
# risk of infection is smaller than eff.inf when exposures are repeated;
# the realized effectiveness against hospitalization follows from it and
# the severity component. This is the number to compare with the 80%
# first-season effectiveness assumed by the Scenario Modeling Hub.
ve_tbl <- do.call(rbind, lapply(c("elderly_vax", "infant_proph"), function(s) {
  a <- if (s == "elderly_vax") "elderly" else "infant"
  r <- res[[s]]
  ve_inf <- ifelse(r$attack_unprot[, a] > 0,
                   1 - r$attack_prot[, a] / r$attack_unprot[, a], NA)
  ve_hosp <- 1 - (1 - ve_inf) * (1 - eff_hosp[a])
  ci <- function(v) {
    m <- mean(v, na.rm = TRUE)
    if (sum(!is.na(v)) < 2) return(sprintf("%.0f", 100 * m))
    se <- sd(v, na.rm = TRUE) / sqrt(sum(!is.na(v)))
    sprintf("%.0f (%.0f, %.0f)", 100 * m, 100 * (m - 1.96 * se), 100 * (m + 1.96 * se))
  }
  eff_inf <- if (a == "elderly") param_base$elderly.vax.eff.inf else param_base$infant.proph.eff.inf
  data.frame(product = labels[s], group = a,
             per_contact_eff_inf = 100 * eff_inf,
             realized_VE_infection = ci(ve_inf),
             eff_hosp = 100 * eff_hosp[a],
             per_exposure_VE_hosp = round(100 * (1 - (1 - eff_inf) * (1 - eff_hosp[a]))),
             realized_VE_hosp = ci(ve_hosp),
             row.names = NULL)
}))
cat("\n=== Product effectiveness among recipients (%): per-contact inputs and realized season-long values ===\n")
print(ve_tbl, right = FALSE)

# Indirect protection: attack rate among UNimmunized infants and older
# adults relative to the no-intervention baseline. Any reduction here is
# transmission blocked by the eff.inf component in other people. Under
# cocooning no infant is immunized, so the infant row is the whole effect.
cat("\n=== Attack rate among unimmunized (%), mean (95% Monte Carlo interval) ===\n")
unprot_tbl <- sapply(res, function(r) sapply(c("infant", "elderly"), function(a) {
  v <- 100 * r$attack_unprot[, a]
  if (length(v) < 2) return(sprintf("%.1f", mean(v)))
  se <- sd(v) / sqrt(length(v))
  sprintf("%.1f (%.1f, %.1f)", mean(v), mean(v) - 1.96 * se, mean(v) + 1.96 * se)
}))
print(unprot_tbl, quote = FALSE)

# Where infants get infected: the share of infant infections acquired from
# household contacts. This is what limits cocooning, which blocks only the
# adult-to-infant part of that share.
cat("\n=== Share of infant infections acquired from household contacts (%) ===\n")
print(sapply(res, function(r) round(100 * r$hh_share_infant, 1)))

# Who infects whom. Every transmission was recorded with set_transmat(),
# with the infector's and recipient's age groups and the layer. Pooled
# over the baseline simulations, this gives the age-by-age transmission
# matrix, the layer split, and an estimate of the reproduction number from
# the seeds: the mean number of secondary infections generated by the
# initial infections, which were placed at random across ages.
tm_none <- do.call(rbind, lapply(seq_len(nsims), function(s) {
  as.data.frame(get_transmat(sims$none, sim = s))
}))
tm_none$infAge <- factor(tm_none$infAge, levels = age_groups)
tm_none$susAge <- factor(tm_none$susAge, levels = age_groups)
waifw <- round(100 * prop.table(table(infector = tm_none$infAge,
                                      recipient = tm_none$susAge)), 1)
cat("\n=== Who infects whom (% of all baseline infections; rows infector, columns recipient) ===\n")
print(waifw)
cat("\nShare of baseline infections caused by each age group (%):\n")
print(round(100 * prop.table(table(tm_none$infAge)), 1))
cat(sprintf("\nShare of baseline infections on the household layer: %.1f%%\n",
            100 * mean(tm_none$layer == 1)))
n_seeds <- round(0.01 * N)
cat(sprintf("Mean secondary infections per seed infection: %.2f\n",
            sum(tm_none$infTime == 1) / (n_seeds * nsims)))

# Infant infections by infector age and layer: the pathway question
inf_tm <- tm_none[tm_none$susAge == "infant", ]
infant_src <- round(100 * prop.table(table(
  infector = inf_tm$infAge,
  layer = factor(inf_tm$layer, levels = 1:2, labels = c("household", "community")))), 1)
cat("\n=== Infant infections by infector age and layer (% of infant infections, baseline) ===\n")
print(infant_src)

# Between-simulation variability: range of the all-ages total across
# simulations, and the CV of season-end cumulative infections by age. The infant stratum is small (about 1% of
# N), so it is the noisiest; this is the reason for the large default N.
hosp_total_sims <- do.call(cbind, lapply(res, function(r) per100k(r, age_groups)))
cat("\n=== All-ages hospitalizations per 100,000: mean and range across simulations ===\n")
print(round(rbind(mean = colMeans(hosp_total_sims),
                  min = apply(hosp_total_sims, 2, min),
                  max = apply(hosp_total_sims, 2, max))))
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
              cocoon = "darkorange", cocoon_random = "goldenrod",
              npi = "firebrick")
lty_scn <- c(none = 1, elderly_vax = 1, infant_proph = 1, both = 1,
             cocoon = 1, cocoon_random = 2, npi = 1)
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
          col = cols_scn[s], lty = lty_scn[s])
  }
}
plot.new()
legend("center", legend = labels, col = cols_scn, lwd = 2, lty = lty_scn,
       bty = "n", cex = 0.8)

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

## --- Plot 3: Hospitalizations averted (with Monte Carlo intervals) and NNI ---
scn_int <- setdiff(names(res), "none")
av_est <- sapply(scn_int, function(s) averted[[s]]$target["est"])
av_lo <- sapply(scn_int, function(s) averted[[s]]$target["lo"])
av_hi <- sapply(scn_int, function(s) averted[[s]]$target["hi"])
short <- c(elderly_vax = "Older-adult vaccine", infant_proph = "Infant antibody",
           both = "Both products", cocoon = "Household cocoon",
           cocoon_random = "Random, same doses", npi = "NPI")
par(mfrow = c(1, 2), mar = c(9, 5, 3, 1), mgp = c(3.5, 1, 0))
yr <- range(c(0, av_est, av_lo, av_hi), na.rm = TRUE)
bp2 <- barplot(av_est, names.arg = short[scn_int], col = cols_scn[scn_int],
               las = 2, cex.names = 0.8,
               ylab = "Averted per 100,000 population",
               main = "Averted in Target Groups",
               ylim = yr + diff(yr) * c(-0.05, 0.15))
abline(h = 0)
if (!all(is.na(av_lo))) {
  arrows(bp2, av_lo, bp2, av_hi, angle = 90, code = 3, length = 0.04)
}
text(bp2, pmax(av_hi, av_est, na.rm = TRUE) + diff(yr) * 0.04,
     sprintf("%.1f", av_est), cex = 0.85, font = 2)
nni_num <- sapply(scn_int, function(s) {
  a <- averted[[s]]$target["est"] * N / 1e5
  if (res[[s]]$doses > 0 && !is.na(a) && a > 0) res[[s]]$doses / a else NA
})
keep <- which(!is.na(nni_num))
if (length(keep) > 0) {
  bp3 <- barplot(nni_num[keep], names.arg = short[scn_int[keep]],
                 col = cols_scn[scn_int[keep]], las = 2, cex.names = 0.8,
                 ylab = "Doses per hospitalization averted",
                 main = "Number Needed to Immunize",
                 ylim = c(0, max(nni_num[keep]) * 1.25))
  text(bp3, nni_num[keep] + max(nni_num[keep]) * 0.05,
       sprintf("%.0f", nni_num[keep]), cex = 0.85, font = 2)
}
