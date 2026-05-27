
##
## RSV: Age-Stratified SEIR Over a Multilayer Network
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##
## Authors: Samuel M. Jenness (Emory University)
## Date: May 2026
##

# Load EpiModel
suppressMessages(library(EpiModel))

# Standard Gallery unit test lines
rm(list = ls())
eval(parse(text = print(commandArgs(TRUE)[1])))

if (interactive()) {
  N <- 5000
  nsims <- 5
  ncores <- 5
  nsteps <- 120
} else {
  N <- 500
  nsims <- 1
  ncores <- 1
  nsteps <- 30
}


# 1. Population Generation ---------------------------------------------------

# Procedural generator for an age-stratified population. The age
# distribution is matched approximately to the United States by sampling
# household compositions and flattening them into per-person ages. The
# household structure supports a plausible joint age distribution (e.g.
# households containing infants always also contain adults), but the
# generated household_id is NOT used as an ERGM constraint and does not
# guarantee that any particular pair of nodes is connected in the family
# layer. The family layer is governed entirely by the nodemix age-mixing
# targets fit below.
generate_population <- function(N, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  hh_types <- list(
    list(p = 0.18, ages = c("adult")),
    list(p = 0.20, ages = c("adult", "adult")),
    list(p = 0.10, ages = c("elderly")),
    list(p = 0.10, ages = c("elderly", "elderly")),
    list(p = 0.04, ages = c("adult", "young", "school")),
    list(p = 0.03, ages = c("adult", "adult", "infant")),
    list(p = 0.05, ages = c("adult", "adult", "young")),
    list(p = 0.05, ages = c("adult", "adult", "young", "school")),
    list(p = 0.10, ages = c("adult", "adult", "school")),
    list(p = 0.05, ages = c("adult", "adult", "school", "school")),
    list(p = 0.05, ages = c("adult", "adult", "school", "elderly")),
    list(p = 0.05, ages = c("adult", "adult", "adult"))
  )
  probs <- sapply(hh_types, function(h) h$p)
  ages <- character(0); hh_id <- integer(0); cur <- 1L
  while (length(ages) < N) {
    idx <- sample.int(length(hh_types), 1, prob = probs)
    m <- hh_types[[idx]]$ages
    ages <- c(ages, m); hh_id <- c(hh_id, rep(cur, length(m))); cur <- cur + 1L
  }
  list(age = ages[1:N], household = hh_id[1:N])
}

set.seed(123)
pop <- generate_population(N, seed = 123)
# Alphabetical age order is what ergm uses internally:
counts <- table(factor(pop$age, levels = c("adult","elderly","infant","school","young")))
cat(sprintf("N = %d. Age counts:\n", N)); print(counts)


# 2. Network Setup -----------------------------------------------------------

# Two age-aware network layers. Both are estimated with nodemix() rather
# than nodefactor() so we can target specific cells of the age-by-age
# mixing matrix directly -- diagonal (within-age) and off-diagonal (cross-
# age) terms together. The cell ordering used by ergm's nodemix `levels2`
# is upper-triangle column-major. With factor levels alphabetical
# (adult, elderly, infant, school, young = 1..5), the canonical pair
# indices are:
#
#   1: (adult, adult)         9: (infant, school)
#   2: (adult, elderly)      10: (school, school)
#   3: (elderly, elderly)    11: (adult, young)
#   4: (adult, infant)       12: (elderly, young)
#   5: (elderly, infant)     13: (infant, young)
#   6: (infant, infant)      14: (school, young)
#   7: (adult, school)       15: (young, young)
#   8: (elderly, school)
#
# We target a sparse subset of cells in each layer (those that carry
# epidemiological signal) and leave the rest free, so the model is
# identifiable and ergm has room to fit without saturation.

nw <- network_initialize(N)
nw <- set_vertex_attribute(nw, "age", pop$age)

# Family layer: long-duration ties, low-degree, with ADULTS as hubs.
# Mean partnership duration equals the simulation horizon, but this is a
# mean, not a hard guarantee -- ties still resimulate each step.
# Constrained cells:
#   1 adult.adult      -- couples / adult-only households
#   4 adult.infant     -- parent-baby
#   7 adult.school     -- parent-child
#  11 adult.young      -- parent-toddler
#   3 elderly.elderly  -- elderly couples
#  10 school.school    -- siblings
fam_md <- c(adult = 3, elderly = 2, infant = 2, school = 3, young = 3)
fam_edges <- round(sum(counts * fam_md[names(counts)]) / 2)
fam_cell_targets <- c(
  adult.adult     = round(counts["adult"]   * 1.0 / 2),
  adult.infant    = round(counts["infant"]  * 2.0),
  adult.school    = round(counts["school"]  * 1.5),
  adult.young     = round(counts["young"]   * 1.7),
  elderly.elderly = round(counts["elderly"] * 1.0 / 2),
  school.school   = round(counts["school"]  * 0.5 / 2)
)
target_fam <- c(fam_edges, as.numeric(fam_cell_targets))
# `levels2` is inlined as a literal vector so the formula remains
# self-contained when ergm re-evaluates it inside netsim's resimulation.
formation_fam <- ~edges + nodemix("age", levels2 = c(1, 4, 7, 11, 3, 10))
coef.diss_fam <- dissolution_coefs(~offset(edges), duration = nsteps)

# Community layer: transient, high-degree, age-assortative.
# Constrained cells:
#   1 adult.adult, 3 elderly.elderly, 10 school.school, 15 young.young
#     (the four major diagonal homophily cells; infant.infant naturally 0)
#   7 adult.school -- a mild cross-tie (parents at school events, etc.)
com_md <- c(adult = 6, elderly = 3, infant = 1, school = 8, young = 5)
com_edges <- round(sum(counts * com_md[names(counts)]) / 2)
com_cell_targets <- c(
  adult.adult     = round(counts["adult"]   * 3.0 / 2),
  elderly.elderly = round(counts["elderly"] * 1.5 / 2),
  school.school   = round(counts["school"]  * 5.0 / 2),
  young.young     = round(counts["young"]   * 2.5 / 2),
  adult.school    = round(counts["school"]  * 0.5)
)
target_com <- c(com_edges, as.numeric(com_cell_targets))
formation_com <- ~edges + nodemix("age", levels2 = c(1, 3, 10, 15, 7))
coef.diss_com <- dissolution_coefs(~offset(edges), duration = 1)

cat(sprintf("\nFamily layer:    edges=%d, cells=%s\n",
            fam_edges, paste(fam_cell_targets, collapse = ",")))
cat(sprintf("Community layer: edges=%d, cells=%s\n",
            com_edges, paste(com_cell_targets, collapse = ",")))

cat("\nFitting family layer ERGM...\n")
est_fam <- netest(nw, formation_fam, target_fam, coef.diss_fam,
                  verbose = FALSE)
cat("Fitting community layer ERGM...\n")
est_com <- netest(nw, formation_com, target_com, coef.diss_com,
                  verbose = FALSE)


# 3. Disease + Intervention Parameters --------------------------------------

source("examples/rsv/module-fx.R")

# RSV natural-history parameters (daily timestep):
#   ei.rate    -- 1/4   mean 4-day latent period
#   ip.rate    -- 1/2   mean 2-day presymptomatic infectious period
#   ir.rate    -- 1/7   mean 7-day symptomatic/asymptomatic infectious period
#   asymp.prob -- 0.3   ~30% of infections never become symptomatic
#   asymp.inf.mult -- 0.5 asymptomatic half as infectious per contact
#
# Per-edge transmission probability (per day):
#   inf.prob.family    -- 0.05  intense close-contact transmission
#   inf.prob.community -- 0.018 less intense casual transmission
#
# Interventions (set per scenario). NOTE on intervention realism: each
# intervention here is implemented as a stylized reduction in
# susceptibility to infection. Real RSV vaccines and monoclonal antibody
# prophylaxis are evaluated against medically attended outcomes -- chiefly
# medically attended lower respiratory tract illness (LRTI),
# hospitalization, and severe disease -- not against infection itself.
# Using a per-edge susceptibility reduction is a pedagogical simplification
# that lets the same downstream pipeline (infections -> hospitalizations
# via age-specific risks) approximate the headline benefit.
#
#   elderly.vax.coverage / elderly.vax.efficacy -- leaky one-shot vaccine.
#     Arexvy / Abrysvo style. The 0.75 default loosely tracks reported
#     efficacy against RSV-associated LRTI in adults 60+. Real-world CDC
#     guidance (2024 update) recommends a single dose for adults 75+ and
#     for adults 50-74 at increased risk of severe disease; "all 65+" here
#     is a simplification.
#   infant.proph.coverage / infant.proph.efficacy -- passive antibody.
#     Nirsevimab style. The 0.70 default loosely tracks reported
#     efficacy against medically attended RSV LRTI in infants. Real-world
#     guidance covers either maternal vaccination or infant monoclonal
#     antibody (nirsevimab or clesrovimab); only the infant-side
#     prophylaxis is modeled here.
#   NPI window: community-layer per-edge inf.prob is reduced by
#     (1 - npi.mask.efficacy) and community edges are thinned by
#     npi.contact.mult, both only between npi.start and npi.end.
#
# References (CDC, accessed Nov 2025):
#   Adult RSV vaccine guidance:
#     https://www.cdc.gov/rsv/hcp/vaccine-clinical-guidance/adults.html
#   Infant/young child guidance:
#     https://www.cdc.gov/rsv/hcp/vaccine-clinical-guidance/infants-young-children.html
#   ACIP GRADE review (older adult RSV vaccines):
#     https://www.cdc.gov/acip/grade/protein-subunit-rsv-vaccines-older-adults.html
#   Nirsevimab MMWR recommendations:
#     https://www.cdc.gov/mmwr/volumes/72/wr/mm7234a4.htm

init <- init.net(i.num = round(0.005 * N))

# Base parameter set holds the defaults for every parameter the modules
# read. Per-scenario overrides are applied via the EpiModel scenarios
# API (create_scenario_list + use_scenario) further below. Intervention
# coverages default to zero and the NPI window defaults to inactive so
# the "none" scenario inherits an all-off baseline directly from this set.
param_base <- param.net(
  inf.prob.family = 0.05,
  inf.prob.community = 0.018,
  asymp.inf.mult = 0.5,
  ei.rate = 1 / 4,
  ip.rate = 1 / 2,
  ir.rate = 1 / 7,
  asymp.prob = 0.3,
  elderly.vax.coverage = 0,
  elderly.vax.efficacy = 0.75,
  infant.proph.coverage = 0,
  infant.proph.efficacy = 0.70,
  npi.start = -1,
  npi.end = -1,
  npi.mask.efficacy = 0.4,
  npi.contact.mult = 0.7
)

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
  verbose = FALSE
)


# 4. Scenarios --------------------------------------------------------------

# Scenarios are defined as a flat data frame (one row per scenario,
# columns matching parameter names in param_base) and converted to a
# scenario list with create_scenario_list(). Inside the loop,
# use_scenario(param_base, scn) returns a parameter object with the
# scenario's row applied to the base defaults.
scenarios.df <- data.frame(
  .scenario.id          = c("none", "elderly_vax", "infant_proph",
                            "both",  "npi"),
  .at                   = 0,
  elderly.vax.coverage  = c(0,       0.7,           0,
                            0.7,     0),
  infant.proph.coverage = c(0,       0,             0.8,
                            0.8,     0),
  npi.start             = c(-1,      -1,            -1,
                            -1,      20),
  npi.end               = c(-1,      -1,            -1,
                            -1,      80)
)
scenarios.list <- create_scenario_list(scenarios.df)

labels <- c(none = "No intervention",
            elderly_vax = "Elderly vaccination",
            infant_proph = "Infant prophylaxis",
            both = "Both (elderly + infant)",
            npi = "NPI (mask + distancing, days 20-80)")

sims <- list()
for (scn in scenarios.list) {
  cat(sprintf("Scenario: %s\n", scn$id))
  sims[[scn$id]] <- netsim(list(est_fam, est_com),
                           use_scenario(param_base, scn), init, control)
}


# 5. Hospitalization Analysis -----------------------------------------------

# Age-specific hospitalization risk per infection.
#
# These values are ILLUSTRATIVE -- chosen to reflect the qualitative
# RSV-NET pattern (highest in infants and older adults, lowest in
# school-age) without being directly fit to any one season's data. CDC
# RSV-NET reports age-specific RSV-associated hospitalization RATES per
# 100,000 population, not per infection; converting them to per-infection
# probabilities requires assumptions about season-specific infection
# attack rates by age that vary year to year. Use these as a teaching
# multiplier rather than a published estimate.
#
# References:
#   CDC RSV Surveillance (RSV-NET):
#     https://www.cdc.gov/rsv/research/rsv-net/dashboard.html
#   Hall CB et al. (2009). The burden of respiratory syncytial virus
#     infection in young children. N Engl J Med 360(6):588-598.
#   Falsey AR et al. (2005). Respiratory syncytial virus infection in
#     elderly and high-risk adults. N Engl J Med 352(17):1749-1759.
hosp_rate <- c(infant = 0.05, young = 0.008, school = 0.002,
               adult = 0.005, elderly = 0.030)

age_groups <- c("infant", "young", "school", "adult", "elderly")
age_pop <- as.numeric(counts[age_groups])
names(age_pop) <- age_groups

# Compute summary metrics for each scenario
cum_inf_by_scn <- list()
doses_by_scn <- numeric(0)
hosp_by_scn <- list()
for (s in names(sims)) {
  df <- as.data.frame(sims[[s]])
  last_t <- max(df$time)
  cum_per_age <- sapply(age_groups, function(a)
    mean(df[df$time == last_t, paste0("cuminf.", a)], na.rm = TRUE))
  cum_inf_by_scn[[s]] <- cum_per_age
  hosp_by_scn[[s]] <- cum_per_age * hosp_rate[age_groups]
  row <- scenarios.df[scenarios.df$.scenario.id == s, ]
  doses_by_scn[s] <- row$elderly.vax.coverage * age_pop["elderly"] +
                     row$infant.proph.coverage * age_pop["infant"]
}

# Summary table: per-age infections (hospitalizations in parens) by scenario
sum_tbl <- data.frame(Group = c(age_groups, "Total hospitalizations"),
                      stringsAsFactors = FALSE)
for (s in names(sims)) {
  cum <- cum_inf_by_scn[[s]]; hosp <- hosp_by_scn[[s]]
  sum_tbl[[s]] <- c(sprintf("%.0f (%.1f)", cum, hosp),
                    sprintf("%.1f", sum(hosp)))
}
cat("\n=== Cumulative infections (hospitalizations) by age ===\n")
print(sum_tbl)

# Per-dose efficiency
cat("\n=== Hospitalizations averted per dose ===\n")
none_hosp <- sum(hosp_by_scn$none)
eff_rows <- list()
for (s in names(sims)) {
  this_hosp <- sum(hosp_by_scn[[s]])
  averted <- none_hosp - this_hosp
  d <- doses_by_scn[s]
  cat(sprintf("%-22s doses=%6.0f  hosp=%6.1f  averted=%+6.1f  per-dose=%s\n",
              labels[s], d, this_hosp, averted,
              if (d > 0) sprintf("%.3f", averted / d) else "--"))
}


# 6. Plots ------------------------------------------------------------------

cols_scn <- c(none = "gray40", elderly_vax = "seagreen",
              infant_proph = "purple", both = "darkblue",
              npi = "firebrick")

## --- Plot 1: Age-stratified cumulative attack rates ---
par(mfrow = c(2, 3), mar = c(3, 3.5, 2, 1), mgp = c(2.2, 0.7, 0))
for (a in age_groups) {
  col_name <- paste0("cuminf.", a)
  max_y <- 0
  for (s in names(sims)) {
    df <- as.data.frame(sims[[s]])
    m <- tapply(df[[col_name]], df$time, mean, na.rm = TRUE) / age_pop[a]
    max_y <- max(max_y, max(m, na.rm = TRUE))
  }
  for (i in seq_along(sims)) {
    s <- names(sims)[i]
    df <- as.data.frame(sims[[s]])
    m <- tapply(df[[col_name]], df$time, mean, na.rm = TRUE) / age_pop[a]
    if (i == 1) {
      plot(as.numeric(names(m)), m, type = "l", lwd = 2,
           col = cols_scn[s], xlab = "Day",
           ylab = "Cumulative attack rate",
           main = paste0(toupper(substr(a, 1, 1)),
                         substr(a, 2, nchar(a)),
                         " (N=", age_pop[a], ")"),
           ylim = c(0, max_y * 1.05))
    } else {
      lines(as.numeric(names(m)), m, lwd = 2, col = cols_scn[s])
    }
  }
}
plot.new()
legend("center", legend = labels, col = cols_scn, lwd = 2,
       bty = "n", cex = 0.95)

## --- Plot 2: Stacked hospitalization burden by age ---
hosp_mat <- sapply(hosp_by_scn, function(h) h)
tot <- colSums(hosp_mat)
par(mfrow = c(1, 1), mar = c(7, 4, 5, 1), mgp = c(2.5, 1, 0))
bp <- barplot(hosp_mat, names.arg = names(sims), las = 2,
              col = c("#3498db", "#f39c12", "#e74c3c",
                      "#27ae60", "#8e44ad"),
              ylab = "Hospitalizations",
              main = "Stacked Hospitalizations by Age (lower = better)",
              ylim = c(0, max(tot) * 1.15))
text(bp, tot + max(tot) * 0.04, sprintf("%.1f", tot),
     cex = 0.9, font = 2)
legend("top", legend = age_groups, horiz = TRUE,
       fill = c("#3498db", "#f39c12", "#e74c3c", "#27ae60", "#8e44ad"),
       bty = "n", cex = 0.9,
       inset = c(0, -0.18), xpd = TRUE)
