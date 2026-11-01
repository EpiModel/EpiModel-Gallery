# Norovirus in a Nursing Home

## Description

A norovirus model of a 96-bed US nursing home over one winter season, built to ask how much each outbreak control measure prevents when transmission runs both from person to person and through contaminated surfaces, and whether the answer depends on how transmission divides between the two routes. It combines three features no other Gallery example combines:

1. **A facility built from clique layers and a shift roster.** Each time step is one 8-hour shift. Residents live in double and single rooms on three units of 32 beds (unit 1 a short-stay rehabilitation unit), and 77 staff work on a weekly roster with two days off a week: nursing assistants (CNAs) assigned to care groups of beds on each shift (4 per unit on days, 3 on evenings, 2 on nights), licensed nurses on each unit, and float staff who work days across all units. Roommates, residents of the same unit (dining and common areas), and staff on the same team are three `netclique()` layers with no ERGM fit; care contacts between residents and their assigned CNA and nurse are built inside the infection module each shift from the roster, so that absences, relief staff, and exclusions change who cares for whom. Residents are discharged and admitted into free beds; new residents join their room and unit through each layer's `arrivals.FUN`.

2. **An environmental reservoir on each unit.** Each unit carries a contamination level, kept between steps in `dat$run`, that infectious people add to (more if they vomit), that decays at 0.18 a day, and that cleaning reduces in scheduled pulses: routine cleaning once a day removes 30%, and hypochlorite during an outbreak with enhanced cleaning removes 60% twice a day. Everyone on the unit picks up from it in proportion to the contamination per person on the unit.

3. **Outbreak control that responds to the facility's own case count.** An outbreak is declared 24 hours after 3 or more people have symptom onset within 72 hours and is declared over after 4 days without a new onset. While it is declared, any of six measures can apply: isolating ill residents, excluding ill staff until 48 hours after symptoms end, cohorting float staff to one unit, closing affected units to admissions, hypochlorite cleaning, and soap and water after care (which cuts what staff hands carry in care and to and from surfaces).

The page has three parts. *Part 1* runs a season with no control measures and checks it against US and European outbreak data. *Part 2* compares each measure, the five main measures together, the same bundle declared earlier, and paid sick leave. *Part 3* recalibrates the model at five environmental shares of transmission, from 10% to 60%, each of which reproduces the same outbreaks, and compares the measures across them.

The annotated tutorial, with the results, is on the [Gallery website](https://epimodel.github.io/EpiModel-Gallery/examples/norovirus-nursing-home/).

## Model Structure

| Status | Description |
|--------|-------------|
| `s` | Susceptible |
| `e` | Infected, incubating (log-normal, median 1.2 days); not infectious |
| `i` | Symptomatic (mean 2 days); 60% vomit, with 1.5 times the infectiousness and shedding |
| `a` | Shedding without symptoms: after symptoms end, or for the 30% who never have them (mean 10 days); 5% as infectious as `i` for residents, 1% for staff |
| `r` | Immune for the rest of the season; 45% of people at the start of the season and on admission |

For an exposure `k` of a susceptible person in a shift, the hazard is `beta * w_k * infness` for contact with an infectious person, where `w_k` is the contact weight for that kind of contact on that shift, and `w_env * E_u / N_u` for the reservoir of unit `u`. Each exposure draws an exponential waiting time; the first one inside the shift, if any, infects and is recorded as the route.

## Modules

| Module | Purpose |
|--------|---------|
| `staffing` | Who is on duty this shift, and which CNA and nurse care for each resident |
| `infect` | Transmission on the three clique layers, in care, by float staff, and from the unit reservoirs, plus scheduled introductions through staff; records person-to-person transmissions with `set_transmat()` |
| `environ` | Deposition, decay, and cleaning of each unit's reservoir |
| `progress` | Natural history; logs each symptom onset |
| `outbreak` | Declaration and end of outbreaks, isolation and exclusion, units closed to admissions |
| `discharge`, `admit` | Discharges by unit length of stay; admissions into free beds, including scheduled introductions by admission, which arrive at a random point in their infection; residents who arrive ill are isolated |

## Running

```bash
# CI mode (three 60-day seasons per scenario, under a minute)
Rscript examples/norovirus-nursing-home/model.R

# Full settings (1,200 seasons per scenario, 36 scenarios, under an hour on ten cores)
Rscript examples/norovirus-nursing-home/model.R "run_full <- TRUE"
```

## Notes

- Introductions are rare, so most seasons would have none. Each simulated season draws its introductions before it runs, conditional on at least one, and season-level results are weighted by the probability of at least one.
- Every season is a single-simulation `netsim()` call run with `future.apply::future_lapply()` and a fixed `future.seed`, so season `s` has the same introductions in every scenario and scenarios are compared in pairs.
- The parameters are illustrative and the model is not calibrated to any facility. The two transmission parameters and the introduction rate are set to reproduce published outbreak sizes, attack rates, and frequency. The introduction rate that reproduces Oregon's outbreak frequency implies fewer staff infections outside work than community incidence suggests; the page discusses why.
