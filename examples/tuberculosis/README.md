# Tuberculosis in Households and the Community

## Description

A tuberculosis (TB) model on a three-layer contact network, built to ask why households carry a small share of TB transmission in high-burden settings despite holding a large share of contact time, and what that means for household contact investigation. It combines three features no other Gallery example combines:

1. **A household clique layer in an open population.** Households are cliques built with `netclique()`. Newborns join the household of an adult aged 15 to 49 through the layer's arrival rule, the dead leave with their edges, and young adults leave home to form new households through `move_to_group()`. Two TERGM layers carry non-household contacts: *regular* contacts with a mean duration of 12 months and *casual* contacts redrawn every month. Each layer's edges carry the contact time adults in the Western Cape spend with household members (17.0 hours a day), repeated contacts (6.0), and non-repeated contacts (24.2) (McCreesh and White 2018).

2. **A TB natural history with age structure.** Adapted from the compartmental model of Rothman et al. (2026): early and late latent infection, reinfection with 79% protection from prior infection, clearance of infection, active TB with diagnosis, self-cure, and death, a six-month treatment course, and relapse. Progression, infectiousness, and TB mortality differ for ages 0 to 4, 5 to 14, and 15 and older, and each case's infectiousness is drawn from a gamma distribution with shape 0.15, so that about 20% of cases cause about 90% of transmission.

3. **Two experiments.** *Part 1* holds contact time and per-contact-hour risk fixed and compares three arrangements of the same contacts: households as they are, regular ties redrawn every month, and household contact time spread over random partners redrawn every month. *Part 2* resumes the baseline and compares household contact investigation without preventive therapy (TPT), with TPT for contacts under 5, and with TPT for contacts of all ages; community-wide screening of adults (ACT3-like, 42% of adults tested each year); and an equal-effort comparator that screens as many randomly chosen people as the household program. The three household arms differ only in who is offered TPT, which isolates its direct effect on contacts.

The annotated tutorial, with the results, is on the [Gallery website](https://epimodel.github.io/EpiModel-Gallery/examples/tuberculosis/).

## Model Structure

| Status | Description |
|--------|-------------|
| `s` | Never infected |
| `e` | Early latent infection (months after infection, highest risk of disease) |
| `l` | Late latent infection (small lifelong risk of reactivation; 21% susceptibility to reinfection) |
| `c` | Infection cleared, spontaneously or by TPT (no reactivation; 21% susceptibility to reinfection) |
| `i` | Active TB, not yet diagnosed (infectious by age and case) |
| `t` | On treatment for six months (not infectious) |
| `r` | Recovered from TB (relapse risk; full susceptibility to reinfection) |

A state with more than one exit uses competing risks at the monthly step. The monthly probability that an edge on layer `k` transmits is `1 - exp(-beta * layer.hours[k] * inf_mult * age_mult * sus)`.

## Modules

| Module | Purpose |
|--------|---------|
| `init_tb` | One-shot setup of TB states by age from a constant annual risk of infection |
| `aging` | Advances age by a month; updates the child/adult attribute the TERGMs mix on |
| `infect` | Transmission on every layer weighted by contact-hours; records every transmission with `set_transmat()` |
| `progress` | Natural history and passive diagnosis with competing risks |
| `screen` | Household contact investigation, community-wide screening, and random screening; records the state of every household contact at each diagnosis |
| `deaths`, `births` | Age-specific background mortality; births at a constant per-capita rate |
| `households` | Leaving home and moving children without an adult, with `move_to_group()` |
| `tally` | Counts by TB state and age, and household structure |

`newborn_household()` is the `arrivals.FUN` for the clique layer.

## Running

```bash
# CI mode (small population, about a minute)
Rscript examples/tuberculosis/model.R

# Full settings (N = 20,000, ten simulations, about half an hour on five cores)
Rscript examples/tuberculosis/model.R "run_full <- TRUE"
```

## Notes

- `tergm`'s default Markov chain per time step is far too short for these layers: about half of the one-month casual ties carried over between months, and the regular layer's children-with-children cell fell about 20% below its target. The example sets `MCMC.burnin.min` and `MCMC.burnin.max` in proportion to `N`, per layer, with `multilayer()`.
- The parameters are illustrative and the model is not calibrated to any country. The setting is a generic high-burden one without HIV.
