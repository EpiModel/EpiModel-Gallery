# Age-Stratified RSV on a Household + Community Network

## Description

A single-season respiratory syncytial virus (RSV) model that combines three EpiModel capabilities no other Gallery example combines:

1. **A household plus community network.** The population is generated household by household from a table of household types, and every household is a fixed clique: a static contact layer with no ERGM and no resimulation, passed to the infection module as an edgelist parameter. On top of it runs a *community* TERGM layer of transient daily contacts in which every cell of the age-by-age mixing matrix is targeted through `nodemix`, with realized mean degree by age verified with `netdx`. Every transmission is recorded with `set_transmat()` together with its layer and both ages, so who-infects-whom is an output.

2. **A five-stratum age structure.** Infants (< 1 year), young children (1-4), school-age children (5-17), adults (18-64), and older adults (65+). Age governs household composition, the community contact profile, per-contact susceptibility (a proxy for prior exposure history), hospitalization risk per infection, and product eligibility. Transmission is seasonally forced with an annual cosine that peaks on day 1.

3. **Age-targeted and household-targeted interventions.** The older-adult vaccine and the infant monoclonal antibody each reduce the per-contact probability of infection (`eff.inf`, which also protects others) and the hospitalization risk given infection (`eff.hosp`, which protects only the recipient); the two combine to 80% effectiveness against hospitalization per exposure, the first-season value assumed by the RSV Scenario Modeling Hub, and the realized season-long effectiveness among recipients is computed from the simulation. Because households are explicit, a cocooning scenario that gives every co-resident of an infant a hypothetical transmission-blocking product is definable, and an equal-dose comparator that gives the same product to people of the same ages chosen at random isolates what the household link contributes. A community-layer NPI is included as a comparison.

The policy question: in one season, which strategy averts the most hospitalizations, how many doses does each hospitalization averted cost (number needed to immunize), and how much of each strategy's effect is indirect? Every intervention effect is reported with a 95% Monte Carlo interval from the between-simulation variance.

The annotated tutorial is on the [Gallery website](https://epimodel.github.io/EpiModel-Gallery/examples/rsv/).

## Model Structure

### Disease Compartments

| Compartment | Status | Description |
|-------------|--------|-------------|
| **S** | `"s"` | Susceptible (per-contact susceptibility scaled by age) |
| **E** | `"e"` | Exposed, latent (mean 4 days) |
| **I_p** | `"i"` + `inf_stage = "ip"` | Presymptomatic infectious (mean 2 days), for infections that will become symptomatic |
| **I_s** | `"i"` + `inf_stage = "is"` | Symptomatic infectious (mean 7 days) |
| **I_a** | `"i"` + `inf_stage = "ia"` | Asymptomatic infectious (mean 7 days, half as infectious); skips the presymptomatic stage |
| **R** | `"r"` | Recovered, immune for the rest of the season |

Whether an infection will be symptomatic is decided when it leaves E (30% asymptomatic). Every stage duration is geometric.

### Age Structure

Contact degrees are realized values at N = 10,000: the household degree follows from the household mix, the community degree from the fitted ERGM.

| Group | Ages | Share of N | Household degree | Community degree | `sus.mult` | Hospitalization risk per infection | Product |
|-------|------|------------|------------------|------------------|------------|-----------------------------------|---------|
| infant | < 1 | 1.2% | 2.8 | 2.3 | 1.00 | 0.030 | monoclonal antibody |
| young | 1-4 | 5.2% | 2.6 | 5.3 | 0.55 | 0.006 | |
| school | 5-17 | 18% | 2.7 | 6.9 | 0.16 | 0.001 | |
| adult | 18-64 | 58% | 1.8 | 5.0 | 0.07 | 0.004 | |
| elderly | 65+ | 18% | 0.9 | 3.1 | 0.13 | 0.045 | vaccine |

Anyone other than an infant who lives with an infant is eligible for the hypothetical cocooning product.

Hospitalization is not a compartment. Expected hospitalizations are computed after the simulation as infections times the per-infection risk, with the risk multiplied by `(1 - eff.hosp)` for immunized people.

## Network Layers

### Household layer (static cliques)

`hh_types` is a named vector: each name lists a household's members by age group and each value is the probability of that household type. `generate_households()` samples households until the population reaches `N` and returns each person's age and household id; `household_edgelist()` connects every pair of co-residents. The resulting edgelist is the whole household layer. It is set on the network as the vertex attribute `hh_id` (so `netsim` carries it to the modules) and passed to `param.net()` as `hh.pairs` (so the infection module can walk it each step). The mix gives a mean household size of 2.3, every infant at least one adult co-resident, about 60% of infants an older sibling, and about 28% of older adults living alone.

A clique layer is used instead of a long-duration ERGM layer because household transmission is closed: the people an infant can infect at home are exactly the people who can infect it. A random-graph family layer with the right degree by age does not have that closure, and household-targeted strategies cannot be defined on it.

### Community layer (TERGM)

`~edges + nodemix("age", levels2 = -1)`, with all 15 mixing cells set from a per-person contact profile that a helper converts to rounded edge-count targets. Targeting the full matrix matters: cells left out of `nodemix` absorb whatever edge count remains from the `edges` target at a uniform per-dyad rate, and because the adult-by-elderly block has far more dyads than any other cross-age block, a sparse specification gives older adults among the highest degrees in the layer. With 14 targeted cells, ergm's default simulated annealing step can fail to match the targets exactly and falls back to slow MCMC estimation; `control.ergm(SAN = control.san(SAN.maxit = 20, SAN.nsteps = 2^21))` lets the dyad-independent model be fit by maximum pseudolikelihood in seconds. Community ties last one day, so the layer reproduces daily contact rates but not the persistence of school and workplace contacts; co-resident pairs are not excluded from it, which adds a negligible number of edges.

## Modules

- `init_attrs`: one-shot setup of `vax_status` (`NA`, `"elderly_vax"`, `"infant_proph"`, `"cocoon"`) and `inf_stage`; places `init.net()` seeds in an infectious substage. Cocooning targets everyone whose `hh_id` matches an infant's; with `cocoon.random = 1` the same number of doses goes to people of the same ages drawn at random.
- `infect`: walks the household edgelist (from `hh.pairs`) and the community edgelist (from `get_edgelist()`) with the same code, applying the layer's transmission probability, the seasonal multiplier, the asymptomatic multiplier, the NPI factors on the community layer, `sus.mult` for the susceptible partner's age, and `eff.inf` for immunized susceptibles. Successful exposures from both layers are pooled, shuffled, and resolved to one infector per newly infected node (the tie-breaking rule of EpiModel's built-in module), then recorded with `set_transmat()` with the layer, both ages, and the infector's infection time.
- `progress`: E to I_p to I_s to R, or E to I_a to R, plus cumulative incident infections by age and among immunized infants and older adults (`cuminf.*`, `cuminf.*.prot`, `n.*.prot`, `n.cocoon`), excluding seeds.

`module.order` is set explicitly in `control.net()` so that infection runs before progression within each step; EpiModel's default would run the user modules first.

## Parameters

| Parameter | Value | Note |
|-----------|-------|------|
| `inf.prob.household` | 0.35 | annual-mean per-contact, per-day probability for a never-infected infant |
| `inf.prob.community` | 0.08 | annual-mean per-contact, per-day probability for a never-infected infant |
| `seas.amp`, `seas.peak` | 0.5, 1 | annual cosine forcing on both layers, peak on day 1 |
| `sus.mult` | 1.00 / 0.55 / 0.16 / 0.07 / 0.13 | infant / young / school / adult / elderly |
| `ei.rate`, `ip.rate`, `ir.rate` | 1/4, 1/2, 1/7 | daily |
| `asymp.prob`, `asymp.inf.mult` | 0.3, 0.5 | |
| `elderly.vax.eff.inf`, `elderly.vax.eff.hosp` | 0.5, 0.6 | combine to 0.80 against hospitalization per exposure |
| `infant.proph.eff.inf`, `infant.proph.eff.hosp` | 0.3, 0.71 | combine to 0.80 against hospitalization per exposure |
| `cocoon.eff.inf` | 0.5 | hypothetical adult product, no severity component |
| `npi.mask.efficacy`, `npi.contact.mult` | 0.4, 0.7 | community layer only, during the NPI window |
| `hh.pairs` | integer matrix | the household edgelist |

The transmission probabilities, forcing amplitude, and `sus.mult` were chosen together so the baseline season lands near the published age gradient of seasonal attack rates (50-70% of infants, 40-60% of 1-4 year olds, 20-30% of school-age children, 7-10% of adults, 3-7% of older adults). All values are illustrative, not fitted. The tutorial page has a provenance table with the basis and status of every parameter.

## Scenarios

| Scenario | Coverage or window |
|----------|--------------------|
| `none` | no intervention |
| `elderly_vax` | 50% of adults 65+ |
| `infant_proph` | 60% of infants |
| `both` | both products |
| `cocoon` | every co-resident of an infant (an upper bound on household targeting) |
| `cocoon_random` | the same number of doses of the cocooning product, to people of the same ages chosen at random |
| `npi` | days 30-90, community contacts cut 30% and remaining contacts transmit at 60% |

Product coverage follows the "usual" 2026-27 assumptions of the RSV Scenario Modeling Hub (Round 4). Eligibility is simplified relative to CDC guidance (nirsevimab or maternal vaccination for infants entering their first season; vaccine for adults 75+ and adults 50-74 at increased risk). No cocooning product is currently recommended for RSV; the scenario is the question an explicit household layer makes answerable.

## Outputs

`model.R` prints the household size distribution, realized mean degree by age and layer, cumulative attack rates by age through the last day, expected hospitalizations per 100,000 by age, the share of infections still in progress at the last day, hospitalizations averted in the target group and in all ages with 95% Monte Carlo intervals, doses and the number needed to immunize, the realized effectiveness of each product among its recipients, attack rates among unimmunized infants and older adults (the indirect effect), the household share of infant infections, the age-by-age transmission matrix and the mean number of secondary infections per seed from the transmission record, and the sources of infant infections by infector age and layer. Plots: age-stratified cumulative attack rates by scenario, hospitalizations per 100,000 stacked by age, and hospitalizations averted with intervals and NNI.

At N = 10,000 with ten simulations the baseline season produces about 100 hospitalizations per 100,000, within the range RSV-NET reports. The older-adult vaccine averts the most hospitalizations in absolute terms; the infant antibody has the lowest NNI. Neither product changes the attack rate among unimmunized people by more than the intervals allow. Even a complete household cocoon averts few infant hospitalizations, because more than half of infant infections come from outside the household and the leaky product protects co-residents by less than its per-contact value over a season of repeated exposure; the random arm with the same doses shows no infant effect. The full run takes about three minutes on a laptop.

## Running the Scripts

`model.R` sources `module-fx.R` from `examples/rsv/` when run from the repository root and from the working directory otherwise. Sourcing it interactively uses the full settings (`N = 10000`, ten simulations, 150 days). `Rscript` uses the small CI settings (`N = 1000`, one simulation, 50 days) whose results are not meant to be interpreted; pass `"run_full <- TRUE"` as the first argument for the full settings:

```bash
Rscript examples/rsv/model.R "run_full <- TRUE"
```

## Relationship to Other RSV Scenario Tools

The [RSV Scenario Modeling Hub](https://rsvscenariomodelinghub.org/) and [R.Scenario.Vax](https://chelsea-hansen.github.io/R.Scenario.Vax/) use age-structured compartmental models calibrated to RSV-NET, with births, maternal immunity, exposure-history immunity, seasonal forcing, and waning across seasons. This example trades those for an explicit contact network with household and individual-level targeting, and is not calibrated. It is a template for building an RSV network model, not a forecasting tool.

## Next Steps

- Births, infant age in months, and maternal vaccination (arrivals module, with new arrivals wired into an existing household).
- Contact persistence: a longer community tie duration at the same mean degree.
- Age-dependent natural history (infectious period and asymptomatic fraction by age).
- Weighted household ties (parent-infant versus sibling-infant contact).
- Waning of natural and product immunity, and multi-season runs.
- Seasonality estimated from data rather than an illustrative cosine.
- Calibration to RSV-NET (issue #58).
- Exposure-history immunity tracked per node rather than by age.
- Cross-layer dependency via `dat.updates` (see the SISMID multilayer tutorial).

## References

See the References section of the [tutorial page](https://epimodel.github.io/EpiModel-Gallery/examples/rsv/) for CDC guidance, the product trials (Hammitt 2022, Simões 2023, Papi 2023, Walsh 2023), the epidemiological sources (Glezen 1986, Hall 1976, Hall 2001, Falsey 2005, Hall 2009, Lessler 2009, Munywoki 2014, Pitzer 2015, Mossong 2008), and the RSV scenario modeling tools (Hansen 2025).

## Author

Samuel M. Jenness, Emory University (http://samueljenness.org/)
