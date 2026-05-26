# SIR with Behavioral Risk Compensation During Illness

## Description

This example demonstrates a methodological lesson: when contact behavior changes over the course of an infection (people reduce their contacts when they feel most ill and return toward baseline as they recover), epidemic models that ignore this within-infection behavior change produce systematically biased estimates of per-act transmissibility and projected intervention impact.

Two structurally identical models are calibrated to the same cumulative attack rate and compared:

- **Naive model.** Contact rate is constant across the entire infectious period.
- **Dynamic model.** Contact rate is reduced during an early (most symptomatic) sub-stage of infection and partially recovers during a late sub-stage.

After calibration, both models are subjected to the same non-pharmaceutical intervention (NPI): an isolation policy that drives the early-stage contact rate to a low, fixed level (household-only contacts). The naive model substantially over-estimates the percent reduction in cumulative incidence, because its untouched baseline contact rate during the symptomatic period was implausibly high.

## Model Structure

### Disease Compartments

| Compartment | Label | Description |
|-------------|-------|-------------|
| Susceptible | **S** | Not infected; at risk of exposure |
| Infectious (early) | **I_e** | Infected, in the early sub-stage (most symptomatic) |
| Infectious (late) | **I_l** | Infected, in the late sub-stage (recovering) |
| Recovered | **R** | Recovered from infection; immune |

The two infectious sub-stages share the standard `status == "i"` label. A parallel `inf.stage` attribute carries `"early"` or `"late"` so that `discord_edgelist()` continues to work without modification and the sub-stage information can be used to look up the per-edge contact multiplier.

### Flow Diagram

```
S  --infection-->  I_e  --el.rate-->  I_l  --lr.rate-->  R
```

New infections always enter `I_e`. Progression is geometric in both sub-stages.

## Modules

### Infection Module (`infect`)

Standard discordant-edge transmission, with one change: the per-edge act count is multiplied by a stage-specific factor keyed to the infected partner's current sub-stage.

```r
stg  <- inf.stage[del$inf]
mult <- ifelse(stg == "early", mult.early, mult.late)

del$transProb <- inf.prob
del$actRate   <- act.rate * mult
del$finalProb <- 1 - (1 - del$transProb)^del$actRate
```

Setting `mult.early = mult.late = 1` recovers a behavior-naive SIR. Setting `mult.early < mult.late < 1` produces a dynamic model with reduced early-stage contacts that recover over time.

### Progression Module (`progress`)

Two stochastic stage transitions per timestep:

- `inf.stage = "early"` to `inf.stage = "late"` at rate `el.rate`.
- `inf.stage = "late"` and `status = "i"` to `status = "r"` at rate `lr.rate`.

Both modules idempotently initialize the `inf.stage` attribute on first call so the example does not depend on module ordering or on a custom `initialize.FUN`.

## Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `inf.prob` | Per-act transmission probability | Calibrated per model |
| `act.rate` | Acts per partnership per timestep | 1 |
| `el.rate` | Per-step rate of early-to-late transition | 1/3 (mean 3 days) |
| `lr.rate` | Per-step rate of late-to-recovered transition | 1/5 (mean 5 days) |
| `mult.early` | Act-rate multiplier when infected partner is in early sub-stage | 0.3 (dynamic), 1.0 (naive) |
| `mult.late` | Act-rate multiplier when infected partner is in late sub-stage | 0.6 (dynamic), 1.0 (naive) |

### Network

| Parameter | Description | Default (interactive) |
|-----------|-------------|-----------------------|
| Population size | Number of nodes | 500 |
| Target edges | Mean concurrent partnerships | 500 (mean degree = 2) |
| Partnership duration | Mean edge duration (timesteps) | 30 |

Treat each timestep as one day.

## Calibration

Both models are calibrated to the same cumulative attack rate (default target 40%) by bisection on `inf.prob`. Because the relationship between `inf.prob` and the final attack rate is monotone on a closed population, bisection converges reliably.

The expected result of calibration:

- The dynamic model needs a substantially higher `inf.prob` than the naive model to reproduce the same attack rate, because its infectious population spends part of the infectious period at reduced contact intensity.
- The ratio of calibrated `inf.prob` values is approximately the inverse of the time-weighted mean of the early and late multipliers. For the defaults: `(0.3 * 3 + 0.6 * 5) / 8 = 0.49`, so the dynamic-to-naive ratio is approximately `1 / 0.49 ~ 2.05`.

## NPI Scenario

Isolation of symptomatic cases is implemented by setting `mult.early` to a low absolute target (`iso.mult = 0.1`, representing household-only contacts during the most symptomatic phase). The late-stage multiplier is left unchanged.

Same intervention, different starting points:

- Naive baseline `mult.early` was 1.0, dropping to 0.1 is a 90% reduction.
- Dynamic baseline `mult.early` was 0.3, dropping to 0.1 is a 67% reduction.

The naive model therefore projects a much larger reduction in cumulative incidence than the dynamic model. The mechanism: the naive model gives itself more "room to intervene" by assuming a higher symptomatic-period baseline than the dynamic model considers realistic.

## Module Execution Order

EpiModel's default custom-module insertion places `progress.FUN` after the built-in `prevalence.FUN`. The example does not depend on a specific ordering: both modules check for and initialize the `inf.stage` attribute on first use.

## Next Steps

- Vary the multipliers to study sensitivity. Move `mult.early` from 0.1 to 1.0 to span "perfect isolation" to "no behavior change" and watch the calibrated `inf.prob` and NPI projections move together.
- Extend the within-infection structure: add a pre-symptomatic stage that is infectious but contact-unchanged (relevant for some respiratory pathogens with substantial pre-symptomatic transmission).
- Layer on a population-level behavioral response: scale `mult.early` and `mult.late` by an indicator of current prevalence to model fear-driven contact reduction in addition to symptom-driven reduction.
- Pair with the [Time-Varying Vaccination](../sir-time-varying-vaccination) example to study how the bias compounds when both the intervention and the underlying behavior vary over time.

## Author

Samuel M. Jenness, Emory University (http://samueljenness.org/)
