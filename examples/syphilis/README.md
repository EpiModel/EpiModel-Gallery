# Syphilis Progression Model

## Description

This example models the transmission, multi-stage progression, diagnosis, and
treatment of syphilis using EpiModel's network-based framework. It extends the
basic SI model by adding multiple disease compartments representing the natural
history of syphilis, plus a treatment and screening module that allows infected
individuals to recover back to the susceptible state.

This is a simplified pedagogical model intended to demonstrate how to build a
multi-stage disease model in EpiModel with stage-dependent transmission,
symptom-driven diagnosis, and population-level screening.

## Model Structure

### Disease Stages

Syphilis progresses through six stages, each represented by a named value of the
`syph.stage` attribute:

| Stage | `syph.stage` | Infectious | Symptomatic | Duration (mean) |
|-------|-------------|-----------|-------------|-----------------|
| Incubating | `"incubating"` | Yes (high) | No | ~4 weeks |
| Primary | `"primary"` | Yes (high) | Possible | ~9 weeks |
| Secondary | `"secondary"` | Yes (high) | Possible | ~17 weeks |
| Early Latent | `"early_latent"` | Yes (low) | No | ~22 weeks |
| Late Latent | `"late_latent"` | No | No | ~29 years |
| Tertiary | `"tertiary"` | No | Yes | Terminal |

Susceptible individuals have `syph.stage = NA` and `status = "s"`.

### Flow Diagram

```mermaid
flowchart LR
    S((S)) -->|"infection"| Inc["Incubating"]
    Inc -->|"~4 wk"| Pri["Primary"]
    Pri -->|"~9 wk"| Sec["Secondary"]
    Sec -->|"~17 wk"| EL["Early Latent"]
    EL -->|"~22 wk"| LL["Late Latent"]
    LL -->|"~29 yr"| Ter["Tertiary"]

    Pri -.->|"symptomatic trt"| S
    Sec -.->|"symptomatic trt"| S
    Ter -.->|"symptomatic trt"| S
    Inc -.->|"screening"| S
    EL -.->|"screening"| S
    LL -.->|"screening"| S

    style Inc fill:#e74c3c,color:#fff
    style Pri fill:#e74c3c,color:#fff
    style Sec fill:#e74c3c,color:#fff
    style EL fill:#f39c12,color:#fff
    style LL fill:#95a5a6,color:#fff
    style Ter fill:#95a5a6,color:#fff
```

Solid arrows show disease progression; dotted arrows show treatment/screening
recovery back to susceptible. Red = infectious (high), orange = infectious
(reduced), gray = non-infectious.

### Transmission

Transmission probability depends on the infector's stage:

- **Incubating, primary, secondary** (`inf.prob.early`): higher transmission
  probability (0.18 per act)
- **Early latent** (`inf.prob.latent`): reduced transmission probability (0.09
  per act)
- **Late latent, tertiary**: not infectious (probability = 0)

The per-partnership transmission rate uses the standard EpiModel formula:
`finalProb = 1 - (1 - transProb)^actRate`.

### Treatment Pathways

There are two routes to treatment:

1. **Symptomatic diagnosis**: Individuals who develop symptoms during the primary
   or secondary stage may seek treatment (probability `early.trt` per timestep).
   Tertiary-stage patients, who are always symptomatic, may receive treatment at
   rate `late.trt`.
2. **Population screening**: Asymptomatic infected individuals may be detected
   through routine screening (probability `scr.rate` per timestep).

Recovery times after treatment initiation:
- Primary and secondary: 1 week
- Screening-detected (all non-tertiary stages): 2 weeks
- Tertiary: 3 weeks

Upon recovery, individuals return to the susceptible state (`status = "s"`,
`syph.stage = NA`).

## Modules

### Infection module (`infect`)

Extends the built-in EpiModel infection module with stage-dependent transmission
probabilities. Also handles initialization of all custom attributes at `at == 2`
(consolidated in one place for clarity). Key custom attributes:

- `syph.stage`: current disease stage (named string or `NA` if susceptible)
- `syph.symp`: symptomatic indicator (0 or 1)
- `infTime`, `priTime`, `secTime`, `elTime`, `llTime`, `terTime`: timestamps for
  stage transitions
- `syph.trt`, `syph.scr`: treatment and screening indicators
- `trtTime`, `scrTime`: treatment and screening timestamps

### Progression module (`progress`)

Simulates transitions between syphilis stages. Each transition is a stochastic
Bernoulli process with a constant hazard rate, so time spent in each compartment
follows a geometric distribution. A minimum 1-timestep delay is enforced before
any stage transition (e.g., `infTime < at` for incubating-to-primary).

Symptomatic status is assigned stochastically at the primary and secondary stages
(probabilities `pri.sym` and `sec.sym`). Tertiary-stage individuals are always
symptomatic.

### Treatment and screening module (`tnt`)

Organized into four sequential phases:

1. **Symptomatic treatment initiation**: Symptomatic primary, secondary, and
   tertiary patients may begin treatment each timestep.
2. **Symptomatic treatment recovery**: Treated patients recover after a
   stage-dependent delay (1 week for early stages, 3 weeks for tertiary).
3. **Screening**: Asymptomatic infected individuals (only) are screened at a
   population-level rate.
4. **Screening-detected recovery**: Screened patients in non-tertiary stages
   recover after a 2-week delay.

All modified attributes (`status`, `syph.stage`, `syph.symp`, `syph.trt`,
`syph.scr`, `trtTime`, `scrTime`) are saved at the end of the function.

## Parameters

### Transmission
| Parameter | Description | Value |
|-----------|-------------|-------|
| `inf.prob.early` | Per-act transmission probability for incubating, primary, and secondary stages | 0.18 |
| `inf.prob.latent` | Per-act transmission probability for early latent stage | 0.09 |
| `act.rate` | Number of acts per partnership per timestep | 2 |

### Stage Progression
| Parameter | Description | Value |
|-----------|-------------|-------|
| `ipr.rate` | Incubating to primary transition rate | 1/4 (~4 week duration) |
| `prse.rate` | Primary to secondary transition rate | 1/9 (~9 week duration) |
| `seel.rate` | Secondary to early latent transition rate | 1/17 (~17 week duration) |
| `elll.rate` | Early latent to late latent transition rate | 1/22 (~22 week duration) |
| `llter.rate` | Late latent to tertiary transition rate | 1/1508 (~29 year duration) |

### Symptoms
| Parameter | Description | Value |
|-----------|-------------|-------|
| `pri.sym` | Probability of developing symptoms per week in primary stage | 0.205 |
| `sec.sym` | Probability of developing symptoms per week in secondary stage | 0.106 |

### Treatment and Screening
| Parameter | Description | Value |
|-----------|-------------|-------|
| `early.trt` | Weekly probability of treatment given symptoms (primary/secondary) | 0.8 |
| `late.trt` | Weekly probability of treatment given symptoms (tertiary) | 1.0 |
| `scr.rate` | Weekly probability of screening (asymptomatic infected population) | 1/52 (~yearly) |

## Output Variables

| Variable | Description |
|----------|-------------|
| `s.num`, `i.num` | Susceptible and infected counts |
| `si.flow` | New infections per timestep |
| `inc.num`, `pr.num`, `se.num`, `el.num`, `ll.num`, `ter.num` | Counts by syphilis stage |
| `sym.num` | Count of symptomatic individuals |
| `ipr.flow`, `prse.flow`, `seel.flow`, `elll.flow`, `llter.flow` | Stage transition flows |
| `rec.flow` | Total recoveries per timestep (treatment + screening) |
| `scr.flow` | New screenings per timestep |
| `scr.num`, `trt.num` | Cumulative screened and currently on treatment |
| `syph.dur`, `syph2.dur`, ..., `syph6.dur` | Mean duration in each stage |

## Next Steps

Good next steps for this example might be:

- Add vital dynamics (births and deaths) to allow endemic equilibrium
- Vary progression rates by an individual attribute (e.g., age or immune status)
- Model reinfection with partial immunity after recovery
- Add partner notification as an additional pathway to treatment

## Authors

Samuel M. Jenness, Yuan Zhao, Emeli Anderson
