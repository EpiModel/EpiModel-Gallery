# SEIR/SEIRS Model: Adding an Exposed State to an SIR

## Description

This example demonstrates how to extend EpiModel's built-in SIR model by adding an **Exposed (E)** compartment, creating an SEIR model. The E compartment represents a latent period during which a person has been infected but is not yet infectious to others. Many infectious diseases have this incubation stage (e.g., influenza, measles, COVID-19), making SEIR one of the most widely used compartmental frameworks in infectious disease epidemiology.

The example also includes an **SEIRS extension** that adds waning immunity: recovered individuals eventually lose protection and return to the susceptible pool. This single change fundamentally alters the long-run dynamics — transforming an epidemic that burns out into one that persists endemically.

This is the foundational "custom compartment" example in the Gallery. The techniques shown here — replacing the built-in infection module, writing a custom progression module, and tracking new compartment counts — are the building blocks for all the more complex models in the collection.

## Model Structure

### Disease Compartments

| Compartment | Label | Description |
|-------------|-------|-------------|
| Susceptible | **S** | Not infected; at risk of exposure |
| Exposed | **E** | Infected but not yet infectious (latent period) |
| Infectious | **I** | Infected and capable of transmitting to contacts |
| Recovered | **R** | Immune after recovery; permanent (SEIR) or temporary (SEIRS) |

### Flow Diagram

```mermaid
flowchart LR
    S["<b>S</b><br/>Susceptible"] -->|"infection<br/>(se.flow)"| E
    E["<b>E</b><br/>Exposed"] -->|"progression<br/>(ei.flow)"| I
    I["<b>I</b><br/>Infectious"] -->|"recovery<br/>(ir.flow)"| R
    R["<b>R</b><br/>Recovered"] -.->|"waning immunity<br/>(rs.flow, SEIRS only)"| S

    style S fill:#3498db,color:#fff
    style E fill:#f39c12,color:#fff
    style I fill:#e74c3c,color:#fff
    style R fill:#27ae60,color:#fff
```

### Transmission

Transmission occurs along **discordant edges** — partnerships where one node is susceptible and the other is infectious. For each discordant pair, the per-timestep probability of transmission accounts for multiple acts per partnership:

```
finalProb = 1 - (1 - inf.prob) ^ act.rate
```

With the default parameters (`inf.prob = 0.5`, `act.rate = 2`), this gives a per-partnership transmission probability of 0.75 per timestep.

**Key SEIR difference:** Newly infected individuals enter the Exposed (`"e"`) state rather than becoming immediately Infectious (`"i"`), as they would in a standard SIR model.

### Disease Progression

Each transition (E→I, I→R, R→S) is modeled as an independent **Bernoulli trial** at each timestep: each eligible individual has a constant probability of transitioning, regardless of how long they have been in their current compartment. This produces **geometrically distributed** waiting times with mean `1/rate`.

## Modules

### Infection Module (`infect`)

Replaces EpiModel's built-in infection module. Extracts the discordant edgelist via `discord_edgelist()`, computes per-partnership transmission probabilities, draws stochastic transmission events, and sets newly infected individuals to status `"e"`. Records the S→E flow as `se.flow`.

### Progression Module (`progress`)

Replaces EpiModel's built-in recovery module. Handles up to three transitions:

1. **E → I** at rate `ei.rate`: exposed individuals become infectious
2. **I → R** at rate `ir.rate`: infectious individuals recover
3. **R → S** at rate `rs.rate` (when `rs.rate > 0`): recovered individuals lose immunity (SEIRS extension)

When `rs.rate = 0`, the R→S block is skipped entirely, giving standard SEIR dynamics. This parameterized design means one module function handles both SEIR and SEIRS without code duplication.

## Parameters

### Transmission

| Parameter | Description | Default |
|-----------|-------------|---------|
| `inf.prob` | Per-act transmission probability | 0.5 |
| `act.rate` | Acts per partnership per timestep | 2 |

### Disease Progression

| Parameter | Description | Default |
|-----------|-------------|---------|
| `ei.rate` | E→I rate (1 / mean latent duration) | 1/50 |
| `ir.rate` | I→R rate (1 / mean infectious duration) | 1/75 |
| `rs.rate` | R→S waning immunity rate (0 = permanent immunity) | 0 (SEIR) or 1/100 (SEIRS) |

### Network

| Parameter | Description | Default |
|-----------|-------------|---------|
| Population size | Number of nodes | 500 |
| Target edges | Mean concurrent partnerships | 150 |
| Target isolates | Nodes with degree 0 | 240 |
| Partnership duration | Mean edge duration (timesteps) | 25 |

> **Note:** These parameters are stylized for pedagogical clarity, not calibrated to a specific pathogen. For influenza, typical values would be `ei.rate ~ 1/2` and `ir.rate ~ 1/5`.

## Module Execution Order

```
resim_nets → infection (infect) → progress → prevalence
```

The built-in `prevalence.net` module runs last and computes `s.num`, `i.num`, and `num` from the updated status attribute. The custom `progress` module additionally tracks `e.num`, `r.num`, and all transition flows.

## Next Steps

- **Add vital dynamics** (births and deaths) to maintain population size over long runs and allow true endemic equilibrium — see [SI with Vital Dynamics](../2018-08-SIwithVitalDynamics)
- **Incorporate vaccination** using all-or-nothing or leaky mechanisms — see [SEIR with AON Vaccination](../2018-10-SEIRwithAONVax) and [SEIRS with Leaky Vaccination](../2018-12-SEIRSwithLeakyVax)
- **Vary progression rates** by individual attributes (e.g., age, immune status) to model heterogeneous disease courses
- **Add stage-dependent infectiousness** where transmission probability depends on disease stage (e.g., higher during early infection) — see the [HIV model](../2019-03-HIV) for this pattern
- **Calibrate to a specific pathogen** by substituting published estimates for latent and infectious period durations

## Authors

Samuel M. Jenness, Emory University (http://samueljenness.org/)

Venkata R. Duvvuri
