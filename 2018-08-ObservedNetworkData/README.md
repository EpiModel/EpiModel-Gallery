# Modeling Epidemics over Observed Networks

## Description

This example demonstrates how to model epidemics over **observed (census) dynamic networks** in EpiModel. The standard EpiModel workflow starts with egocentrically observed network data, fits a temporal ERGM (exponential random graph model) with the generative network effects of interest, then simulates from that model fit over time. This example takes an alternative approach: using a dynamic network census where all nodes and edges are observed over a series of discrete time steps.

This is the only gallery example that bypasses EpiModel's ERGM estimation and simulation pipeline entirely. Instead, the observed `networkDynamic` object is placed directly into the simulation data structure, and `resimulate.network` is set to `FALSE` so the network is used as-is rather than being regenerated at each timestep. This requires custom initialization and infection modules.

The example uses a simulated dataset (`concurrencyComparisonNets`) from the `networkDynamicData` package, treated as if it were an observed network census. Two scenarios are demonstrated: a basic SI model with constant transmission probability, and a two-stage disease model with time-varying transmission probability and network visualization.

## Model Structure

### Disease Compartments

| Compartment | Label | Description |
|-------------|-------|-------------|
| Susceptible | **S** | Not infected; at risk of infection |
| Infectious | **I** | Infected; can transmit to susceptible partners |

### Flow Diagram

```mermaid
flowchart LR
    S["<b>S</b><br/>Susceptible"] -->|"infection<br/>(si.flow)"| I["<b>I</b><br/>Infectious"]

    style S fill:#3498db,color:#fff
    style I fill:#e74c3c,color:#fff
```

This is a simple SI model with no recovery. In a closed population (no vital dynamics), prevalence increases monotonically toward saturation.

### Observed Network Data

Unlike other gallery examples, this model does not call `netest()` to fit an ERGM. Instead, the observed `networkDynamic` object contains edge spells (onset/terminus times for each partnership) recorded over ~100 discrete time steps. Key characteristics:

- **1000 nodes**, undirected network
- **Edge spells** define when each partnership is active
- The existing `status.active` vertex attribute is removed before simulation, since we simulate our own epidemic
- **No network resimulation**: the observed edge structure is fixed throughout the simulation

## Modules

### Initialization Module (`init_obsnw`)

Replaces EpiModel's built-in `initialize.net`, which expects a `netest` object (fitted ERGM). Since observed networks have no model fit, this custom module handles initialization directly:

1. Creates the master data list via `create_dat_object()`
2. Places the observed `networkDynamic` object directly into `dat$run$nw[[1]]`
3. Initializes core node attributes and randomly assigns initial infections
4. Optionally stores time-varying disease status on the network (when `track.nw.attr = TRUE`), enabling network visualization with `plot(sim, type = "network", col.status = TRUE)`
5. Calls `prevalence.net()` to record initial epidemic statistics

### Infection Module (`infect_obsnw`)

A streamlined transmission module that works with the observed network. It supports two modes, selected automatically based on which parameters are provided:

- **Simple mode** (when `inf.prob` is set): constant per-act transmission probability for all infected individuals
- **Time-varying mode** (when `inf.prob.stage1` is set): transmission probability depends on infection duration. During the primary stage (`dur.stage1` timesteps), transmission occurs at `inf.prob.stage1`; afterward, at `inf.prob.stage2`.

In both modes, the module uses `discord_edgelist()` to identify susceptible-infectious partnerships at each timestep, computes per-partnership transmission probabilities adjusted for the act rate, and stochastically determines new infections.

## Parameters

### Example 1: Basic SI

| Parameter | Description | Value |
|-----------|-------------|-------|
| `inf.prob` | Per-act transmission probability | 0.5 |
| `act.rate` | Acts per partnership per timestep | 1 |

### Example 2: Time-Varying Transmission

| Parameter | Description | Value |
|-----------|-------------|-------|
| `inf.prob.stage1` | Transmission probability during primary stage | 0.05 |
| `inf.prob.stage2` | Transmission probability during secondary stage | 0.15 |
| `dur.stage1` | Duration of primary stage (timesteps) | 5 |
| `act.rate` | Acts per partnership per timestep | 1 |
| `track.nw.attr` | Track time-varying status on network for visualization | TRUE |

## Caveats: Observed Network Boundaries

The observed network has edge activity recorded over a finite window (~100 timesteps). Edges that are active at the last observed time remain active indefinitely (a `networkDynamic` convention). This means:

- Nothing prevents simulating past the observation window
- However, results become meaningless: the frozen edge set no longer reflects real dynamics, so prevalence plateaus and incidence drops to zero
- **Always match `nsteps` to the number of observed time steps in the network**

The interactive section of the model script demonstrates this by simulating 200 timesteps on a ~100-step network.

## Next Steps

- Model a different disease type (SIS, SIR) by adding recovery or immunity modules
- Add vital dynamics (births, deaths) to study longer-term epidemics -- see [SI with Vital Dynamics](../2018-08-SIwithVitalDynamics)
- Use a different dataset from the `networkDynamicData` package
- Take a static network dataset (e.g., from the `ergm` package) and fit a TERGM on it with an assumed edge dissolution rate, bridging the observed and model-based approaches

## Author

Samuel M. Jenness, Emory University (<http://samueljenness.org/>)
