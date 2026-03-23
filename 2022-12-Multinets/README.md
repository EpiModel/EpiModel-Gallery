# Epidemics with Multiple (Multilayer) Networks

## Description

This example demonstrates how to simulate an epidemic on a **multilayer
network** in EpiModel. A multilayer network consists of two or more network
layers that share the same node set but have different edge sets — representing
distinct types of relationships. For example, the two layers might represent
main partnerships and casual contacts, or sexual contacts and needle-sharing
contacts.

The key feature of multilayer models is **cross-layer dependency**: the degree
(number of edges) in one layer can influence edge formation in the other. This
models the realistic constraint that individuals have finite relational capacity
— being active in one type of partnership reduces availability for the other.

This example uses a simple SI model in a closed population to isolate the
multilayer network mechanics from disease model complexity. No custom
`module-fx.R` is needed — it uses EpiModel's built-in SI modules.


## Model Structure

### Network Layers

| Property | Layer 1 | Layer 2 |
|----------|---------|---------|
| **Interpretation** | Longer-duration partnerships | Shorter-duration partnerships |
| **Formation terms** | `edges + nodematch("race") + nodefactor("deg1+.net2")` | `edges + degree(1) + nodefactor("deg1+.net1")` |
| **Target edges** | 90 | 75 |
| **Mean degree** | 0.36 | 0.30 |
| **Mean duration** | 100 time steps | 75 time steps |
| **Layer-specific terms** | Race homophily (`nodematch`) | Degree-1 count (`degree(1)`) |
| **Cross-layer term** | `nodefactor("deg1+.net2")` | `nodefactor("deg1+.net1")` |

### Cross-Layer Dependency

The `nodefactor("deg1+.netX")` terms create negative degree correlation across
layers. The target statistic for these terms (10 out of 90 or 75 total edges)
is deliberately small relative to total edges, producing strong negative
correlation: most nodes that are active in one layer are inactive in the other.

The cross-layer dependency is maintained during simulation by the
`network_layer_updates()` callback function, passed via the `dat.updates`
argument to `control.net()`. This function runs between layer resimulations to
update the degree attributes that each layer's formation model depends on:

| `network` value | Timing | Action |
|-----------------|--------|--------|
| 0 | Before layer 1 resimulation | Update `deg1+.net2` from current layer 2 edges |
| 1 | Before layer 2 resimulation | Update `deg1+.net1` from current layer 1 edges |
| 2 | After layer 2 resimulation | Recalculate `deg1+.net2` for summary statistics |

### Network Setup with `san()`

Because the two layers depend on each other's degree distribution, fitting the
ERGM models requires reasonable starting values for the cross-layer degree
attributes. The `san()` function (simulated annealing from the `ergm` package)
generates edge sets consistent with target statistics before estimation. This
is a two-step process:

1. Generate layer 1 edges consistent with target mean degree and race homophily
2. Generate layer 2 edges accounting for the cross-layer constraint from step 1


## Modules

This example uses **built-in EpiModel SI modules only**. No custom
`module-fx.R` is required. The multilayer functionality is entirely in the
network estimation and simulation layers, not in the epidemic modules.


## Parameters

### Epidemic

| Parameter | Description | Value |
|-----------|-------------|-------|
| `inf.prob` | Per-act transmission probability | 0.5 |
| `act.rate` | Acts per partnership per time step | 2 |

Both are deliberately high to generate a clearly visible epidemic, since the
focus of this example is network mechanics rather than epidemiological realism.

### Network Formation

| Parameter | Layer 1 | Layer 2 |
|-----------|---------|---------|
| Target edges | 90 | 75 |
| Race homophily (`nodematch`) | 60 | — |
| Degree-1 count (`degree(1)`) | — | 120 |
| Cross-layer `nodefactor` | 10 | 10 |

### Network Dissolution

| Parameter | Layer 1 | Layer 2 |
|-----------|---------|---------|
| Mean duration (time steps) | 100 | 75 |


## Key EpiModel Functions for Multilayer Models

| Function / Argument | Purpose |
|----------------------|---------|
| `netsim(list(est.1, est.2), ...)` | Passing a list of `netest` objects signals a multilayer model |
| `dat.updates` in `control.net()` | Callback function for cross-layer attribute updates |
| `multilayer()` in `nwstats.formula` | Specifies per-layer network statistics to track |
| `tergmLite = TRUE` | Uses lightweight edgelist representation (required for `dat$run$el[[]]` access in the update callback) |
| `resimulate.network = TRUE` | Required for multilayer models with cross-layer dependency |


## Analysis

The script produces:
1. **Network diagnostics** for each layer (formation statistics tracking targets)
2. **SI compartment plot** showing susceptible and infected counts over time
3. **Prevalence curve** over time
4. **Incidence curve** (new infections per time step)
5. **Summary table** with final prevalence, cumulative infections, and peak incidence


## Next Steps

Multilayer network functionality in EpiModel opens many possible extensions:

- **Layer-specific transmission**: different `inf.prob` or `act.rate` by edge
  type (e.g., lower transmission in casual vs. main partnerships)
- **More than two layers**: add a third network layer (e.g., needle-sharing)
- **Vital dynamics**: add arrivals and departures with custom modules
- **Different disease models**: extend to SIR, SEIR, or disease-specific models
- **Heterogeneous populations**: add more nodal attributes beyond the binary
  `race` attribute used here

For a more detailed tutorial on multilayer networks, see the
[SISMID EpiModel course materials](https://epimodel.github.io/sismid/).


## Authors

Chad Klumb (University of Washington) and Samuel Jenness (Emory University)
