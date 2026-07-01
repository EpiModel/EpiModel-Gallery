# Multilayer Networks: Independent Layers

## Description

This example is the entry point to **multilayer network** modeling in EpiModel.
A multilayer network consists of two or more network layers that share the same
node set but have different edge sets, representing distinct types of
relationships. For example, the two layers might represent main partnerships
(steady, long-lasting) and casual partnerships (more transient). An epidemic can
spread over edges in either layer.

The layers here are **independent**: a person's number of partners in one layer
does not constrain the number in the other. That keeps the mechanics minimal, so
the only genuinely new idea is that you fit each layer as an ordinary
single-layer network model and then hand `netsim()` a **list** of them. There is
no simulated annealing, no cross-layer nodal attribute, and no update callback.
Those are required only when the layers depend on each other, which is the
companion [Multilayer Networks: Cross-Layer Dependency](../multinets/) example.

If you have built a single-layer `edges` model with `netest()` and `netsim()`,
you already know almost everything here. This example uses a simple SI model in a
closed population to isolate the multilayer network mechanics from disease model
complexity. No custom `module-fx.R` is needed; it uses EpiModel's built-in SI
modules.


## Model Structure

### Network Layers

The key picture to hold in mind is one set of people with two separate edge sets
drawn over them. A node can have several main partners and no casual partners, or
the reverse. Because the layers are independent, those two partner counts are
unrelated.

| Property | Layer 1 (main) | Layer 2 (casual) |
|----------|----------------|-------------------|
| **Interpretation** | Steady, long-lasting partnerships | Transient, higher-turnover contacts |
| **Formation terms** | `edges + nodematch("race")` | `edges + degree(1)` |
| **Target edges** | 90 | 75 |
| **Mean degree** | 0.36 | 0.30 |
| **Mean duration** | 200 time steps | 20 time steps |
| **Layer-specific terms** | Race homophily (`nodematch`) | Degree-1 count (`degree(1)`) |
| **Cross-layer term** | none (independent) | none (independent) |

### What Independence Buys You

Each layer's formation model references only its own structure. Layer 2's
formula (`edges + degree(1)`) never mentions layer 1, so the two ERGMs can be fit
and resimulated in complete isolation. This is what removes the need for the
cross-layer machinery used in the companion example:

| Step in the cross-layer example | Needed here? |
|---------------------------------|--------------|
| `san()` bootstrap of starting degree attributes | No |
| Cross-layer `nodefactor()` formation terms | No |
| `dat.updates` callback between layer resimulations | No |
| Pass a list of `netest` objects to `netsim()` | Yes |


## Modules

This example uses **built-in EpiModel SI modules only**. No custom `module-fx.R`
is required. The multilayer functionality lives entirely in the network
estimation and simulation layers, not in the epidemic modules.


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
| Race homophily (`nodematch`) | 60 | none |
| Degree-1 count (`degree(1)`) | none | 120 |

### Network Dissolution

| Parameter | Layer 1 | Layer 2 |
|-----------|---------|---------|
| Mean duration (time steps) | 200 | 20 |


## Key EpiModel Functions for Multilayer Models

| Function / Argument | Purpose |
|----------------------|---------|
| `netsim(list(est1, est2), ...)` | Passing a list of `netest` objects signals a multilayer model; list order sets the layer numbering |
| `resimulate.network = TRUE` | Redraw the dynamic layers each time step |
| `multilayer()` in `nwstats.formula` | Optional per-layer diagnostic statistics |
| `print(sim, network = k)` / `plot(sim, network = k, ...)` | Inspect a specific layer |


## Analysis

The script produces:
1. **Network diagnostics** for each layer (formation statistics tracking targets)
2. **Per-layer formation diagnostics** during the epidemic simulation
3. **SI compartment plot** showing susceptible and infected counts over time
4. **Prevalence curve** over time
5. **Incidence curve** (new infections per time step)
6. **Summary table** with final prevalence, cumulative infections, and peak incidence

The short-duration casual layer tends to drive early spread, while the
long-duration main layer sustains transmission over time. Because transmission
can occur on either layer, the overall epidemic reflects contributions from both.


## Next Steps

- **Cross-layer dependency**: make a person's activity in one layer reduce their
  activity in the other (finite relational capacity, negative degree
  correlation). This adds a cross-layer formation term, a one-time `san()`
  bootstrap, and an update callback. See
  [Multilayer Networks: Cross-Layer Dependency](../multinets/).
- **Layer-specific transmission**: different `inf.prob` or `act.rate` by edge
  type (e.g., lower transmission on casual ties).
- **More than two layers**: add a third layer (e.g., a one-off contact layer with
  `duration = 1`) and pass `netsim(list(est1, est2, est3), ...)`.
- **Different disease models**: extend to SIR or SEIR with a custom progression
  module, see [Adding an Exposed State](../seir-exposed-state/).


## Authors

Samuel Jenness (Emory University)
