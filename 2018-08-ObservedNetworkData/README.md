# Modeling Epidemics over Observed Networks

## Description

This example shows how to model epidemics over observed dynamic networks in EpiModel. The standard approach for modeling epidemics over networks in EpiModel is to start with egocentrically observed network data, fit a temporal ERGM with the generative network effects of interest, then simulate from that model fit over time. An alternative approach is possible when one has a dynamic network census: an observed network with all nodes and edges observed over a series of discrete time steps. One needs to be quite careful to understand the structure and potential limitations of this type of data, since the opportunities for missingness are present; but we will ignore those issues assume we have perfect network data.

### Modules

The **initialization module** (function = `new_init_mod`) handles the process of setting up the epidemic simulation. The default module contained in the function `initialize.net` serves as our starting point for writing our own new module, `net.init.mod`. The default module does lots of work, such as simulating an initial network from the model fit, that we do not need to do with an observed network. Here, there are three key steps for the initialization module: set up the master data list, `dat`, including its core data structures; initialize infection among the nodes in the network; and use the `get_prev.net` function to record summary epidemiological statistics. With this modules, whereas `x` would ordinarily be the fitted network model from `netest`, now it will just be the networkDynamic object detailed above.

The **infection module** (function = `new_infect_mod`) handles the process of disease transmission over the observed network. The infection module must also be updated from the built-in version because of some features within the default function that depend on having a fitted network model. So the default module function, `infection.net`, serves as the basis of our new module, `my.inf.mod`. It is a stripped down version of the default that provides a much clearer picture of the processes within the module, but it is not general enough to handle all the epidemic modeling cases supported within EpiModel (e.g., time-vary infection probabilities or simulating epidemics over bipartite networks). The key element within both the default and this updated module is the `discord_edgelist` function that examines the current state of the network at `at`, and returns a matrix of disease discordant pairs (dyads over which disease may be transmitted).

### Parameters

The epidemic model parameters are basic here because we're not changing any of the core epidemiology from a simple SI model.

-   `inf.prob`: the probability that an infection will occur given an act between a susceptible and infected node

### Example 2: Adding Networking Tracking and Time-Varying Risk

As an update to the first example above, we updated the modules to track individual-level disease status over time and implement time-varying infection risk. Tracking changing disease status allows us to plot the network at various time points in the simulation, similar to the built-in epidemic models in EpiModel. We removed this functionality for the primary example here but it is easy to add back in; it requires using the `activate.vertex.attribute.active` function from the `networkDynamic` package. Second, we implemented a time-varying transmission probability that is a function of the duration of infection of the infected vertex in a disease-discordant dyad. In this example, the primary stage of disease lasts 5 time steps and has a corresponding transmission probability of 5% per act, while the secondary disease stage lasts for the remainder of the disease duration (infinite in this closed-population SI model) and has a corresponding transmission probability of 15%.

## Next Steps

Next steps for this example might be to model a different disease type or add other processes such as vital dynamics, or use a different dataset in the `networkDynamicData` package, or take a static network dataset (e.g., from the `ergm` package) and fit a TERGM on it with an assumed edge dissolution rate.

## Author

Samuel M. Jenness, Emory University (<http://samueljenness.org/>)
