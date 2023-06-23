# Epidemics with Multiple Networks

## Description

This example shows how to use core EpiModel functionality to include multiple,
interacting networks in a single epidemic model. The node set is required to be
the same for all networks. The edge formation and dissolution models may vary
from one network to the next, and the formation model for one network may depend
on the edge states in the other networks, represented via edge-dependent nodal
attributes. We give an example with two networks, including cross-network
dependency. The two network models are fit successively and then passed jointly
to `netsim`.

### Modules

We use built-in modules only, to illustrate their capabilities.

### Parameters

The cross-network dependency requires that we maintain the edge-dependent nodal
attributes during the `netsim` run via an appropriately specified `dat.updates`
argument to `control.net`.

Additionally, certain arguments to `control.net` may vary from one network to
the next via the `multilayer` specification, and we illustrate this
functionality here for the `nwstats.formula` argument.

## Next Steps

Multinetwork functionality is new to core EpiModel, and there are many possible
extensions. In the future it may be possible to have transmission probabilities
vary by edge type, for example.

## Author

Chad Klumb, University of Washington
(borrowing some code and comments from the SEIR example)
