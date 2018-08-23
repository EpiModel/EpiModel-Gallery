# SI Model with Age Mixing and Vital Dynamics

## Description
This example shows how to model a relatively simple SI epidemic over a dynamic network, but with vital dynamic processes for aging, births, and deaths. The modules implement a different approach to mortality, in which the stochastic death process is a function of both an age-specific mortality rate and disease-induced mortality. The network model uses an `absdiff` term to specify a parametric form of age mixing, and that age attribute that is initialized on the network is then pulled into the epidemic modules. 

### Modules
The **aging module** (function = `aging`) sets up the age attribute at the initial time step by pulling from the network object. Then it will subsequently update the age attribute as updating the age in increments of a week at each time step.

The **death module** (function = `dfunc`)  simulates mortality as a function of an age-specific mortality rate and disease-induced mortality. The module relies on the fact that age specific mortality rates are passed in as a vector of rates in which the position of in the vector corresponds to `age+1`. For those eligible (alive) individuals whose disease status is `"i"`, the subsequent mortality rates are multiplied by the value of the `mortality.disease.mult` parameter.

The **birth module** (function = `bfunc`) implements a more simplified birth process compared to the built-in birth module, mainly for greater legibility. Note that the `age` attribute must be set on both the `attr` list and the `nw` object for new incoming nodes. 

### Parameters
The epidemic model parameters are basic here because we're not changing any of the core epidemiology from a simple SI model.

* `inf.prob`: the probability that an infection will occur given an act between a susceptible and infected node. 
* `mortality.rates`: a vector of mortality rates, where the position on the vector correspond to the mortality rate of persons of `age+1` (for example, the first rate is for persons less than one year of age). 
* `mortality.disease.mult`: the multiplier in age-specific mortality rates for persons with a disease status of `"i"`. 
* `birth.rate`: a scalar for the rate of births per person per week. 

## Author
Samuel M. Jenness, Emory University (http://samueljenness.org/)
