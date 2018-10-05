# SEIR Model with All or Nothing Vaccination

## Description
This example shows how to model an all or nothing vaccination intervenion on a SEIR epidemic over a dynamic network with vital dynamic processes for population entrance (ex. birth) and exit (ex. death). 

EpiModel includes an integrated SIR model, but here we show how to model an SEIR disease like influenza. The `E` compartment in this disease is an exposed state in which a person has been infected but is not infectious to others. Many infectious diseases have this latent non-infectious stage, and in general it provides a general framework for transmission risk that is dependent on one’s stage of disease progression. 

The birth module implements a stochastic entrance process that acts as a function of a standard birth rate. The death module implements a stochastic exit process that acts as a function of either a standard mortality rate or a disease-induced mortality rate.

An all or nothing vaccine is implemented with parameters ω (omega), the fraction of both the population at time 0 and the population entering into the model that is vaccinated, and χ (chi), the fraction of both the population at time 0 and the population entering into the model that is vaccinated and that is also protected from disease infection, but not disease exposure, by the vaccine. The implication of this is that if an individual is vaccinated and is also protected from the vaccine, then they might be exposed to a disease but will never become infectious and able to transmit the disease to others.

### Modules
The **death module** (function = `dfunc`)  simulates mortality as a function of an age-specific mortality rate and disease-induced mortality. The module relies on the fact that a standard mortality rate is passed in by the module user as a rate - `mortality.rate` - in the epidemic model parameter settings. For those eligible (alive) individuals whose disease status is `"i"`, the subsequent mortality rates are multiplied by the value of the `mortality.disease.mult` parameter.

The **birth module** (function = `bfunc`) has been extended from the birth model included in the github SEIR model example in order to include new logic around implementation of an all or nothing vaccine intervention. 

### Parameters
The epidemic model parameters include those needed for the SEIR model and those .

* `inf.prob`: the probability that an infection will occur given an act between a susceptible and infected node. 
* `act.rate` the number of acts per partnership per unit time 
* `ei.rate` the rate of exposed persons moving to the infectious state (1/average duration spent in `E`) 
* `ir.rate` the rate of infectious persons moving to the recovered state (1/average duration spent in `I`)
* `mortality.rate`: a scalar for the standard mortality rate of the population. For disease status of `"i"`, this rate is multiplied by the `mortality.disease.mult` explained below.
* `mortality.disease.mult`: the multiplier acting on mortality rate for persons with a disease status of `"i"`. 
* `birth.rate`: a scalar for the rate of births per person per week.
* `vaccine.rate`: equivalent to ω (omega) discussed in the module description above. A scalar for the proportion of both the population at time 0 and the population entering into the model that is vaccinated
* `protection.rate`: equivalent to χ (chi) discussed in the module description above. A scalar for the proportion of both the population at time 0 and the population entering into the model that is both vaccinated and also protected from disease infection, but not disease exposure, by the vaccine

## Authors
Samuel M. Jenness, Emory University (http://samueljenness.org/)
Venkata R. Duvvuri
Connor M. Van Meter, Emory University