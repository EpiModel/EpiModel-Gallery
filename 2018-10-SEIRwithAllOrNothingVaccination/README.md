# SEIR Model with All or Nothing Vaccination

## Description
This example shows how to model an all or nothing vaccination intervenion on a SEIR epidemic over a dynamic network with vital dynamic processes for population entrance (e.g. birth) and exit (e.g. death). 

EpiModel includes an integrated SIR model, but here we show how to model an SEIR disease like influenza. The `E` compartment in this disease is an exposed state in which a person has been infected but is not infectious to others. Many infectious diseases have this latent non-infectious stage, and in general SEIR modeling provides a general framework for transmission risk that is dependent on one’s stage of disease progression. 

The birth module implements a stochastic entrance process that acts as a function of a standard birth rate. The birth module has been modified to additionally simulate vaccination and vaccine protection processes. 
The death module implements a stochastic exit process that acts as a function of either a standard mortality rate or a disease-induced mortality rate.

### Modules
The **death module** (function = `dfunc`)  simulates mortality as a function of an age-specific mortality rate and disease-induced mortality. The module relies on the fact that a standard mortality rate is passed in by the module user as a rate - `mortality.rate` - in the epidemic model parameter settings. For those eligible (alive) individuals whose disease status is `"i"`, the subsequent mortality rates are multiplied by the value of the `mortality.disease.mult` parameter.

The **birth module** (function = `bfunc`) has been extended from the birth model included in the github SEIR model example in order to include new logic around implementation of an all or nothing vaccine intervention.
Vaccination and protection is simulated through three methods:

1. Initialization - The starting population of individuals in the model may have received vaccination and may have vaccine protection as a result. Vaccination and vaccine protection rates can be manipulated using the vaccination.rate.initialization and protection.rate.initialization model parameters.
2. Progression - Unvaccinated individuals in the model have the opportunity to become vaccinated (e.g. through vaccination campaigns). However, only those individuals that receive vaccination when they are susceptible are considered vaccine protected. For each timestep in the model, vaccination of unvaccinated individuals is simulated followed by protection of those vaccinated individuals who are susceptible based on the `vaccination.rate.progression` and `protection.rate.progression` model parameters.
3. Birth - Individuals "born" into the population have the opportunity to become vaccinated (e.g. simulating the practice of vaccinating newborns). Vaccine protection may be conferred as a result of the vaccination. Rates of vaccination and vaccine protection in individuals born into the network is determined by the `vaccination.rate.births` and `protection.rate.births` parameters.

A key assumption of this model is that individuals may not be vaccinated more than once. As a result, individuals that are vaccinated but fail to confer vaccine protection will not have the opportunity to become vaccine protected through subsequent vaccinations.

### Parameters
The epidemic model parameters include those needed for the SEIR model and those pertaining to vaccination and vaccine protection rates.

* `inf.prob`: the probability that an infection will occur given an act between a susceptible and infected node. 
* `act.rate` the number of acts per partnership per unit time 
* `ei.rate` the rate of exposed persons moving to the infectious state (1/average duration spent in `E`) 
* `ir.rate` the rate of infectious persons moving to the recovered state (1/average duration spent in `I`)
* `mortality.rate`: a scalar for the standard mortality rate of the population. For disease status of `"i"`, this rate is multiplied by the `mortality.disease.mult` explained below.
* `mortality.disease.mult`: the multiplier acting on mortality rate for persons with a disease status of `"i"`. 
* `birth.rate`: a scalar for the rate of births per person per week.
* `vaccination.rate`: equivalent to ω (omega) discussed in the module description above. A scalar for the proportion of the population at time 0, the population entering into the model, and the population of unvaccinated individuals that may become vaccinated.
* `protection.rate`: equivalent to χ (chi) discussed in the module description above. A scalar for the proportion of the population at time 0, the population entering into the model, and the population of newly vaccinated, susceptible individuals that confer vaccine protection upon vaccination.
* `vaccination.rate.initialization`: A scalar for the proportion of the population at time 0 who are vaccinated.
* `protection.rate.initialization`: A scalar for the proportion of the population at time 0 who are confered vaccine protection after becoming vaccinated.
* `vaccination.rate.progression`: A scalar for the proportion of the unvaccinated population who become vaccinated throughout the simulation progression.
* `protection.rate.progression`: A scalar for the proportion of the susceptible and newly vaccinated population who become vaccinated through the vaccination progression method and who are confered vaccine protection.
* `vaccination.rate.births`: A scalar for the proportion of the individuals who enter into the network vaccinated.
* `protection.rate.births`: A scalar for the proportion of individuals who enter into the network vaccinated and are confered vaccine protection.

## Authors
Samuel M. Jenness, Emory University (http://samueljenness.org/)
Venkata R. Duvvuri
Connor M. Van Meter, Emory University
