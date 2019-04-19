# HIV Model

## Description
This example shows how to model a simple, one-mode HIV transmission process over a dynamic network with vital dynamic processes for population arrival and departure. 

EpiModel includes an integrated SI model, but here we show how to model a one-mode SI disease with four distinct infectious - I - sub-compartments, like HIV.   The `E` compartment in this disease is an exposed state in which a person has been infected but is not infectious to others. Many infectious diseases have this latent non-infectious stage, and in general SEIRS modeling provides a general framework for transmission risk that is dependent on one's stage of disease progression. For more background on the SEIRS model, see [SEIR Model: Adding an Exposed State to an SIR](https://github.com/statnet/EpiModel-Gallery/tree/master/2018-08-AddingAnExposedState "EpiModel Gallery - SEIR Model"). Note that while the gallery example was originally built to model SEIR, the model was extended to include SEIRS modelling capabilities.

The arrival module implements a stochastic entrance process that acts as a function of a standard arrival rate. The arrival module has been modified to additionally simulate a stochastic leaky vaccination and vaccine protection processes. In a leaky vaccine model, the number of susceptible individuals who are conferred some degree of vaccine protected is a product of the fraction vaccinated $\omega$ omega) and the fraction of the vaccinated who are conferred protection by the vaccine $\chi$ (chi). The degree of protection conferred to vaccine protected individuals is the vaccine efficacy $\psi$ (psi). Susceptible individuals who are conferred vaccine protection in a leaky vaccine model remain susceptible to infection from the disease but with a reduced force of infection - likelihood of becoming infected - compared with unvaccinated individuals. This reduced force of infection through vaccine protection may be represented by: (1 - vaccine efficacy) * (probability of infection). This contrasts with the all-or-nothing vaccine model where vaccine protected individuals were conferred complete immunity to the disease and could not become infected. 

An assumption of this leaky vaccine model is that individuals may not be vaccinated more than once. As a result, individuals who are vaccinated but are not conferred vaccine protection will not have the opportunity to become vaccine protected through subsequent vaccinations. An additional assumption is that natural disease experience does not result in reduced force of infection after the individual has progressed from recovered to susceptible. Finally, the model assumes that even if a vaccine protected individual becomes infected, once they progress through the compartments of the model and again become susceptible they will retain their original vaccine protection and have a reduced force of infection.

The departure module implements a stochastic exit process that acts as a function of either a standard departure rate or a disease-induced departure rate.

### Modules
The **departure module** (function = `dfunc`)  simulates departure as a function of a disease-induced departure rate. The module relies on the fact that a standard departure rate is passed in by the module user as a rate - `departure.rate` - in the epidemic model parameter settings. For those eligible (active) individuals whose disease status is `"i"`, the departure rate is multiplied by the value of the `departure.disease.mult` parameter representing an increased likelihood of model departure for infected indviduals.

The **arrival module** (function = `afunc`) has been extended from the birth model included in the github SEIR model example in order to include new logic around implementation of a leaky vaccine intervention.
Vaccination and protection is simulated through three methods:

1. Initialization - The starting population of individuals in the model may have received vaccination and may have some degree of vaccine protection as a result. Vaccination and vaccine protection rates can be user-defined using the vaccination.rate.initialization and protection.rate.initialization model parameters.
2. Progression - Unvaccinated individuals in the model have the opportunity to become vaccinated (e.g. through vaccination campaigns). However, only those individuals that receive vaccination when they are susceptible may confer some degree of vaccine protection. For each timestep in the model, vaccination of unvaccinated individuals is simulated followed by protection of those vaccinated individuals who are susceptible based on the `vaccination.rate.progression` and `protection.rate.progression` model parameters.
3. Arrival - Individuals arriving into the population have the opportunity to become vaccinated. The vaccination may confer some degree of vaccine protection. Rates of vaccination and vaccine protection in individuals born into the network is determined by the `vaccination.rate.arrivals` and `protection.rate.arrivals` parameters.

### Parameters
The epidemic model parameters include those needed for the SEIRS model and those pertaining to vaccination and vaccine protection rates.

* `inf.prob`: the probability that an infection will occur given an act between a susceptible and infected node. 
* `act.rate` the number of acts per partnership per unit time 
* `ei.rate` the rate of exposed persons moving to the infectious state (1/average duration spent in `E`) 
* `ir.rate` the rate of infectious persons moving to the recovered state (1/average duration spent in `I`)
* `mortality.rate`: a scalar for the standard mortality rate of the population. For disease status of `"i"`, this rate is multiplied by the `mortality.disease.mult` explained below.
* `mortality.disease.mult`: the multiplier acting on mortality rate for persons with a disease status of `"i"`. 
* `birth.rate`: a scalar for the rate of births per person per week.
* `vaccination.rate.initialization`: A scalar for the proportion of the population at time 0 who are vaccinated.
* `protection.rate.initialization`: A scalar for the proportion of the population at time 0 who are confered vaccine protection after becoming vaccinated.
* `vaccination.rate.progression`: A scalar for the proportion of the unvaccinated population who become vaccinated throughout the simulation progression.
* `protection.rate.progression`: A scalar for the proportion of the susceptible and newly vaccinated population who become vaccinated through the vaccination progression method and who are confered vaccine protection.
* `vaccination.rate.arrivals`: A scalar for the proportion of the individuals who enter into the network vaccinated.
* `protection.rate.arrivals`: A scalar for the proportion of individuals who enter into the network vaccinated and are confered vaccine protection.
* `vaccine.efficacy`: A scalar for the degree of protection resulting from vaccination and vaccine protection, resulting in a reduced force of infection.

## Authors
Connor M. Van Meter, Emory University
