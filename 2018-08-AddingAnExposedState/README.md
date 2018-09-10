# SEIR Model: Adding an Exposed State to an SIR

## Description
EpiModel includes an integrated SIR model, but here we show how to model an SEIR disease like influenza. The `E` compartment in this disease is an exposed state in which a person has been infected but is not infectious to others. Many infectious diseases have this latent non-infectious stage, and in general it provides a general framework for transmission risk that is dependent on one’s stage of disease progression.

### Modules
The built-in **infection module** (function = `infect`) is designed to handle a wide variety of specifications to the integrated network models. The main processes within the module are to extract a _discordant edgelist,_ which is a matrix of ID numbers of active dyads in the network in which one member of the dyad is disease-susceptible and the other is infected. Given epidemic parameters for the probability of infection per act and the number of acts per unit time, the per-partnership transmission rate is calculated as an exponential function of the two. The key component that required updating in this module was to specify the newly infected person’s infection status as `"e"`, whereas the default value transitions persons to `"i"`.

The **disease progression module** (function = `progress`) will simulate the transition from exposed to infectious, and also the transition from infectious to recovered in the model. Instead both progression events are individual-level stochastic processes, in which there is a constant hazard of transition. The times spent in each disease compartment therefore follow a geometric distribution. The code block here is for the transition between exposed to infectious states. For each eligible person for transition, the event is modeled as a stochastic process following a Bernoulli distribution with the rate parameter, `ei.rate`. Persons who do transition to the infectious state have their individual-level status attribute updated to the `"i"` value, at which point they are now capable of infecting others.

### Parameters
The listing of main epidemic model parameters is as follows: 

* `inf.prob` is the probability that an infection will occur given an act between a susceptible and infected node 
* `act.rate` is the number of acts per partnership per unit time 
* `ei.rate` is the rate of exposed persons moving to the infectious state (1/average duration spent in `E`) 
* `ir.rate` is the rate of infectious persons moving to the recovered state (1/average duration spent in `I`)


### Extension #1: Adding an R --> S Transition (SEIRS Model)
In an extension contributed by Venkata R. Duvvuri, the SEIR model was expanded into an SEIRS model by adding an additional transition from the recovered state back into the susceptible state. The new parameter added:

* `rs.rate` is the rate of recovered persons who loss their immunity after some time due to waning of immunity moving to the susceptible state (1/average duration spent in `R`)

## Next Steps
Good next steps for this example might be to incorporate a vaccination or intervention strategy, or to vary the rate of disease progression based on an additional attribute of persons in the network.

## Authors
Samuel M. Jenness, Emory University (http://samueljenness.org/)
Venkata R. Duvvuri
