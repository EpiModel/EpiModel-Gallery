# Syphilis Progressing Model 

## Description
Here we show how to model the transmission and progression of syphilis among men who have sex with men using EpiModel SI model by adding additional compartments to model the multiple stages of syphilis progression which include infectious incubation, primary, secondary stages, a low infectious early latent stage and a non-infectious late latent stage. In addition, we also include diagnosis of symptomatic individuals and screening of general populations.

### Modules
The built-in **infection module** (function = `infect`) is designed to handle a wide variety of specifications to the integrated network models. The main processes within the module are to extract a _discordant edgelist,_ which is a matrix of ID numbers of active dyads in the network in which one member of the dyad is disease-susceptible and the other is infected. Given epidemic parameters for the probability of infection per act and the number of acts per unit time, the per-partnership transmission rate is calculated as an exponential function of the two. The key component that required updating in this module was to specify the different transmission probability based on the disease stages.

The **disease progression module** (function = `progress`) will simulate the transitions between different syphilis stages: incubating (newly infected), primary, secondary, early latent, late latent, and tertiary. All progression events are individual-level stochastic processes, in which there is a constant hazard of transition. The times spent in each disease compartment therefore follow a geometric distribution. For each eligible person for transition, the event is modeled as a stochastic process following a Bernoulli distribution with the rate parameter, `ipr.rate`, `prse.rate`, `seel.rate`,`elll.rate`, `llter.rate`. Persons who do transition to the infectious state have their individual-level status attribute updated to the `"i"` value and syphilis stages updated accordingly, at which point they are capable of infecting others with corresponding transmission probabilities.

The **Diagnosis, screening and treatment module** (function = `tnt`) will simulate the treatment of symptomatic patients and screening of asymptomatic patients at different syphilis stages.

### Parameters
The listing of main epidemic model parameters is as follows: 

* `inf.prob1` is the probability that an infection will occur given an act between a susceptible and infected node at incubating, primary and secondary syphilis stages (0.18) 
* `inf.prob2` is the probability that an infection will occur given an act between a susceptible and infected node at early latent (0.09)
* `act.rate` is the number of acts per partnership per unit time 
* `ipr.rate` is the rate of incubating stage moving to the primary stage (1/average duration spent in `incubating`, which is around 28 days) 
* `prse.rate` is the rate of primary stage moving to the secondary stage (1/average duration spent in `primary`, which is around 63 days)
* `seel.rate` is the rate of secondary stage moving to the early latent stage (1/average duration spent in `secondary`, which is around 119 days)
* `elll.rate` is the rate of early latent stage moving to the late latent stage (1/average duration spent in `early latent`, which is around 154 days)
* `llter.rate` is the rate of late latent stage moving to the tertiary stage (1/average duration spent in `late latent`, which is around 29 years)
* `pri.sym` is the probability of showing symptomatic for each week during primary stage (0.205)
* `sec.sym` is the probability of showing symptomatic for each week during secondary stage (0.106)
* `early.trt` is the probability of receiving treatment given symptoms in primary and secondary stages (0.8)
* `late.trt` is the probability of receiving treatment given symptoms in the tertiary 
state (1.0)
* `scr.rate` is the probability of get screened as general population (yearly screening on average 1/52)

## Next Steps
Good next steps for this example might be to vary the rate of disease progression based on an additional attribute of persons in the network.

## Authors

Samuel M. Jenness, Emory University (http://samueljenness.org/)

Yuan Zhao
