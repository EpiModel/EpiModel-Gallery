# **SEIRS Model: Adding a transition from Recovered State to Susceptible state to SEIR model** 

#**Description**
EpiModel includes SI, SIS, SIR and SEIR models. The addition of transition from recovered (R) compartment to susceptible (S) compartment explains that the recovered population will become susceptible as they loss their immunity that is gained from a natural infection (like influenza) or a vaccination (like influenza) after some period due to waning of immunity.

# **Modules**
The built-in **infection module** (function = infect) is designed to handle a wide variety of specifications to the integrated network models. The main processes within the module are to extract a discordant edgelist, which is a matrix of ID numbers of active dyads in the network in which one member of the dyad is disease-susceptible and the other is infected. Given epidemic parameters for the probability of infection per act and the number of acts per unit time, the per-partnership transmission rate is calculated as an exponential function of the two. The key component that required updating in this module was to specify the rate of transition from recovered compartment to susceptible compartment after spending some time in recovered compartment. 

The disease **progression module** (function = progress) will simulate the transition from exposed to infectious, the transition from infectious to recovered and also transition from recovered to susceptible in the model. These progression events are individual-level stochastic processes, in which there is a constant hazard of transition. The times spent in each disease compartment therefore follow a geometric distribution. The code block here is for the transition between recovered to susceptible states. For each eligible person for transition, the event is modeled as a stochastic process following a Bernoulli distribution with the rate parameter, rs.rate. Persons who do transition to the susceptible state have their individual-level status attribute updated to the "r" value, at which point they are now losing their immunity gained from the natural infection or from the vaccination of susceptible for infection.

# **Parameters**
The listing of main epidemic model parameters is as follows:
•	inf.prob is the probability that an infection will occur given an act between a susceptible and infected node
•	act.rate is the number of acts per partnership per unit time
•	ei.rate is the rate of exposed persons moving to the infectious state (1/average duration spent in E)
•	ir.rate is the rate of infectious persons moving to the recovered state (1/average duration spent in I)
•	rs.rate is the rate of recovered persons who loss their immunity after some time due to waning of immunity moving to the susceptible state (1/average duration spent in R) 

# **Author**
Venkata R. Duvvuri (venkata.r.duvvuri@gmail.com)
