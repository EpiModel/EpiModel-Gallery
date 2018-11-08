# Modeling Information Diffusion Process in a Social Network Using SI Model

## Description
In this example we build an information diffusion model based on exsiting SI process in EpiModel to monitor dissemination of new ideas inside a social network. Instead of information diffusion (infection process) happening between any discordant nodes, we monitor two scenarios in which either it can only happen when suscpetible individuals have more than a minimum degree of partnerships with infectious people, or the probability of idea transmission is the logistic function of number of susceptible individuals'partnership with infectious individuals. New argument includes minimum degree of infection, log odds of probability with 1 degree increase of partnership, degree at 50% transmission probability. No entirely new modules are needed, but infection process module is edited: 

## Modules
The **infection module** (function = `infect_mod`) includes the following changes from the base EpiModel infection module (`infection.net`):

* Book-keeping the degree of discordant edges for susceptible nodes

* The infection probability is only assigned to susceptible nodes with more than minimum degree of discordant edges, otherwise is 0

The **infection module** (function = `infect_newmod`) includes the following changes from the base EpiModel infection module (`infection.net`):

* Book-keeping the degree of discordant edges for susceptible nodes

* The infection probability is a logistic function of degree of discordant edges for susceptible nodes with user specified paramters
 
### Parameters

The new or altered epidemic model parameters are:

* `min.degree`: the minimum degree of discordant edges for the infection to occur 

* `log_a`:

* `log_b`:

# Worked example

In the worked example in `model.R`, we consider a case in which strain 1 is highly infectious but short-lived (`inf.prob` = 0.5, `rec.rate` = 0.05), while strain 2 is much less infectious but longer-lived (`inf.prob.st2` = 0.01, `rec.rate.st2` = 0.005). The initial frequency of the two strains is equal (`pct.st2` = 0.5).

We first compare two scenarios: one in which all relationships are equally likely, and one in which all individuals are limited to one partner at a time. These two dynamic network models possess the same expected mean degree, relational durations, act rate, and all other parameters.

Our simulations demonstrate that strain 1 dominates strain 2 in the completely random network, but strain 2 dominates in the strict-monogamy network. Indeed, strain 1 goes extinct in the latter case. We leave it to the reader to deduce for themselves the logic behind this pattern.

We then consider intermediate levels of relational concurrency to understand how and where this transition from one dominant strain to the other occurs. We find that Strain 1's absolute prevalence increases monotonically with greater concurrency (with some stochasticity), while strain 2 increases and then decreases again, although with less variation overall. The two cross over at around 70 individuals (out of 1000) practicing concurrency on average at any point in time.

## Author

Samuel Jenness, Emory University  

Yuan Zhao, Emory University   

