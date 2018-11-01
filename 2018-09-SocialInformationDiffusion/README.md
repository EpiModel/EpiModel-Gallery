# Modeling Two Competing Strains in an SIS Epidemic

## Description
In this example we build an SIS model that considers two strains of a pathogen that differ in their infectiousness and their recovery rates. New arguments include the infection probability for the second strain, the recovery rate for the second strain, and the proportion of the initial infections that entail the second strain. No entirely new modules are needed, but two of the built-in modules (infection and recovery) are edited: 

## Modules
The **infection module** (function = `infection.2strains`) includes the following changes from the base EpiModel infection module (`infection.net`):

* Existing functionality to pass different transmission probabilities by direction across two modes is turned off, and warnings added if users invoke this functionality

* Transmission probabilities within each discordant relationship are calculated based on the strain of the infected partner (and the time since infection, if time-varying parameters are used)

* Newly infected individuals have their strain assigned from the infected partner

* Book-keeping is added for incidence and prevalence by strain
 
The **recovery module** (function = `recov.2strains`) is based off of the recovery module in the TestAndTreatIntevention example in the EpiModel Gallery, which is in turn based off of the built-in recovery module (`recovery.net`). Three key changes are found here:

* Strain is assigned to the initially-infected population

* Recovery rates are calculated based on an infected individual's strain

* Book-keeping is added for recoveries by strain


### Parameters
The new or altered epidemic model parameters are:

* `inf.prob`: the probability that an infection will occur given an act between a susceptible individual and one who is infected with *strain 1* 
* `inf.prob.st2`: the probability that an infection will occur given an act between a susceptible individual and one who is infected with *strain 2*
* `rec.rate`: the rate of recovery for those infected with *strain 1*
* `rec.rate.st2`: the rate of recovery for those infected with *strain 2* 
* `pct.st2`: the probability that an initially infected individual is infected with *strain 2*

# Worked example
In the worked example in `model.R`, we consider a case in which strain 1 is highly infectious but short-lived (`inf.prob` = 0.5, `rec.rate` = 0.05), while strain 2 is much less infectious but longer-lived (`inf.prob.st2` = 0.01, `rec.rate.st2` = 0.005). The initial frequency of the two strains is equal (`pct.st2` = 0.5).

We first compare two scenarios: one in which all relationships are equally likely, and one in which all individuals are limited to one partner at a time. These two dynamic network models possess the same expected mean degree, relational durations, act rate, and all other parameters.

Our simulations demonstrate that strain 1 dominates strain 2 in the completely random network, but strain 2 dominates in the strict-monogamy network. Indeed, strain 1 goes extinct in the latter case. We leave it to the reader to deduce for themselves the logic behind this pattern.

We then consider intermediate levels of relational concurrency to understand how and where this transition from one dominant strain to the other occurs. We find that Strain 1's absolute prevalence increases monotonically with greater concurrency (with some stochasticity), while strain 2 increases and then decreases again, although with less variation overall. The two cross over at around 70 individuals (out of 1000) practicing concurrency on average at any point in time.

## Author
Steven M. Goodreau, University of Washington (http://faculty.washington.edu/goodreau)
