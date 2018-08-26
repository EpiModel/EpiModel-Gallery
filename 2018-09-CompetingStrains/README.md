# Modeling Two Competing Strains in an SIS Epidemic

## Description
This example demonstrates how to XXX. Basic idea. New arguments: infection probabilities for each strain, mortality rate for each strain, recovery rate for each strain.

### Modules
The **infection module** (function = `infection.2strains`) has two changes from the base Epimodel infection module.  The first is that the transmission probability depends on the strain of the infected partner. The second is that newly infected individuals wiull have their strain assigned from the infected partner.

The **recovery module** (function = `recov.2strains`). XXX.

The **mortality module** (function = `deaths.2strains`). XXX.

XXX set the initial strains
XXX bookkeeping 
XXX plotting

### Parameters
The epidemic model parameters are:

* `inf.prob`: the probability that an infection will occur given an act between a susceptible node and a node who is infected with *strain 1* 
* `inf.prob.st2`: the probability that an infection will occur given an act between a susceptible node and a node who is infected with *strain 2*
* `rec.rate`: the rate of recovery for those infected with *strain 1*
* `rec.rate.st2`: the rate of recovery for those infected with *strain 2* 
* `di.rate`: the disease-specific mortality rate for those infected with *strain 1* 
* `di.rate.st2`: the disease-specific mortality rate for those infected with *strain 2* 

## Next Steps
XXX

## Author
Steven M. Goodreau, University of Washington (http://faculty.washington.edu/goodreau)