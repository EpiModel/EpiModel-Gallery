# Test and Treat Intervention for an SIS Epidemic

## Description
This example demonstrates how to build an testing and treatment intervention for an SIS disease such as a bacterial infection like Gonorrhea. The main elements of this intervention itself are handled by the new `tnt` module while the recovery process (back to susceptibility) is handled with the updated (from the base, built-in) disease recovery module. 

### Modules
The **testing module** (function = `tnt`) handles the process of disease diagnosis on a routine interval for both infected and uninfected persons. One is eligible to test if currently in an undiagnosed state. The diagnosis state is reset back to `0` through either recovery from infection, or in the absence of recovery an interval of time passing (i.e., if infected, one never received treatment or treatment didn't work). 

The **recovery module** (function = `recov`) handles the process of recovery that now depends on the diagnosis status. In this model structure, one recovers from disease faster if diagnosed because receiving treatment. The main update from the built-in recovery module was to add this heterogenous recovery that was dependent on diagnosis status.

### Parameters
The epidemic model parameters are:

* `inf.prob`: the probability that an infection will occur given an act between a susceptible and infected node 
* `act.rate`: the number of acts per partnership per unit time 
* `rec.rate`: the rate of recovery for those who **are not** diagnosed 
* `rec.rate.tx`: the rate of recovery for those who **are** diagnosed
* `test.rate`: the rate at which everyone in the population who is not currently in a diagnosed state tests per time step 
* `test.dur`: the duration (in time steps) that the diagnosis status applies before getting reset back to `0`

## Next Steps
Next steps for this example might be to separate out the testing and treatment steps, or to vary the rate of testing or treatment based on an additional attribute of persons in the network.

## Author
Samuel M. Jenness, Emory University (http://samueljenness.org/)
