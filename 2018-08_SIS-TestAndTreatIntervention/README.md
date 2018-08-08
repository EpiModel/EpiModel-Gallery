# Test and Treat Intervention for an SIS Epidemic

## Description
This example demonstrates how to build an testing and treatment intervention for an SIS disease such as a bacterial infection like Gonorrhea. The main elements of this intervention itself are handled by the new `tnt` module while the recovery process (back to susceptibility) is handled with the updated (from the base, built-in) disease recovery module. 

The TnT intervention parameters are a rate of testing and a "duration" (in weeks) that the test applies. The rate of testing corresponds to the proportion of test-eligible persons who will, on average, test at the current time step. One is eligible to test if one has a current status of undiagnosed. That status gets reset either if one is infected and then recovers from disease, or automatically after a certain number of weeks, the test duration, because the diagnosis, either negative or positive, wouldn't last forever.

The recovery module handles the process of recovery that now depends on the diagnosis status. One may recover from disease presumably faster if diagnosed because all diagnosed persons would receive treatment.

## Next Steps
A good next step for this Gallery example might be to separate out the testing and treatment steps, or to vary the rate of testing or treatment based on an additional attribute of persons in the network.
