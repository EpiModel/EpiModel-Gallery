# Simple Cost-effectiveness Model (Simple SI model with cost and utility tracking)

## Description
This example shows how to use EpiModel to conduct a cost-effectiveness analysis in the context of an SI epidemic within a dynamic population undergoing birth, aging, and death. The two competing strategies compared in this cost-effectiveness analysis are 1) a baseline scenario with no intervention in place and 2) a universal prophylaxis intervention where the probability of infection per discordant sex act is halved. It is assumed that the clinical care accrued by healthy individuals is less than that of infected individuals, and by reducing the rate of infection, the prophylaxis intervention both improves health (in terms of quality adjusted life years, or QALYs) and reduces clinical care costs. However, the implementation costs of the intervention are substantial and must be weighed against the program's benefits. Given these trade-offs, we seek to answer the question of whether the prophylaxis intervention represents a cost-effective investment of resources from a population-level health care perspective.

There are a number of adjustments that must be made to the estimation of the initial network as well as the simulation modules in order to conduct a cost-effectiveness analysis. These changes are described here only broadly, but further details can be found within comments surrounding the example code.

Firstly, the distinction must be made between egos who are sexually active (attribute variable `active.s`) versus those who are alive (attribute variable `active`). Individuals who are sexually inactive will not participate in the sexual network, but they should still accrue clinical care costs and QALYs until death. Individuals over age 65 are assumed to be "sexually retired", and the fitted sexual network must prevent the formation of sexual partnerships with such individuals.

Secondly, we must carefully consider the time horizon of our analysis and specify the time-step at which tracking of costs and effects begins, which usually aligns with the initiation of the intervention(s) in question. We also must specify a final time-step for the analytic time horizon, when all processes related to disease transmission, population entry, and sexual dynamics are stopped. In order to fully capture the health and costs of individuals still living at the end of this main time horizon, we continue the simulation and cost/QALY tracking until all individuals have died. These residual costs and QALYs of individuals still living at the end of the main time horizon are known as "end horizon" effects. In the code that accompanies this example, the simulation is only run only for a brief section of the end horizon and not until all individuals have died in order to reduce computing time.

Thirdly, depending on the complexity of the functions defining costs and effects for simulated individuals, a cost-effectiveness module may need to be added to the simulation. Under certain circumstances, it is possible to calculate all costs and effects externally to the simulation using only the `sim$epi` output. Both methods are employed for the purpose of demonstration in `model.R`.

### Modules
The **cost-effectiveness module** (function = `costeffect`) calculates the costs and QALYs accrued by individuals along with relevant intervention costs at each time-step. Individual costs and QALYs are functions of attributes within the model.

The **aging module** (function = `aging`) sets up the age attribute at the initial time step by pulling from the network object. Then it will subsequently update the age attribute in increments of a week at each time step for individuals who are still alive (attribute `active == 1`).

The **death module** (function = `dfunc`)  simulates mortality as a function of an age-specific mortality rate. Forces egos into sexual retirement at age 65. Egos who sexually retire but do not die continue to be tracked in cost-effectiveness module.

The **arrival module** (function = `afunc`) implements an updated population arrival process. This module is disabled upon reaching the end horizon.

The **infection module** (function = `ifunc`) implements SI disease transmission. Transmission is disabled for lingering partnerships involving a sexually inactive ego, and this module is disabled upon reaching the end horizon.

The **network resimulation module** (function = `resimfunc`) resimulates the sexual network at each timestep as individuals sexually retire, die, or enter the population. This module is identical to the default module except that it is disabled upon reaching the end horizon, greatly accelerating computational speed.

### Parameters
The epidemic model parameters are basic here because we're not changing any of the core epidemiology from a simple SI model.

* `cea.start`: the time-step at which the accounting for costs and effects begins to take place.

* `end.horizon`: the time-step at which the main analytic time horizon ends and the accounting of "end horizon" effects begins.

* `d_r`: the annual rate at which costs and QALYs are discounted according to when they took place relative to `cea.start`. The general rationale for discounting in CEA is that costs and benefits that are deferred have lower value than those that are realized immediately. 

* `inter.eff`: 1 - the multiplicative reduction in probability that an infection will occur given an act in a sero-discordant partnership under the prophylaxis intervention. 

* `inter.start`: the time-step at which the prophylaxis intervention begins.

* `inter.cost`: the total cost of the intervention, spread evenly between weekly time-steps from `inter.start` to `end.horizon`.

* `sus.cost`: the weekly cost ($) of clinical care for a healthy (susceptible) individual.

* `inf.cost`: the weekly cost ($) of clinical care for an infected individual.

* `sus.qaly`: the QALYs accrued by a healthy (susceptible) individual over one year.

* `inf.qaly`: the QALYs accrued by a infected individual over one year.

* `age.decrement`: the additive reduction in QALYs accrued over one year for each year of the individual's age.

## Author
Greg Knowlton, University of Minnesota
