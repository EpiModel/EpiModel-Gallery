# HIV Model

## Description
This example shows how to model a simple, one-mode HIV transmission process over a dynamic network with vital dynamic processes for population arrival and departure.  
For a detailed diagram of the model, please reference the accompanying HIV Model Diagram PowerPoint.

EpiModel includes an integrated SI model, but here we show how to model a one-mode SI disease, in this case a simple HIV model, with four distinct infectious - I - sub-compartments, based on the deterministic transmission model suggested by Granich et al in their 2006 Lancet paper [Universal voluntary HIV testing with immediate antiretroviral therapy as a strategy for elimination of HIV transmission: a mathematical model](https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(08)61697-9/fulltext "Granich et al HIV Model").  
The four infectious phases of the proposed model include:  
**Acute phase** - occurs immediately after transmission and is characterized by high viral load and relatively high infectivity  
**Chronic phases** - chronic 1 and chronic 2, characterized by similarly low relative infectivity but distinguished by declining CD4+ count in the chronic 2 phase after reaching an apex in the chronic 1 phase  
**Final phase** - declining CD4+ counts with rising viral load lead to an increase in relative infectivity a decline in survival  

The EpiModel progression module has been updated to allow individuals to stochastically transition between these phases.  

The Granich model includes an ART treatment intervention in which each infectious phase subcompartment contains a corresponding compartment for individuals in those infectious phases being treated with ART. Disease progression is slowed for individuals receiving ART treatment.  
The progression module has also incorporated a stochastic ART treatment intervention process in which simulated individuals in each of the four infectious phases of HIV and not already on ART may be stochastically selected for ART treatment. Individuals receiving ART treatment have decreased infectiousness and have slower transition to subsequent HIV infectious phases than similar individuals not receiving ART treatment. Similar to how individuals not already receiving ART may be stochastically selected to start receiving ART treatment, individuals receiving ART may be stochastically selected to discontinue ART treatment. Separate subcompartments have been created for infectious individuals receiving ART treatment.

The arrival module implements a stochastic entrance process that acts as a function of a standard arrival rate.

The departure module implements a stochastic exit process that acts as a function of either a standard departure rate or a disease-induced departure rate.

### Modules
The **infection module** simulates HIV transmission from persons living with HIV to persons not living with HIV. Disease transmission is a function of the number of acts between persons, the phase progression of the individual living with HIV - acute- and final-phase individuals have higher likelihood of transmission than chronic-phase individuals, whether or not the individual living with HIV is receiving ART treatment, and the per act transmission probability.

The **progression module** simulates both disease phase progression for an individual living with HIV and ART treatment and discontinuance processes. After disease transmission, individuals progress from the acute phase of HIV to the chronic 1 phase, fro chronic 1 to the chronic 2 phase, and from chronic 2 to the final phase. Individuals living with HIV and not receiving ART may be selected to receive ART treatment. Once on ART, disease progression is slowed and transmission probability is reduced. Similarly, the module simulates individuals discontinuing their ART treatment in which case their disease progression returns to pre-ART rates and transmission probability is no longer reduced.

The **departure module** (function = `dfunc`)  simulates departure from the model as a function of a disease-induced departure rate. The standard departure rate is passed in by the module user as a rate - `departure.rate` - in the epidemic model parameter settings. For those eligible (active) individuals whose disease status is `"i"`, the departure rate is multiplied by the value of the `departure.disease.mult` parameter representing an increased likelihood of model departure for infected indviduals. Persons living with final-phase HIV depart the model at different, often higher, rates than other individuals; as a result, a separate departure process has been established for individuals living with final-stage HIV mimicking disease phase progression.

The **arrival module** (function = `afunc`) simulates arrival into the model as a function of the network size.


### Parameters
The epidemic model parameters include those needed for establishing the infrastructure for the simple HIV model transmisison and progression and those pertaining to ART treatment and discontinuance.

* `inf.prob.chronic`: the probability that transmission may occur given an act between a susceptible individual and an individual living with chronic-phase HIV  
* `relative.inf.prob.acute`: a scalar multiplier for the increased relative transmissibility that may occur given an act between a susceptible individual and an individual living with acute-phase HIV  
* `relative.inf.prob.final`: a scalar multiplier for the increased relative transmissibility that may occur given an act between a susceptible individual and an individual living with final-phase HIV  
* `relative.inf.prob.ART`: a scalar multiplier for the decreased relative transmissibility that may occur given an act between a susceptible individual and an individual living with HIV but receiving ART treatment  
* `act.rate`: the number of acts per partnership per unit time  
* `AcuteToChronic1.Rate`: the rate of persons living with acute HIV moving to the chronic 1 phase (1/average duration spent in the acute phase of HIV)  
* `Chronic1ToChronic2.Rate`: the rate of persons living with chronic HIV, first phase, moving to the second phase (1/average duration spent in the chronic 1 phase of HIV)  
* `Chronic2ToFinal.Rate`: the rate of persons living with chronic HIV, second phase, moving to the final phase (1/average duration spent in the chronic 2 phase of HIV)  
* `FinalToDepart.Rate`: the rate of persons living with final phase HIV departing the model (1/average duration spent in the final phase of HIV)  
* `ART.Treatment.Rate`: the rate in which persons living with HIV and not already on ART receive ART treatment  
* `ART.Discontinuance.Rate`: the rate in which persons living with HIV and receiving ART discontinue their treatment  
* `ART.Progression.Reduction.Rate`: a scalar multiplier for the decreased relative progression to subsequent HIV infectious phases for individuals receiving ART treatment  
* `arrival.rate`: a scalar for the rate of arrivals per person per week  
* `departure.rate`: the standard departure rate of the population. For disease status of `"i"`, this rate is multiplied by the `departure.disease.mult` explained below  
* `departure.disease.mult`: a scalar multiplier for the increased relative risk of departure for persons with a disease status of `"i"'  


## Authors
Samuel M. Jenness, Yuan Zhao, Connor Van Meter
