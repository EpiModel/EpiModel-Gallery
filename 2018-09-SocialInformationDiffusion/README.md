
Modeling Information Diffusion Process in a Social Network Using SI Model
=========================================================================

### Description

In this example we build an information diffusion model based on exsiting SI process in EpiModel to monitor dissemination of new ideas inside a social network. Instead of information diffusion (infection process) happening between any discordant nodes, we monitor two scenarios in which:

1.  The transmission needs a threshold of degree of discordant edgelist: it can only happen when suscpetible individuals have more than minimum degree of partnerships with infectious individuals;

2.  The probability of transmission is a function of degree of discordant partnerships.
    New argument includes minimum degree of discordant edgelist, log odds ratio of transmission with 1 degree increase of discordant partnership, baseline log odds when susceptible individuals have 0 discordant partners. No entirely new modules are needed, but infection process module is edited:

### Modules

**Scenario 1:**
The **infection module** (function = `infect_mod`) includes the following changes from the base EpiModel infection module (`infection.net`):

-   Book-keeping the degree of discordant edges for susceptible nodes

-   The infection probability is only assigned to susceptible nodes with more than minimum degree of discordant edges, otherwise is 0

**Scenario 2:**
The **infection module** (function = `infect_newmod`) includes the following changes from the base EpiModel infection module (`infection.net`):

-   Book-keeping the degree of discordant edges for susceptible nodes

-   The infection probability is a logistic function of degree of discordant edges for susceptible nodes with user specified paramters

### Parameters

The new or altered epidemic model parameters are:

**Scenario 1**:

-   `min.degree`: the minimum degree of discordant relationship for the infection to occur

**Scenario 2**:

-   `beta0`: baseline log odds of transmission when susceptible individuals have 0 degree of discordant partnership, usually negative as transmission probability is usually 0 when one has 0 degree of discordant partnership (suggest &lt;-3).

-   `beta1`: the increase of log odds with 1 degree increase of discordant partnership, usually set as positive as increase of degree can increase transmission probability.

![](coefs.png)

### Worked example

In the work example 1 in `model.R`, we consider a situation in which people are densely connected and susceptible individuals only accept new ideas if they have more than 3 partnerships with infected individuals (min.degree=3). Similar to other SI model, the infection process is slow initially and increases as more people become infected. However, the si.flow is lower than normal SI model and not as smooth as infectious probability jumps at minimum degree.

In the work example 2 in `model.R`, we consider a situation in which baseline transmission probability is almost 0 when susceptible individuals have 0 discordant partnership (log odds beta0=-7), and transmission probability increases relatively fast with the increase of the degree (log odds ratio beta1=0.5). The process is similar to infection process in normal SI model, incidence increases initially as more people are infected then decreases as the susceptible population decreases.

Author
------

Samuel Jenness, Emory University

Yuan Zhao, Emory University
