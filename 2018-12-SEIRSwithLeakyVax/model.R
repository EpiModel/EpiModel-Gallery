##
## SEIRS Model with Vital Dynamics and a Leaky Vaccine Implementation
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##
## Authors: Connor M. Van Meter
## Date: January 2018
##

## Load EpiModel
suppressMessages(library(EpiModel))

# Standard Gallery unit test lines
rm(list = ls())
eval(parse(text = print(commandArgs(TRUE)[1])))

if (interactive()) {
  nsims <- 5
  ncores <- 5
  nsteps <- 500
} else {
  nsims <- 1
  ncores <- 1
  nsteps <- 50
}

# Network model estimation ------------------------------------------------

# Initialize the network
n <- 500
nw <- network_initialize(n)

# Define the formation model: edges
formation = ~edges

# Input the appropriate target statistics for each term
mean_degree <- 0.8
edges <- mean_degree * (n/2)

# Input the appropriate target statistics for each term
target.stats <- c(edges)

#Set departure rate
departure_rate <- 0.008

# Parameterize the dissolution model
coef.diss <- dissolution_coefs(dissolution = ~offset(edges),
                               duration = 50, d.rate = departure_rate)
coef.diss

# Fit the model
est <- netest(nw, formation, target.stats, coef.diss)

# Model diagnostics
dx <- netdx(est, nsims = nsims, ncores = ncores, nsteps = nsteps)

print(dx)
plot(dx)


# Epidemic model simulation -----------------------------------------------

# Model parameters
param <- param.net(inf.prob = 0.5,
                   arrival.rate = 0.01,
                   departure.rate = departure_rate,
                   departure.disease.mult = 2,
                   act.rate = 1,
                   ei.rate = 0.05,
                   ir.rate = 0.05,
				           rs.rate = 0.05,
                   vaccination.rate.initialization = 0.05,
                   protection.rate.initialization = 0.8,
                   vaccination.rate.progression = 0.05,
                   protection.rate.progression = 0.8,
                   vaccination.rate.arrivals = 0.6,
                   protection.rate.arrivals = 0.8,
                   vaccine.efficacy = 0.8)

# Initial conditions
init <- init.net(i.num = 20)

# Read in the module functions
if (interactive()) {
  source("2018-12-SEIRSwithLeakyVax/module-fx.R", echo = TRUE)
} else {
  source("module-fx.R")
}

# Control settings
control <- control.net(type = NULL,
                       nsteps = nsteps,
                       nsims = nsims,
                       ncores = ncores,
                       infection.FUN = infect,
                       progress.FUN = progress,
                       departures.FUN = dfunc,
                       arrivals.FUN = afunc,
                       resimulate.network = TRUE,
                       verbose = TRUE,
                       module.order = c("resim_nets.FUN", "departures.FUN",
                                        "arrivals.FUN", "infection.FUN",
                                        "progress.FUN", "prevalence.FUN"))

# Run the network model simulation with netsim
sim <- netsim(est, param, init, control)
print(sim)


# Examine the data from the simulation
df <- as.data.frame(sim)
df

#Epidemic plot of SEIRS compartment counts, entrances, and exits over simulation
par(mar = c(3,3,1,1), mgp = c(2,1,0))
plot(sim, y = c("e.num", "i.num", "r.num"),
     legend = TRUE)

# Standardized Incidence and Prevalence
sim <- mutate_epi(sim, ir.rate = se.flow / s.num,
                       prev = e.num / num)
plot(sim, y = c("ir.rate", "prev"), legend = TRUE)
