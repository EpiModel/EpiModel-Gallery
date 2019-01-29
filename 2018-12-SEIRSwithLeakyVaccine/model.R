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


# Network model estimation ------------------------------------------------

# Initialize the network
n <- 100
nw <- network.initialize(n, directed = FALSE)

# Define the formation model: edges
formation = ~edges

# Input the appropriate target statistics for each term
mean_degree <- 0.8
edges <- mean_degree * (n/2)

# Input the appropriate target statistics for each term
target.stats <- c(edges)

#Set departure rate
departure_rate = 0.008

# Parameterize the dissolution model
coef.diss <- dissolution_coefs(dissolution = ~offset(edges),
                               duration = 20, d.rate = departure_rate)
coef.diss

# Fit the model
est <- netest(nw, formation, target.stats, coef.diss)

# Model diagnostics
dx <- netdx(est, nsims = 10, nsteps = 520)

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
                   vaccine.efficacy = 0.8
)

# Initial conditions
init <- init.net(i.num = 20)

# Read in the module functions
source("module-fx.R", echo = TRUE)

# Control settings
control <- control.net(nsteps = 520,
                       nsims = 1,
                       ncores = 1,
                       infection.FUN = infect,
                       progress.FUN = progress,
                       recovery.FUN = NULL,
                       arrivals.FUN = afunc,
                       departures.FUN = dfunc,
                       delete.nodes = TRUE,
                       depend = TRUE,
                       verbose = TRUE)

# Run the network model simulation with netsim
sim <- netsim(est, param, init, control)
print(sim)


##################################################################

# Examine the data from the simulation
df <- as.data.frame(sim)
df

#Data frame for SEIR-V compartment counts
df2 <- df[, c("time", "num", "s.num", "se.flow", "e.num", "ei.flow", "i.num",
              "ir.flow", "r.num", "rs.flow", "a.num", "a.flow", "d.num",
              "d.flow")]
df2

#Data frame for vaccination flow and vaccination flow breakdown
#by vaccination method
df3 <- df[, c("vac.flow", "vac.init.flow", "vac.prog.flow", "vac.arrival.flow",
              "vac.num", "vac.init.num", "vac.prog.num", "vac.arrival.num")]
df3

#Data frame for vaccination protection flow and vaccination protection breakdown
#by vaccination method
df4 <- df[, c("prt.flow", "prt.init.flow", "prt.prog.flow", "prt.arrival.flow",
              "prt.num", "prt.init.num", "prt.prog.num", "prt.arrival.num")]
df4

#Epidemic plot of SEIRS compartment counts, entrances, and exits over simulation
par(mar = c(3,3,1,1), mgp = c(2,1,0))
plot(sim, y = c("s.num","e.num","i.num","r.num", "a.num", "d.num", "num"),
     mean.col = 1:7, mean.lwd = 1, mean.smooth = FALSE,
     qnts = 1, qnts.col = 1:7, qnts.alpha = 0.25, qnts.smooth = FALSE,
     legend = TRUE)

#Cumulative Incidence and Prevalence Plots of the Vaccine Model
sim <- mutate_epi(sim, ci = se.flow / s.num, prev = e.num / num)
plot(sim, y = c("ci", "prev"), mean.lwd = 1, mean.smooth = TRUE, legend = TRUE)


##################################################################
#VACCINE MODEL COMPARISON
##################################################################

# Update one or more vaccine parameters, simulate the new model,
# and compare with the original model

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
                   vaccine.efficacy = 0.1
)

# Initial conditions
init <- init.net(i.num = 20)

# Read in the module functions
source("module-fx.R", echo = TRUE)

# Control settings
control <- control.net(nsteps = 520,
                       nsims = 1,
                       ncores = 1,
                       infection.FUN = infect,
                       progress.FUN = progress,
                       recovery.FUN = NULL,
                       arrivals.FUN = afunc,
                       departures.FUN = dfunc,
                       delete.nodes = TRUE,
                       depend = TRUE,
                       verbose = TRUE)

# Run the network model simulation with netsim
sim2 <- netsim(est, param, init, control)
print(sim2)

##Compare incidence and prevalence of simulation 1 to simulation 2

#Calculate cumulative incidence and prevalence of simulation 2
sim2 <- mutate_epi(sim2, ci2 = se.flow / s.num, prev2 = e.num / num)

par(mfrow = c(1,1))
plot(sim, y = c("ci", "prev"), mean.lwd = 1, mean.smooth = TRUE, legend = TRUE)
plot(sim2, y = c("ci2", "prev2"), mean.lwd = 1, mean.smooth = TRUE, add = TRUE,
     mean.col = c("steelblue", "firebrick"), qnts.col = c("steelblue",
                                                          "firebrick"))
