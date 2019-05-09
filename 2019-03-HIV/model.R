##
## HIV Transmission Model
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##
## Authors: Connor Van Meter, Samuel M. Jenness, Yuan Zhao
## Date: March 2019
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
n <- 1000
nw <- network.initialize(n, directed = FALSE)

# Define the formation model: edges
formation = ~edges

# Input the appropriate target statistics for each term
mean_degree <- 0.8
edges <- mean_degree * (n/2)

target.stats <- c(edges)

#Set departure rate
departure_rate = 0.003

# Parameterize the dissolution model
coef.diss <- dissolution_coefs(dissolution = ~offset(edges),
                               duration = 52, d.rate = departure_rate)
coef.diss

# Fit the model
est <- netest(nw, formation, target.stats, coef.diss)

# Model diagnostics
dx <- netdx(est, nsims = 1, nsteps = 735)

print(dx)
plot(dx, plots.joined = FALSE, qnts.alpha = 0.8)


# Epidemic model simulation -----------------------------------------------

# Model parameters
param <- param.net(inf.prob.chronic = 0.01,
                   relative.inf.prob.acute = 10,
                   relative.inf.prob.AIDS = 5,
                   relative.inf.prob.ART = 0.05,
                   act.rate = 4,
                   AcuteToChronic1.Rate = 1/12,
                   Chronic1ToChronic2.Rate = 1/260,
                   Chronic2ToAIDS.Rate = 1/260,
                   AIDSToDepart.Rate = 1/104,
                   ART.Treatment.Rate = 0.10,
                   ART.Discontinuance.Rate = 0.05,
                   ART.Progression.Reduction.Rate = 0.5,
                   arrival.rate = 0.005,
                   departure.rate = departure_rate,
                   departure.disease.mult = 2)

# Initial conditions
start_prevalence = 0.05
init <- init.net(i.num = round(start_prevalence * n))


# Read in the module functions
source("module-fx.R", echo = TRUE)

# Control settings
control <- control.net(nsteps = nsteps,
                       nsims = nsims,
                       ncores = ncores,
                       infection.FUN = infect,
                       progress.FUN = progress,
                       recovery.FUN = NULL,
                       arrivals.FUN = afunc,
                       departures.FUN = dfunc,
                       depend = TRUE,
                       verbose = TRUE)

# Run the network model simulation with netsim
sim <- netsim(est, param, init, control)
print(sim)


# Examine the data from the simulation
df <- as.data.frame(sim)
df


## Plot outcomes

# SI Compartment Counts
par(mar = c(2,2,1,1), mgp = c(2,1,0))
plot(sim, y = c("s.num","i.num"),
     mean.col = 1:2, mean.lwd = 1, mean.smooth = TRUE,
     qnts.col = 1:2, qnts.alpha = 0.25, qnts.smooth = TRUE,
     legend = TRUE)

# HIV status sub-compartment counts
par(mar = c(5,5,1,1), mgp = c(2,1,0))
plot(sim, y = c("acute.num", "chronic1.num", "chronic2.num",
                "AIDS.num"),
     mean.col = 1:4, mean.lwd = 1, mean.smooth = TRUE,
     legend = TRUE)


# Standardized Incidence and Prevalence
sim <- mutate_epi(sim, ir.rate = si.flow / s.num,
                  prev = i.num / num)
par(mfrow = c(1,2))
plot(sim, y = "prev", main = "Prevalence")
plot(sim, y = "ir.rate", main = "Incidence")
