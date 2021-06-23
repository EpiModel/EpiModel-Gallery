
##
## Simple SI Model with Vital Dynamics
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##
## Author: Samuel M. Jenness (Emory University)
## Date: August 2018
##

# Load EpiModel
suppressMessages(library(EpiModel))

# Standard Gallery unit test lines
rm(list = ls())
eval(parse(text = print(commandArgs(TRUE)[1])))

if (interactive()) {
  nsims <- 5
  ncores <- 5
  nsteps <- 256
} else {
  nsims <- 1
  ncores <- 1
  nsteps <- 52
}


# Vital Dynamics Setup ----------------------------------------------------

ages <- 0:85

# Rates per 100,000 for age groups: <1, 1-4, 5-9, 10-14, 15-19, 20-24, 25-29,
#                                   30-34, 35-39, 40-44, 45-49, 50-54, 55-59,
#                                   60-64, 65-69, 70-74, 75-79, 80-84, 85+
# source: https://www.statista.com/statistics/241572/death-rate-by-age-and-sex-in-the-us/
departure_rate <- c(588.45, 24.8, 11.7, 14.55, 47.85, 88.2, 105.65, 127.2,
                    154.3, 206.5, 309.3, 495.1, 736.85, 1051.15, 1483.45,
                    2294.15, 3642.95, 6139.4, 13938.3)
# rate per person, per week
dr_pp_pw <- departure_rate / 1e5 / 52

# Build out a mortality rate vector
age_spans <- c(1, 4, rep(5, 16), 1)
dr_vec <- rep(dr_pp_pw, times = age_spans)
data.frame(ages, dr_vec)

par(mar = c(3,3,2,1), mgp = c(2,1,0), mfrow = c(1,1))
barplot(dr_vec, col = "steelblue1", xlab = "age", ylab = "Departure Rate")


# Network Model Estimation ------------------------------------------------

# Initialize the network
n <- 500
nw <- network_initialize(n)

# Set up ages
ageVec <- sample(ages, n, replace = TRUE)
nw <- set_vertex_attribute(nw, "age", ageVec)

# Define the formation model: edges
formation <- ~edges + absdiff("age")

# Input the appropriate target statistics for each term
mean_degree <- 0.8
edges <- mean_degree * (n/2)
avg.abs.age.diff <- 1.5
absdiff <- edges * avg.abs.age.diff

target.stats <- c(edges, absdiff)

# Parameterize the dissolution model
coef.diss <- dissolution_coefs(~offset(edges), 60, mean(dr_vec))
coef.diss

# Fit the model
est <- netest(nw, formation, target.stats, coef.diss)

# Model diagnostics
dx <- netdx(est, nsims = nsims, ncores = ncores, nsteps = nsteps,
            nwstats.formula = ~edges + absdiff("age") + isolates + degree(0:5))
print(dx)
plot(dx)


# Epidemic model simulation -----------------------------------------------

# Epidemic model parameters
param <- param.net(inf.prob = 0.15,
                   departure.rates = dr_vec,
                   departure.disease.mult = 2,
                   arrival.rate = mean(dr_vec))

# Initial conditions
init <- init.net(i.num = 50)

# Read in the module functions
if (interactive()) {
  source("2018-08-SIwithVitalDynamics/module-fx.R", echo = TRUE)
} else {
  source("module-fx.R")
}

# Control settings
control <- control.net(type = NULL,
                       nsims = nsims,
                       ncores = ncores,
                       nsteps = nsteps,
                       aging.FUN = aging,
                       departures.FUN = dfunc,
                       arrivals.FUN = afunc,
                       infection.FUN = infection.net,
                       resim_nets.FUN = resim_nets,
                       resimulate.network = TRUE,
                       verbose = FALSE)

# Run the network model simulation with netsim
sim <- netsim(est, param, init, control)
print(sim)

# Plot outcomes
par(mfrow = c(1,2))
plot(sim, main = "State Prevalences", popfrac = TRUE)
plot(sim, main = "State Sizes", sim.lines = TRUE,
     qnts = FALSE, mean.smooth = FALSE)

par(mfrow = c(1, 2))
plot(sim, y = "num", main = "Population Size", qnts = 1, ylim = c(450, 550))
plot(sim, y = "meanAge", main = "Mean Age", qnts = 1, ylim = c(35, 50))

par(mfrow = c(1, 2))
plot(sim, y = "d.flow", mean.smooth = TRUE, qnts = 1, main = "Departures")
plot(sim, y = "a.flow", mean.smooth = TRUE, qnts = 1, main = "Arrivals")

# Examine the data
df <- as.data.frame(sim, out = "mean")
head(df, 25)
