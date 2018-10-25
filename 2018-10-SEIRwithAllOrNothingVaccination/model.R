##
## SEIR Model with Vital Dynamics and an All or Nothing Vaccine Implementation
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##
## Authors: Samuel M. Jenness, Venkata R. Duvvuri, Connor M. Van Meter
## Date: October 2018
##

## Load EpiModel
suppressMessages(library(EpiModel))

# Standard Gallery unit test lines
rm(list = ls())
eval(parse(text = print(commandArgs(TRUE)[1])))


# Network model estimation ------------------------------------------------

# Initialize the network
n <- 20
nw <- network.initialize(n, directed = FALSE)


##Set up vaccination and protection charateristics for initial nodes of network

#Example starting condition below is that all nodes of the initialized network
#are "unvaccinated" (as opposed to "vaccinated") and "vulnerable" (as opposed to "protected").
#Note: Vaccination and protection can only be conferred at birth - extension of the birth module
# vaccinatedVec <- sample(rep("unvaccinated", n))
# protectedVec <- sample(rep("vulnerable", n))

#Note: "Vaccinated" and "protected" are attribute names used in the module
# nw <- set.vertex.attribute(nw, "vaccinated", vaccinatedVec)
# nw <- set.vertex.attribute(nw, "protected", protectedVec)

# Define the formation model: edges
formation = ~edges

# Input the appropriate target statistics for each term
mean_degree <- 0.8
edges <- mean_degree * (n/2)

# Input the appropriate target statistics for each term
target.stats <- c(edges)

#Set mortality rate
mr_rate = 0.000

# Parameterize the dissolution model
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 40, d.rate = mr_rate)
coef.diss

# Fit the model
est <- netest(nw, formation, target.stats, coef.diss)

# Model diagnostics
dx <- netdx(est, nsims = 20, nsteps = 3,
            nwstats.formula = ~edges)
print(dx)
plot(dx, plots.joined = FALSE, qnts.alpha = 0.8)


# Epidemic model simulation -----------------------------------------------

# Model parameters
param <- param.net(inf.prob = 0,
                   mortality.rate = mr_rate,
                   mortality.disease.mult = 2,
                   act.rate = 1,
                   ei.rate = 0,
                   ir.rate = 0,
                   #vaccination.rate = 0.4,
                   #protection.rate = 0.8,
                   vaccination.rate.initialization = 0.4,
                   protection.rate.initialization = 0.5,
                   vaccination.rate.progression = 0.5,
                   protection.rate.progression = 0.5,
                   vaccination.rate.births = 0.5,
                   protection.rate.births = 0.5,
                   birth.rate = 0.01
                   )

# Initial conditions
init <- init.net(i.num = 4)

# Read in the module functions
source("module-fx.R", echo = TRUE)

# Control settings
control <- control.net(nsteps = 50,
                       nsims = 1,
                       infection.FUN = infect,
                       progress.FUN = progress,
                       recovery.FUN = NULL,
                       births.FUN = bfunc,
                       deaths.FUN = dfunc,
                       delete.nodes = TRUE,
                       depend = TRUE,
                       verbose = TRUE)

# Run the network model simulation with netsim
sim <- netsim(est, param, init, control)
print(sim)


##PLOTS AND TESTING

# Plot outcomes
par(mar = c(3,3,1,1), mgp = c(2,1,0))
plot(sim, y = c("s.num","e.num","i.num","r.num", "vaccinated.num", "protected.num", "b.flow", "d.flow"),
     mean.col = 1:8, mean.lwd = 1, mean.smooth = FALSE,
     qnts = 1, qnts.col = 1:8, qnts.alpha = 0.25, qnts.smooth = FALSE,
     legend = TRUE)

plot(sim, y = c("vaccinated.num"),
     mean.col = 1:8, mean.lwd = 1, mean.smooth = FALSE,
     qnts = 1, qnts.col = 1:8, qnts.alpha = 0.25, qnts.smooth = FALSE,
     legend = TRUE)

plot(sim, y = c("b.flow","se.flow", "ei.flow", "ir.flow"),
     mean.col = 1:6, mean.lwd = 1, mean.smooth = TRUE,
     qnts.col = 1:6, qnts.alpha = 0.25, qnts.smooth = TRUE,
     #ylim = c(0, 3),
     legend = TRUE)

plot(sim, y = c("vaccinated.flow", "protected.flow"),
     mean.col = 1:6, mean.lwd = 1, mean.smooth = TRUE,
     qnts.col = 1:6, qnts.alpha = 0.25, qnts.smooth = TRUE,
     #ylim = c(0, 3),
     legend = TRUE)

par(mfrow = c(1, 2))
plot(sim, y = "num", main = "Population Size", qnts = 1, ylim = c(450, 550))
plot(sim, y = "meanAge", main = "Mean Age", qnts = 1, ylim = c(35, 50))

par(mfrow = c(1, 2))
plot(sim, y = "d.flow", mean.smooth = TRUE, qnts = 1, main = "Deaths")
plot(sim, y = "b.flow", mean.smooth = TRUE, qnts = 1, main = "Births")

# Examine the average across simulations at beginning, middle, end
df <- as.data.frame(sim)
head(df, 25) #First 25 steps (weeks)
df[c(2, 50, 100), ] #Weeks 2, 50, 100

df2 <- data.frame(df$time, df$num, num_2=df$s.num+df$e.num+df$i.num+df$r.num, df$s.num, df$e.num, df$i.num, df$r.num, df$d.flow, df$b.flow, df$se.flow)
tail(df2,1)
df2

# Convert model to a data frame for further analysis
# Default conversion is means across simulations
df <- as.data.frame(sim)
head(df, 10)

# Extracting individual values also possible
df <- as.data.frame(sim, out = "vals")
head(df, 10)
tail(df, 10)

# Extract the full dynamic network for further analysis
nw1 <- get_network(sim, sim = 1)
nw1

# Temporal edgelist
nwdf <- as.data.frame(nw1)
tail(nwdf, 25)


