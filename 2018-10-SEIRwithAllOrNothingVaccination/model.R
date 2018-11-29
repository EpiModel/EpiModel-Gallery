##
## SEIR Model with Vital Dynamics and an All or Nothing Vaccine Implementation
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##
## Authors: Samuel M. Jenness, Venkata R. Duvvuri, Connor M. Van Meter
## Date: November 2018
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

#Set mortality rate
mr_rate = 0.008

# Parameterize the dissolution model
coef.diss <- dissolution_coefs(dissolution = ~offset(edges),
                               duration = 10, d.rate = mr_rate)
coef.diss

# Fit the model
est <- netest(nw, formation, target.stats, coef.diss)

# Model diagnostics
dx <- netdx(est, nsims = 3, nsteps = 52)

print(dx)
plot(dx)


# Epidemic model simulation -----------------------------------------------

# Model parameters
param <- param.net(inf.prob = 0.5,
                   birth.rate = 0.01,
                   mortality.rate = mr_rate,
                   mortality.disease.mult = 2,
                   act.rate = 1,
                   ei.rate = 0.05,
                   ir.rate = 0.05,

                   #Uncomment below for same rate of vaccination and protection across births, initialized network, and network progression
                   # vaccination.rate = 0.5,
                   # protection.rate = 0.5

                   #Keep uncommented for different rates of vaccination and protection across births, initialized network, and network progression
                   vaccination.rate.initialization = 0.05,
                   protection.rate.initialization = 0.8,
                   vaccination.rate.progression = 0.05,
                   protection.rate.progression = 0.8,
                   vaccination.rate.births = 0.6,
                   protection.rate.births = 0.8
)

# Initial conditions
init <- init.net(i.num = 20)

# Read in the module functions
source("module-fx.R", echo = TRUE)

# Control settings
control <- control.net(nsteps = 52,
                       nsims = 4,
                       ncores = 4,
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


##################################################################

# Examine the data from the simulation
df <- as.data.frame(sim)
df

df2 <- df[, c("time", "num")] # ...

#Data frame for SEIR-V compartment counts
df2 <- data.frame(df$time, df$num, df$s.num, df$e.num, df$i.num, df$r.num, df$v.num, df$b.num, df$b.flow, df$d.num, df$d.flow)
df2

#Data frame for vaccination flow and vaccination flow breakdown by vaccination method
df3 <- data.frame(df$vac.flow, df$vac.init.flow, df$vac.prog.flow, df$vac.birth.flow, df$vac.num, df$vac.init.num, df$vac.prog.num, df$vac.birth.num)
df3

#Data frame for vaccination protection flow and vaccination protection breakdown by vaccination method
df4 <- data.frame(df$prt.flow, df$prt.init.flow, df$prt.prog.flow, df$prt.birth.flow, df$prt.num, df$prt.init.num, df$prt.prog.num, df$prt.birth.num)
df4

#Epidemic plot of SEIR-V compartment counts, entrances, and exits over simulation
par(mar = c(3,3,1,1), mgp = c(2,1,0))
plot(sim, y = c("s.num","e.num","i.num","r.num", "v.num", "b.num", "d.num", "num"),
     mean.col = 1:8, mean.lwd = 1, mean.smooth = FALSE,
     qnts = 1, qnts.col = 1:8, qnts.alpha = 0.25, qnts.smooth = FALSE,
     legend = TRUE)
