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
n <- 100
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
mr_rate = 0.005

# Parameterize the dissolution model
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 10, d.rate = mr_rate)
coef.diss

# Fit the model
est <- netest(nw, formation, target.stats, coef.diss)

# Model diagnostics
dx <- netdx(est, nsims = 3, nsteps = 52,
            nwstats.formula = ~edges)
print(dx)
plot(dx, plots.joined = FALSE, qnts.alpha = 0.8)


# Epidemic model simulation -----------------------------------------------

# Model parameters
param <- param.net(inf.prob = 0.08,
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
init <- init.net(i.num = 4)

# Read in the module functions
source("module-fx.R", echo = TRUE)

# Control settings
control <- control.net(nsteps = 52,
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
# par(mar = c(3,3,1,1), mgp = c(2,1,0))
# plot(sim, c("s.num","e.num","i.num","r.num", "v.num", "b.num", "num"), type="epi",
#      mean.col = 1:7, mean.lwd = 1, mean.smooth = FALSE,
#      qnts = 1, qnts.col = 1:7, qnts.alpha = 0.25, qnts.smooth = FALSE,
#      legend = TRUE)
#
# plot(sim, c("v.num"), type="epi",
#      mean.col = 1:8, mean.lwd = 1, mean.smooth = FALSE,
#      qnts = 1, qnts.col = 1:8, qnts.alpha = 0.25, qnts.smooth = FALSE,
#      legend = TRUE)

# plot(sim, y = c("b.flow","se.flow", "ei.flow", "ir.flow"),
#      mean.col = 1:6, mean.lwd = 1, mean.smooth = TRUE,
#      qnts.col = 1:6, qnts.alpha = 0.25, qnts.smooth = TRUE,
#      #ylim = c(0, 3),
#      legend = TRUE)
#
# plot(sim, y = c("vaccinated.flow", "protected.flow"),
#      mean.col = 1:6, mean.lwd = 1, mean.smooth = TRUE,
#      qnts.col = 1:6, qnts.alpha = 0.25, qnts.smooth = TRUE,
#      #ylim = c(0, 3),
#      legend = TRUE)
#
# par(mfrow = c(1, 2))
# plot(sim, y = "num", main = "Population Size", qnts = 1, ylim = c(450, 550))
# plot(sim, y = "meanAge", main = "Mean Age", qnts = 1, ylim = c(35, 50))
#
# par(mfrow = c(1, 2))
# plot(sim, y = "d.flow", mean.smooth = TRUE, qnts = 1, main = "Deaths")
# plot(sim, y = "b.flow", mean.smooth = TRUE, qnts = 1, main = "Births")

# Examine the average across simulations at beginning, middle, end
df <- as.data.frame(sim)
head(df, 25) #First 25 steps (weeks)

#Data frame for SEIR-V numbers and vaccination flow
df2 <- data.frame(df$time, df$num, num_2=df$s.num+df$e.num+df$i.num+df$r.num+df$v.num, df$s.num, df$e.num, df$i.num, df$r.num, df$v.num)
df2

#Testing vaccination and protection vectors ('All-Or-Nothing Vaccine Model - Data Profiling.xlsx')
#Check average over two timesteps
# df_t2 <- subset(df,df[,2] == 2)
# means_t2 <- c(mean(df_t2$s.num),mean(df_t2$e.num),mean(df_t2$i.num),mean(df_t2$r.num),mean(df_t2$v.num), mean(df_t2$num))
# means_t2
# colMeans(df_t2)

#Check average over three timesteps
# df_t3 <- subset(df,df[,2] == 3)
# means_t3 <- c(mean(df_t3$s.num),mean(df_t3$e.num),mean(df_t3$i.num),mean(df_t3$r.num),mean(df_t3$v.num), mean(df_t3$num))
# means_t3
# colMeans(df_t3)


#
# # Convert model to a data frame for further analysis
# # Default conversion is means across simulations
# df <- as.data.frame(sim)
# head(df, 10)
#
# # Extracting individual values also possible
# df <- as.data.frame(sim, out = "vals")
# head(df, 10)
# tail(df, 10)
#
# # Extract the full dynamic network for further analysis
# nw1 <- get_network(sim, sim = 1)
# nw1
#
# # Temporal edgelist
# nwdf <- as.data.frame(nw1)
# tail(nwdf, 25)
