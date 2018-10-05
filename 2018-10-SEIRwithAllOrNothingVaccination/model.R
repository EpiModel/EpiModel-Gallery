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
n <- 50
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
mr_rate = 0.001

# Parameterize the dissolution model
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 60, d.rate = mr_rate)
coef.diss

# Fit the model
est <- netest(nw, formation, target.stats, coef.diss)

# Model diagnostics
dx <- netdx(est, nsims = 2, nsteps = 520,
            nwstats.formula = ~edges)
print(dx)
plot(dx, plots.joined = FALSE, qnts.alpha = 0.8)


# Epidemic model simulation -----------------------------------------------

# Model parameters
param <- param.net(inf.prob = 0.8,
                   mortality.rate = mr_rate,
                   mortality.disease.mult = 2,
                   act.rate = 3,
                   ei.rate = 0.001,
                   ir.rate = 0.001,
                   vaccine.rate = 0.1,
                   vaccine.efficacy = 0.1,
                   birth.rate = 0.001
                   )

# Initial conditions
init <- init.net(i.num = ceiling(0.1*n))

# Read in the module functions
source("C:/Users/conno/OneDrive/Documents/EpiModel Lab/SEIR with Vital Dynamics and All or Nothing Vaccination - Module.R", echo = TRUE)

# Control settings
control <- control.net(nsteps = 52 * 2,
                       nsims = 1,
                       infection.FUN = infect,
                       progress.FUN = progress,
                       recovery.FUN = NULL,
                       births.FUN = bfunc,
                       deaths.FUN = dfunc,
                       delete.nodes = TRUE,
                       depend = TRUE,
                       verbose = TRUE,
                       save.transmat = TRUE)

# Run the network model simulation with netsim
sim <- netsim(est, param, init, control)
print(sim)

# Plot outcomes
par(mar = c(3,3,1,1), mgp = c(2,1,0))
plot(sim, y = c("s.num","e.num","i.num","r.num", "vaccinated.num", "protected.num"),
     mean.col = 1:6, mean.lwd = 1, mean.smooth = FALSE,
     qnts = 1, qnts.col = 1:6, qnts.alpha = 0.25, qnts.smooth = FALSE,
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
head(nwdf, 25)

# A transmission matrix contains the time-ordered chain of transmissions
tm1 <- get_transmat(sim, sim = 1)
head(tm1, 10)

# Plotting with ggplot
df <- as.data.frame(sim, out = "vals")
df.mean <- as.data.frame(sim)

nw <- get_network(sim)
nw
nw <- color_tea(nw, verbose = FALSE)

slice.par <- list(start = 1, end = 50, interval = 1,aggregate.dur = 1, rule = "any")
render.par <- list(tween.frames = 10, show.time = FALSE)
plot.par <- list(mar = c(0, 0, 0, 0))
compute.animation(nw, slice.par = slice.par, verbose = TRUE)

nw <- color_tea(nw, old.var = "testatus", old.sus = "s", old.inf = "i",
          old.rec = "r", new.var = "ndtvcol", new.sus, new.inf, new.rec,
          verbose = FALSE)

render.d3movie(
  nw,
  render.par = render.par,
  plot.par = plot.par,
  vertex.col = c("blue","red"),
  edge.col = "darkgrey",
  vertex.border = "lightgrey",
  displaylabels = FALSE,
  filename = paste0(getwd(), "/movie.html"))
