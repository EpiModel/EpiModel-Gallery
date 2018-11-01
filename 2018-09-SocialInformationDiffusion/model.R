##
## Social Information Diffusion Model: Adding a minimum number of degree to S in an SI network
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##
## Author: Yuan Zhao
## Date: September 2018
##

# Load EpiModel
suppressMessages(library(EpiModel))
# Standard Gallery unit test lines
rm(list = ls())
eval(parse(text = print(commandArgs(TRUE)[1])))

# Network model estimation ------------------------------------------------

# Initialize the network
nw <- network.initialize(500, directed = FALSE)

# Define the formation model: edges + isolates (number with degree of 0)
formation = ~edges + isolates

# Input the appropriate target statistics for each term
target.stats <- c(300, 50)

# Parameterize the dissolution model
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 10)
coef.diss

# Fit the model
est <- netest(nw, formation, target.stats, coef.diss)

# Model diagnostics
dx <- netdx(est, nsims = 8, ncores = 8, nsteps = 500,
            nwstats.formula = ~edges + isolates + degree(0:5))
print(dx)
plot(dx, plots.joined = FALSE, qnts.alpha = 0.8)


# Epidemic model simulation -----------------------------------------------

# Model parameters
param <- param.net(inf.prob = 0.5, act.rate = 2,
                   ei.rate = 0.01, ir.rate = 0.01, min.degree=3)

# Initial conditions
init <- init.net(i.num = 150)

# Read in the module functions
source("module-fx.R", echo = TRUE)

# Control settings
control <- control.net(nsteps = 300,
                       nsims = 5,
                       ncores = 1,
                       infection.FUN = infect_mod
)

# Run the network model simulation with netsim
sim <- netsim(est, param, init, control)
print(sim)

# Plot outcomes
par(mar = c(3,3,1,1), mgp = c(2,1,0))
plot(sim,
     mean.col = 1:4, mean.lwd = 1, mean.smooth = FALSE,
     qnts = 1, qnts.col = 1:4, qnts.alpha = 0.25, qnts.smooth = FALSE,
     legend = TRUE)

plot(sim, y = c("si.flow"),
     mean.col = 1:4, mean.lwd = 1, mean.smooth = TRUE,
     qnts.col = 1:4, qnts.alpha = 0.25, qnts.smooth = TRUE,
     ylim = c(0, 3), legend = TRUE)

# Average across simulations at beginning, middle, 
df <- as.data.frame(sim)
df[c(3, 4, 5), ]
df[c(2, 100, 300), ]



# Plot network
par(mar = c(3,3,1,1), mgp = c(2,1,0))
plot(sim,type="network",at=10,sims="mean",
     col.status=TRUE, main="Prevalence at t10")
plot(sim,type="network",at=11,sims="mean",
     col.status=TRUE, main="Prevalence at t11")
plot(sim,type="network",at=200,sims="mean",
     col.status=TRUE, main="Prevalence at t200")

if (interactive() == TRUE) {
  library("ndtv")  
  nw <- get_network(sim)
  nw <- color_tea(nw, verbose = FALSE)
  slice.par <- list(start = 1, end = 25, interval = 1, 
                    aggregate.dur = 1, rule = "any")
  render.par <- list(tween.frames = 10, show.time = FALSE)
  plot.par <- list(mar = c(0, 0, 0, 0))
  compute.animation(nw, slice.par = slice.par, verbose = TRUE)
  render.d3movie(
    nw,
    render.par = render.par,
    plot.par = plot.par,
    vertex.col = "ndtvcol",
    edge.col = "darkgrey",
    vertex.border = "lightgrey",
    displaylabels = FALSE,
    filename = paste0(getwd(), "/movie.html"))
  
}
