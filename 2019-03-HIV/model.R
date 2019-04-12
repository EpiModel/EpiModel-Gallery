##
## HIV Progression Model (One-Mode)
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##
## Authors: Samuel M. Jenness, Yuan Zhao, Connor Van Meter
## Date: March 2019
##

## Load EpiModel
suppressMessages(library(EpiModel))

# Standard Gallery unit test lines
rm(list = ls())
eval(parse(text = print(commandArgs(TRUE)[1])))

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
departure_rate = 0.008

# Parameterize the dissolution model
coef.diss <- dissolution_coefs(dissolution = ~offset(edges),
                               duration = 10, d.rate = departure_rate)
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
                   relative.inf.prob.final = 5,
                   relative.inf.prob.ART = 1/2,
                   arrival.rate = 0.005,
                   departure.rate = departure_rate,
                   departure.disease.mult = 2,
                   act.rate = 4,
                   AcuteToChronic1.Rate = 1/3,
                   Chronic1ToChronic2.Rate = 1/60,
                   Chronic2ToFinal.Rate = 1/60,
                   FinalToEnd.Rate = 1/24,
                   ART.Treatment.Rate = 0.10,
                   ART.Discontinuance.Rate = 0.05,
                   ART.Progression.Reduction.Rate = 0.5
                  )

# Initial conditions
init <- init.net(i.num = round(0.181 * n))


# Read in the module functions
source("module.R", echo = TRUE)

# Control settings
control <- control.net(nsteps = 735,
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


## Plot outcomes
#
# plot(sim, y = c("s.num","s.flow","i.num","i.flow"),
#      mean.col = 1:2, mean.lwd = 1, mean.smooth = TRUE,
#      qnts.col = 1:2, qnts.alpha = 0.25, qnts.smooth = TRUE,
#      legend = TRUE)
#
# plot(sim, y = c("s.num","s.flow","acute.num","acute.flow","chronic1.num","chronic1.flow","chronic2.flow","final.num","final.flow"),
#      mean.col = 1:5, mean.lwd = 1, mean.smooth = TRUE,
#      qnts.col = 1:5, qnts.alpha = 0.25, qnts.smooth = TRUE,
#      legend = TRUE)
#
# plot(sim, y = c("scr.flow"),
#      legend = TRUE)
#
# plot(sim, y = c("syph.dur","syph2.dur","syph3.dur","syph4.dur","syph5.dur","syph6.dur"),
#      mean.col = 1:6, mean.lwd = 1, mean.smooth = TRUE,
#      qnts.col = 1:6, qnts.alpha = 0.25, qnts.smooth = TRUE,
#      ylim = c(0,40),legend = TRUE)
#
# plot(sim, y = c("sym.num"),
#      mean.col = 1, mean.lwd = 1, mean.smooth = TRUE,
#      qnts.col = 1, qnts.alpha = 0.25, qnts.smooth = TRUE,
#      legend = TRUE)
#
# # Average across simulations at beginning, middle, end
# df <- as.data.frame(sim)
# df[c(2, 100, 400), ]
