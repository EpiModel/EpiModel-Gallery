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
                   relative.inf.prob.final = 5,
                   relative.inf.prob.ART = 1/2,
                   act.rate = 4,
                   AcuteToChronic1.Rate = 1/3,
                   Chronic1ToChronic2.Rate = 1/60,
                   Chronic2ToFinal.Rate = 1/60,
                   FinalToDepart.Rate = 1/24,
                   ART.Treatment.Rate = 0.10,
                   ART.Discontinuance.Rate = 0.05,
                   ART.Progression.Reduction.Rate = 0.5,
                   arrival.rate = 0.005,
                   departure.rate = departure_rate,
                   departure.disease.mult = 2)

# Initial conditions
start_prevalence = 0.181
init <- init.net(i.num = round(start_prevalence * n))


# Read in the module functions
source("module-fx.R", echo = TRUE)

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



# Examine the data from the simulation
df <- as.data.frame(sim)
df

#Data frame for overall SI, arrival and departure counts
df2 <- df[, c("time", "num",
              "s.num", "i.num", "si.flow",
              "a.flow",
              "d.flow")]
df2

#Data frame for HIV status sub-compartment counts
df3 <- df[, c("time","num",
              "s.num","i.num","si.flow",
              "acute.num","acute.flow",
              "chronic1.num","chronic1.flow",
              "chronic2.num","chronic2.flow",
              "final.num","final.flow",
              "a.flow",
              "d.flow")]
df3

#Data frame for detailed HIV status and ART status sub-compartment counts
df4 <- df[, c("time","num",
              "s.num","i.num","si.flow",
              "acute.ART.net.flow", "chronic1.ART.net.flow",
              "chronic2.ART.net.flow", "final.ART.net.flow",
              "ART.net.flow",
              "acute.num","acute.ART.num", "acute.NoART.num",
              "acute.flow",
              "chronic1.num","chronic1.ART.num","chronic1.NoART.num",
              "chronic1.flow","chronic1.ART.flow","chronic1.NoART.flow",
              "chronic2.num","chronic2.ART.num","chronic2.NoART.num",
              "chronic2.flow","chronic2.ART.flow","chronic2.NoART.flow",
              "final.num","final.ART.num","final.NoART.num",
              "final.flow","final.ART.flow","final.NoART.flow",
              "a.flow",
              "d.flow",
              "depart.standard.flow", "depart.standard.ART.flow",
              "depart.standard.NoART.flow",
              "depart.final.flow", "depart.final.ART.flow",
              "depart.final.NoART.flow")]
df4



## Plot outcomes

#SI Compartment Counts
par(mar = c(2,2,1,1), mgp = c(2,1,0))
plot(sim, y = c("s.num","i.num"),
     mean.col = 1:2, mean.lwd = 1, mean.smooth = TRUE,
     qnts.col = 1:2, qnts.alpha = 0.25, qnts.smooth = TRUE,
     legend = TRUE)

#HIV status sub-compartment counts
par(mar = c(5,5,1,1), mgp = c(2,1,0))
plot(sim, y = c("s.num","acute.num",
                "chronic1.num","chronic2.num",
                "final.num"),
     mean.col = 1:5, mean.lwd = 1, mean.smooth = TRUE,
     qnts.col = 1:5, qnts.alpha = 0.25, qnts.smooth = TRUE,
     legend = TRUE)


# Standardized Incidence and Prevalence
sim <- mutate_epi(sim, ir.rate = si.flow / s.num,
                  prev = i.num / num)
plot(sim, y = c("ir.rate", "prev"), legend = TRUE)
