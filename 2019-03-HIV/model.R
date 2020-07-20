
##
## HIV Transmission Model
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##
## Authors: Connor Van Meter, Samuel M. Jenness, Yuan Zhao, Emeli Anderson
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
nw <- network_initialize(n)

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
dx <- netdx(est, nsims = nsims, ncores = ncores, nsteps = nsteps)

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
                   ART.Treatment.Rate = 0.01,
                   ART.Discontinuance.Rate = 0.005,
                   ART.Progression.Reduction.Rate = 0.5,
                   arrival.rate = 0.002,
                   departure.rate = departure_rate)

# Initial conditions
start_prevalence = 0.05
init <- init.net(i.num = round(start_prevalence * n))


# Read in the module functions
if (interactive()) {
  source("2019-03-HIV/module-fx.R", echo = TRUE)
} else {
  source("module-fx.R")
}

# Control settings
control <- control.net(type = NULL,
                       nsteps = nsteps,
                       nsims = nsims,
                       ncores = ncores,
                       infection.FUN = infect,
                       progress.FUN = progress,
                       arrivals.FUN = afunc,
                       departures.FUN = dfunc,
                       resimulate.network = TRUE,
                       verbose = TRUE,
                       module.order = c("resim_nets.FUN", "progress.FUN",
                                        "infection.FUN", "arrivals.FUN",
                                        "departures.FUN", "prevalence.FUN"))

# Run the network model simulation with netsim
sim <- netsim(est, param, init, control)
print(sim)


# Examine the data from the simulation
df <- as.data.frame(sim)
df

df2 <- df[, c("num",
              "s.num", "i.num", "acute.ART.num", "acute.NoART.num", "chronic1.ART.num",
              "chronic1.NoART.num", "chronic2.ART.num", "chronic2.NoART.num",
              "AIDS.ART.num", "AIDS.NoART.num")]
df2
df2[c(1, 2, 0.2*nsteps, 0.4*nsteps, 0.6*nsteps, 0.8*nsteps, nsteps), ]


## Plot outcomes

# SI Compartment Counts
par(mar = c(3,3,3,1), mgp = c(2,1,0))
plot(sim, y = c("s.num","i.num"),
     mean.col = 1:2, mean.lwd = 1, mean.smooth = TRUE,
     qnts.col = 1:2, qnts.alpha = 0.25, qnts.smooth = TRUE,
     legend = TRUE)

# HIV and ART status compartment counts
par(mar = c(3,3,3,1), mgp = c(2,1,0))
plot(sim, y = c("s.num", "acute.ART.num", "acute.NoART.num", "chronic1.ART.num",
                "chronic1.NoART.num", "chronic2.ART.num", "chronic2.NoART.num",
                "AIDS.ART.num", "AIDS.NoART.num"),
     mean.col = 1:9, mean.lwd = 1, mean.smooth = TRUE)
legend("topright", legend = c("s.num", "acute.ART.num", "acute.NoART.num",
                              "chronic1.ART.num", "chronic1.NoART.num",
                              "chronic2.ART.num", "chronic2.NoART.num",
                              "AIDS.ART.num", "AIDS.NoART.num"), lty = 1,
       cex = 0.5, col = c(1:9))

# HIV and ART status compartment counts w/out susceptible compartment
par(mar = c(3,3,3,1), mgp = c(2,1,0))
plot(sim, y = c("acute.ART.num", "acute.NoART.num", "chronic1.ART.num",
                "chronic1.NoART.num", "chronic2.ART.num", "chronic2.NoART.num",
                "AIDS.ART.num", "AIDS.NoART.num"),
     mean.col = 1:8, mean.lwd = 1, mean.smooth = TRUE)
legend("topleft", legend = c("acute.ART.num", "acute.NoART.num",
                             "chronic1.ART.num", "chronic1.NoART.num",
                             "chronic2.ART.num", "chronic2.NoART.num",
                             "AIDS.ART.num", "AIDS.NoART.num"), lty = 1,
       cex = 0.5, col = c(1:8))

# Standardized Incidence and Prevalence
sim <- mutate_epi(sim, ir.rate = acute.flow / s.num,
                  prev = i.num / num)
par(mar = c(3,3,3,1), mgp = c(2,1,0), mfrow = c(1,2))
plot(sim, y = "prev", main = "Prevalence", ylab = "")
plot(sim, y = "ir.rate", main = "Incidence", ylab = "")

# ART Treatment Prevalence
sim <- mutate_epi(sim, ART.num =
                    acute.ART.num + chronic1.ART.num +
                    chronic2.ART.num + AIDS.ART.num)
sim <- mutate_epi(sim, ART.prev = ART.num / i.num)
par(mfrow = c(1,1))
plot(sim, y = "ART.prev", main = "ART Treatment Prevalence", ylab = "")
