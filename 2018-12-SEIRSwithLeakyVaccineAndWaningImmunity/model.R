##
## SEIR Model with Vital Dynamics and an All or Nothing Vaccine Implementation
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##
## Authors: Connor M. Van Meter
## Date: November 2018
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
mean_degree <- 3
edges <- mean_degree * (n/2)

# Input the appropriate target statistics for each term
target.stats <- c(edges)

#Set mortality rate
mr_rate = 0.0085/52

# Parameterize the dissolution model
coef.diss <- dissolution_coefs(dissolution = ~offset(edges),
                               duration = 20, d.rate = mr_rate)
coef.diss

# Fit the model
est <- netest(nw, formation, target.stats, coef.diss)

# Model diagnostics
dx <- netdx(est, nsims = 1, nsteps = 1040)

print(dx)
plot(dx)


# Epidemic model simulation -----------------------------------------------

# Model parameters
param <- param.net(inf.prob = 0.8, #Up to 80% secondary attack rate (CDC)
                   birth.rate = 0.0135/52, #US birth rate in 2016 was 13.5 per 1000 population
                   mortality.rate = mr_rate, #US mortality rate in 2016 was 850 per 100,000
                   mortality.disease.mult = 1.005, #0.5% mortality rate in infants under 6 months
                   act.rate = 1,
                   ei.rate = 0.875, #mean latent period of pertussis: 8 days
                   ir.rate = 0.4, #mean infectious period of pertussis: 14-21 days (2.5 weeks considered)
                   rs.rate = 0.0048, #infection-acquired immunity lasts between 4-20 years (4 years considered here - 208 weeks)
                   vaccination.rate.initialization = 0.172, #17.2% Tdap coverage in adults reported in 2013
                   protection.rate.initialization = 0.7, #Estimate Tdap protects about 7/10 people who receive it
                   vaccination.rate.progression.disease.experienced = 0.0005,
                   vaccination.rate.progression.disease.naive = 0.0001,
                   protection.rate.progression = 0.7,
                   vaccination.rate.births = 0.875, #DTaP coverage in children 19-35 months, 2017
                   protection.rate.births = 0.9, #9 out of 10 children fully protected after 5th DTaP dose
                   leaky.degree.of.protection.max = 0.872, #84% vaccine effectiveness at >8 years since last vaccination
                   leaky.degree.of.protection.min = 0.528, #41% vaccine effectiveness at >8 years since last vaccination
                   time.to.vaccine.protection.decay = 416, #Vaccine effectiveness decays from 0.84 to 0.41 over 8 years
                   decay.error.tolerance = 0.01 #Used for determining the exponential decay constant, lambda
)

# Initial conditions
init <- init.net(i.num = 20)

# Read in the module functions
source("module-fx.R", echo = TRUE)

# Control settings
control <- control.net(nsteps = 520,
                       nsims = 1,
                       ncores = 1,
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

#Review transmission probabilities for vaccine protected
# and non-vaccine protected individuals
x <- c(1:time.to.vaccine.protection.decay)
t = rep(0,time.to.vaccine.protection.decay)
lambda = -log(decay.error.tolerance)/time.to.vaccine.protection.decay
y = (1 - exp(-lambda*(x - t)))*((1 - leaky.degree.of.protection.min) - (1 - leaky.degree.of.protection.max)) + (1 - leaky.degree.of.protection.max)
plot(x, y, xlab = "Time in Weeks", ylab = "Transmission Probability", ylim = c(0,1.5), type = 'l', col = "blue")
lines(x, rep(inf.prob, time.to.vaccine.protection.decay), type = 'l', lty = 2, col = "red")
legend("topleft", legend = c("Vaccine-Protected Transmission Probability", "Non Vaccine-Protected Transmission Probability"), col = c("blue", "red"), lty = 1:2, cex = 0.8)

#Examine the data from the simulation
df <- as.data.frame(sim)
df

#Data frame for SEIR-V compartment counts
df2 <- df[, c("time", "num", "s.num", "e.num", "i.num", "r.num", "v.num", "b.num",
              "b.flow", "d.num", "d.flow")]
df2

#Data frame for vaccination flow and vaccination flow breakdown
#by vaccination method
df3 <- df[, c("vac.flow", "vac.init.flow", "vac.prog.flow", "vac.birth.flow",
              "vac.num", "vac.init.num", "vac.prog.num", "vac.birth.num")]
df3

#Data frame for vaccination protection flow and vaccination protection breakdown
#by vaccination method
df4 <- df[, c("prt.flow", "prt.init.flow", "prt.prog.flow", "prt.birth.flow",
              "prt.num", "prt.init.num", "prt.prog.num", "prt.birth.num")]
df4

#Epidemic plot of SEIR-V compartment counts, entrances, and exits over simulation
par(mar = c(3,3,1,1), mgp = c(2,1,0))
plot(sim, y = c("s.num","e.num","i.num","r.num", "v.num", "b.num", "d.num", "num"),
     mean.col = 1:8, mean.lwd = 1, mean.smooth = FALSE,
     qnts = 1, qnts.col = 1:8, qnts.alpha = 0.25, qnts.smooth = FALSE,
     legend = TRUE)

#Cumulative Incidence and Prevalence Plots of the Vaccine Model
sim <- mutate_epi(sim, ci = se.flow / s.num, prev = e.num / num)
plot(sim, y = c("ci", "prev"), mean.lwd = 1, mean.smooth = TRUE, legend = TRUE)


##################################################################
#VACCINE MODEL COMPARISON
##################################################################

# Update one or more vaccine parameters, simulate the new model,
# and compare with the original model

# Model parameters
param <- param.net(inf.prob = 0.5,
                   birth.rate = 0.01,
                   mortality.rate = mr_rate,
                   mortality.disease.mult = 2,
                   act.rate = 1,
                   ei.rate = 0.05,
                   ir.rate = 0.05,
                   vaccination.rate.initialization = 0.05,
                   protection.rate.initialization = 0.3,
                   vaccination.rate.progression = 0.05,
                   protection.rate.progression = 0.3,
                   vaccination.rate.births = 0.2,
                   protection.rate.births = 0.3
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
sim2 <- netsim(est, param, init, control)
print(sim2)

##Compare incidence and prevalence of simulation 1 to simulation 2

#Calculate cumulative incidence and prevalence of simulation 2
sim2 <- mutate_epi(sim2, ci2 = se.flow / s.num, prev2 = e.num / num)

par(mfrow = c(1,1))
plot(sim, y = c("ci", "prev"), mean.lwd = 1, mean.smooth = TRUE, legend = TRUE)
plot(sim2, y = c("ci2", "prev2"), mean.lwd = 1, mean.smooth = TRUE, add = TRUE,
     mean.col = c("steelblue", "firebrick"), qnts.col = c("steelblue", "firebrick"))
