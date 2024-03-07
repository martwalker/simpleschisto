#########################################
## main driver file
########################################

# clear workspace
rm(list=ls())

# load packages
source("packages.R")

# load utility functions
source("funcs.R")

# load default parameters
source("parms.R")

## set transmission conditions & calculate initial value
parameters["R0"] <- 2.8
parameters["k"] <- .3
parameters["N1N2"] <- .01
W0 <- funcs$findstable3(par=parameters)

## set MDA parameters
parameters["stop.t"] <- 100
parameters["coverage"] <- 0.6
parameters["start.tx"] <- 50
parameters["n.tx"] <- 14
parameters["dotx"] <-1
out <- funcs$runmod(parameters, W0=W0)
plot(out$W~out$time, type="l")
plot(out$E~out$time, type="l")
plot(out$prev~out$time, type="l")
plot(out$kdyn~out$time, type="l")
plot(out$RE~out$time, type="l")
plot(out$y~out$time, type="l")
plot(out$prevh~out$time, type="l")


# create LHS
#source("cluster/createLHS.R")
#source("createLHS.R")

# run stability analysis
#source("cluster/calcequib.R")
#source("calcequib.R")

# process stability analysis output
#source("cluster/processequib.R")
#source("processequib.R")

# run treatment scenario 
#source("cluster/runtx.R")
#source("runtx.R")

# process dynamic model output
#source("cluster/processdyn.R")
#source("processdyn.R")

# clear the stability analysis (not used)
#source("cleanstability.R")

# join stability analysis (static) with dynamic data
source("joinstabilitydyn.R")

# calculate PRCCs
source("PRCC.R")


# load dynamic model output
#dyndat <- readRDS("dyndat.rds")

# load times to elimination output
#timedat <- readRDS("timedat.rds")

# make plots
#source("graphpar.R")




