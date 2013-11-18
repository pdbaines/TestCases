# load library
library(LNSSimple)

useC <- FALSE # TRUE
verbose <- TRUE
useMulti <- FALSE
maxit <- 30

# true values
set.seed(1234)

n <- 200
true.taus <- c(1e-17,5e-17)
true.betas <- c(0.5,3)
A <- rep(1e19,n)
b <- rep(10,n)

# simulate data
if (verbose){
	cat("Simulating data...\n")
}
dat <- sim.dat.pareto.pois.b(n,true.taus,true.betas,A,b)
if (verbose){
	cat("done. Fitting model...\n")
}

# fit the broken power law model with the knowledge of B
# ~680 seconds on my aging MBP with useC=TRUE
# Takes a lot longer with useC=FALSE (pure R code).

f1time <- system.time({
	fit1 <- broken.power.B(Y=dat$Y, A=A, B=2, b=b, Nlim=maxit, NM.maxit=150, NM.funevals=300, display=verbose, useC=useC)
})

if (verbose){
	cat(paste0("Finished fit 1, took ",round(f1time["elapsed"],3)," seconds...\n"))
}

# fit the broken power law model with automatic selection of B (take quite a bit of time)
# This is slow -- to make quicker, reduce maxit to a smaller number...
f2time <- system.time({
	fit2 <- broken.power(Y=dat$Y, A=A, b=b, Nlim=maxit, NM.maxit=150, NM.funevals=300, display=verbose, PPmulti=useMulti, useC=useC)
})

if (verbose){
	cat(paste0("Finished fit 2, took ",round(f2time["elapsed"],3)," seconds...\n"))
}

save.image(paste0("results_C_",useC,".RData"))






