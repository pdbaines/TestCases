## Bayesian Logistic Regression

Simple MCMC code to illustrate Bayesian Logistic Regression, with 
a multivariate normal prior for the regression parameter.
Two algorithms are available for sampling, a Metropolis
algorithm with a multivariate-normal random walk proposal
(and automated scale tuning for a fixed correlation matrix),
and a Metropolis-within-Gibbs algorithm (again, with
automated tuning for each individual parameter).

Application is provided for the breast cancer dataset
from the UCI machine learning repository.

Requires:

    library(mvtnorm)
    library(coda)
    library(MASS)

Run using:

    source("BLR_data.R")

