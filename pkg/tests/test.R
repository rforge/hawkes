library(hawkes)

# Simple test that the log-likelihood function is working.
l <- likelihoodHawkes(1, 2, 3, 1:10)
show(l)