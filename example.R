source("PowerLaw.R")
source("tests.R")
source("power.R")
source("read_data.R")
source("simulation.R")
source("size.R")
source("bootstrap_tests.R")

# generate synthetic power law data first
#-----------------------------------------
kmin=1 # lower cut off >=1
kmax=10000 # upper cut off
beta=2 #power exponent
massFunc=powerLawDensity(beta,kmin,kmax) #mass function of discrete power law


n=2000 # number of observations
set.seed(1012021)
simulation=rmultinom(n=1,size=n,prob=massFunc) #simulate from power law

# apply asymptotic test to simulated frequency data
#-----------------------------------------


# Attention, kmax will be calculated from kmin and length of the frequency
# If you need larger kmax parameter then simply add sufficiently number of zeros
# to the right side of the frequency vector.

# The parameter of the asymptotic tests are:
# alpha - significance level 
# frequency - vector of observed frequencies
# kmin - lower bound of power law 
result=asymptotic_test(alpha=0.05,
                       frequency = simulation,
                       kmin=kmin)

# the results contains:
# min_eps - the minimum tolerance parameter for which H0 can be rejected
# distance -Euclidean distance between empirical CDF and CDF of the closest power law
# beta - minimum distance estimate of power law exponent
# sample_size - sample size

print(result)

# apply bootstrap test to simulated frequency data
#-----------------------------------------

# Attention, kmax will be calculated from kmin and length of the frequency
# If you need larger kmax parameter then simply add sufficiently number of zeros
# to the right side of the frequency vector.

# The parameter of the asymptotic tests are:
# alpha -  significance level
# frequency - vector of observed frequencies
# kmin - lower bound of power law
# nSimulation - number of bootstrap samples

result=bootstrap_test(alpha=0.05,
                      frequency = simulation,
                      nSimulation = 1000,
                      kmin=kmin)

# the results contains:
# min_eps - the minimum tolerance parameter for which H0 can be rejected
# distance -Euclidean distance between empirical CDF and CDF of the closest power law
# beta - minimum distance estimate of power law exponent
# sample_size - sample size

print(result)

#------------------------------------------------------------------------------------
# Usually we do not observe the frequency data directly. 
# Instead we observe the counting data only and also do not know the upper and lower bound
# of the power low. In this case we need to transform the counting date to frequencies and
# we also may perform multiple equivalence tests for different values of upper 
# and lower bounds. Here we show the convenient function "multiple_test" to perform these
# tasks efficiently.

# First we transform our frequency data t counting data:

counting=c()
for (i in c(kmin:kmax)){
  m=simulation[i-kmin+1]
  counting=c(counting,rep(i,m))
}

# Next we define parameters for the multiple equivalence tests:

parameter=list(scale=1, #scaling maybe necessary to make computations feasible
               alpha=0.05, #significance level 
               nSimulation=1000, #number of bootstrap replications  
               test=asymptotic, # test should be asymptotic or bootstrap
               kmins=c(1,2,5,10), # possible lower bounds of power law
               kmaxs=c(5000,10000,20000), # possible upper bounds of power law
               counting=counting) # counting data  where the upper and lower bounds are unknown
result=multiple_test(parameter)

# Table of the estimated beta's (minimum distance estimate of power law exponent): 
result$beta

# Table of Euclidean distances between empirical CDF and CDF of the closest power law
result$distance

# Table of sample sizes
result$sample_size

# Table of  min_eps - the minimum tolerance parameters for which H0 can be rejected
result$min_eps

