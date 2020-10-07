source("PowerLaw.R")
source("tests.R")
source("power.R")
source("read_data.R")
source("simulation.R")
source("size.R")
source("bootstrap_tests.R")

# generate syntetic power law data first
#-----------------------------------------
kmin=1 # lower cut off >=1
kmax=10000 # upper cut off
beta=2 #power exponent
massFunc=powerLawDensity(beta,kmin,kmax) #mass function of discrete power law


n=2000 # number of samples
set.seed(1012021)
simulation=rmultinom(n=1,size=n,prob=massFunc) #simulate from power law

# apply asymptotic test to simulated frequency data
#-----------------------------------------


# Attention, kmax equals length(frequency)+kmin-1
result=asymptotic_test(alpha=0.05,
                       frequency = simulation,
                       kmin=kmin)
print(result)

# apply bootstrap test to simulated frequency data
#-----------------------------------------

