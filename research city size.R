source("PowerLaw.R")
source("tests.R")
source("power.R")
source("read_data.R")
source("simulation.R")
source("size.R")
source("bootstrap_tests.R")

#vector of city sizes in Germany
citySize=readVector("C:\\Users\\Ostrovski\\Google Drive\\Writing\\PowerLaw\\CitySize\\list_ge.csv")

#different k_min values
kmins=c(10e3,20e3,30e3,40e3,50e3)

#different k_max values, first value coincides with largest city
kmaxs=c(max(citySize),5e6,10e6,20e6)

#scale for the populatuion measurement
#it is necessary for the computational feasibility.
scale=1000
#significance level
alpha=0.05

#carry out multiple asymptotic tests for the power law
#an asymptotic test is computed for each combination of k_min and k_max 
#use all but one CPU cores

result=multiple_test(alpha,citySize,kmins,kmaxs, scale, bootstrap = TRUE,
                     nSimulation = 1000)

#write test results
write.table(result$beta,paste("beta_",scale,".csv"))
write.table(result$distance,paste("distance_",scale,".csv"))
write.table(result$min_eps,paste("min_eps_",scale,".csv"))
write.table(result$sample_size,paste("sample_size_",scale,".csv"))

#MLE Estimator of beta
result=multiple_MLE(alpha,citySize,kmins,kmaxs,scale)

write.table(result$beta,paste("MLE_beta_",scale,".csv"))
write.table(result$sample_size,paste("MLE_sample_size_",scale,".csv"))

#compute test power at the power law points
###########################################
kmin=20e3
kmax=10e6
scale=10e3
nSamples=1000
n=662


# asymptotic test
test<-function(counting,kmin,kmax, scale){
  asymptotic_test(alpha=0.05,frequency = counting,kmin,kmax,scale, tol=0.001)
}

for (beta in c(2.1, 2.2, 2.3, 2.4, 2.5)) {
  size=sizeAtPowerLaw(n,kmin,kmax,scale,beta,nSamples)
  write.table(size, paste("size",beta,".csv"))
}


# compute test power at boundary points
###########################################

# asymptotic test 

pw=boundaryPower(n,eps,kmin,kmax,scale,beta,alpha, boundaryPointType = 1,
                 bootstrap = TRUE, nSimulation = 1000,tol=0.001, 
                 adjEps=adjEps)
write.table(pw, "powerLawStress.csv")

# pw=boundaryPower(n,eps,kmin,kmax,scale,beta,alpha, boundaryPointType = 2,
#                  bootstrap = FALSE, nSimulation = 1000,tol=0.001, 
#                  adjEps=adjEps)
# write.table(pw, "uniformRandomStress.csv")
 
#MLE at the power law
####################################
kmin=20e3
kmax=10e6
scale=10e3
nSamples=1000
n=662

test<-function(counting,kmin,kmax, scale){
  res=powerLawMLE(counting,kmin/scale,kmax/scale,1,3)
  return(res$minimum)
}

for (beta in c(2.1, 2.2, 2.3, 2.4, 2.5)) {
  size=sizeAtPowerLaw(n,kmin,kmax,scale,beta,nSamples)
  write.table(size, paste("size",beta,".csv"))
}




