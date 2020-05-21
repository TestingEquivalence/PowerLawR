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

parameter=list(scale=1e3,alpha=0.05, tol=0.001,nSimulation=1000, nDirections = 100, 
               test="bootstrap2", kmins=kmins, kmaxs=kmaxs, counting=citySize)

#carry out multiple tests for the power law
#given test is computed for each combination of k_min and k_max 
#use all but one CPU cores

result=multiple_test(parameter)

#write test results
write.table(result$beta,paste("beta.csv"))
write.table(result$distance,paste("distance.csv"))
write.table(result$min_eps,paste("min_eps.csv"))
write.table(result$sample_size,paste("sample_size.csv"))

#MLE Estimator of beta
result=multiple_MLE(parameter)

write.table(result$beta,paste("MLE_beta.csv"))
write.table(result$sample_size,paste("MLE_sample_size.csv"))

#compute test power at the power law points
###########################################
parameter=list(kmin=20e3,kmax=10e6,scale=1e3,nSamples=1000,n=662,
          alpha=0.05, tol=0.001,nSimulation=1000, nDirections = 100, test="asymptotic",
          eps=0.08)  


for (beta in c(2.1,2.2,2.3,2.4,2.5)) {
  parameter$beta=beta
  size=sizeAtPowerLaw(parameter)
  write.table(size, paste("size",beta*100,".csv"))
}


# compute test power at boundary points
###########################################
parameter=list(kmin=20e3,kmax=10e6,scale=1e3,nSamples=1000,n=662, beta=2.3,
               alpha=0.05, tol=0.001,nSimulation=1000, nDirections = 100,
               test="MLE",eps=0.12,epsAdj=1)  


for (eps in c(0.08,0.10,0.12)){
  parameter$boundaryPointType=1
  parameter$eps=eps
  pw=boundaryPower(parameter)
  write.table(pw, paste("powerLawStress",eps*100,".csv"))
}

for (eps in c(0.08,0.10,0.12)){
  parameter$boundaryPointType=2
  parameter$eps=eps
  pw=boundaryPower(parameter)
  write.table(pw, paste("uniformRandomStress",eps*100,".csv"))
}
 




