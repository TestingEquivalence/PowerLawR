source("PowerLaw.R")
source("tests.R")
source("power.R")
source("read_data.R")
source("simulation.R")
source("size.R")
source("bootstrap_tests.R")

#vector of word counts in contemporary american english corpora
words=readVector("C:\\data\\words.csv")

parameter=list(scale=1000,alpha=0.05, tol=0.001,nSimulation=1000,  
               test=asymptotic, 
               kmins=c(1000,2000,3000, 5000, 10000), 
               kmaxs=c(max(words),2e6,5e6,10e6), 
               counting=words)

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
parameter=list(kmin=1000, kmax=5e6, scale=1000, nSamples=1000, n=1984,
               alpha=0.05, tol=0.001, nSimulation=100, 
               test=asymptotic)  


for (beta in c(1.9, 2.0, 2.05, 2.1)) {
  parameter$beta=beta
  size=sizeAtPowerLaw(parameter)
  write.table(size, paste("size",beta*100,".csv"))
}
 
# compute test power at boundary points
###########################################
parameter=list(kmin=700, kmax=5e6, scale=100, nSamples=1000, n=2839, 
               beta=2.05, alpha=0.05, tol=0.001, nSimulation=100,
               test=asymptotic, epsAdj=1)  


for (eps in c(0.20)){
  parameter$eps=eps
  pw=boundaryPower(parameter)
  write.table(pw, paste("uniformRandomStress",eps*100,".csv"))
}
