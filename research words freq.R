source("PowerLaw.R")
source("tests.R")
source("power.R")
source("read_data.R")
source("simulation.R")
source("size.R")
source("bootstrap_tests.R")

#vector of city sizes in Germany
words=readVector("C:\\Users\\Ostrovski\\Google Drive\\Writing\\PowerLaw\\WordCounts\\words.csv")

#different k_min values
kmins=c(500,600,700,800,1000)
#kmins=c(1000,2000,5000, 10000)

#different k_max values, first value coincides with most frequent word
kmaxs=c(max(words),2e6,5e6,10e6)

#scale for the words frequency
#it is necessary for the computational feasibility.
scale=10
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
  asymptotic_test(alpha=0.05,frequency = counting,kmin,kmax,scale, 
                  tol=0.001)
}

#bootstrap test 1
test<-function(counting,kmin,kmax,scale,eps){
  bootstrap_test1(alpha=0.05,frequency = counting,
                  kmin,kmax,scale,nSimulation=1000,tol=0.001)
}

#bootstrap test 2
test<-function(counting,kmin,kmax,scale,eps){
  bootstrap_test2(frequency =counting ,kmin,kmax,scale,
                  nSimulation=1000,nDirections = 100 ,eps,tol=0.001)
}

for (beta in c(2.1, 2.2, 2.3, 2.4, 2.5)) {
  size=sizeAtPowerLaw(n,kmin,kmax,scale,beta,nSamples)
  write.table(size, paste("size",beta,".csv"))
}


# compute test power at boundary points
###########################################
kmin=20e3
kmax=10e6
beta=2.3
scale=10e3
nSamples=1000
n=662
epsAdj=1

# asymptotic test 

test<-function(counting,kmin,kmax,scale,eps){
  res=asymptotic_test(alpha=0.05,frequency = counting,kmin,kmax,scale, 
                      tol=0.001)
  return(res[1]<eps)
}

#bootstrap test 1
test<-function(counting,kmin,kmax,scale,eps){
  res=bootstrap_test1(alpha =0.05,frequency = counting,
                      kmin,kmax,scale,nSimulation=1000,tol=0.001)
  return(res[1]<eps)
}

#bootstrap test 2
test<-function(counting,kmin,kmax,scale,eps){
  pval=bootstrap_test2(frequency =counting ,kmin,kmax,scale,
                       nSimulation=1000,nDirections = 100 ,eps,tol=0.001)
  return(pval<0.05)
}

for (eps in c(0.08,0.10,0.12)){
  pw=boundaryPower(n,eps,kmin,kmax,scale,beta,nSamples, boundaryPointType=1,epsAdj)
  write.table(pw, paste("powerLawStress",eps*100,".csv"))
}

for (eps in c(0.08,0.10,0.12)){
  pw=boundaryPower(n,eps,kmin,kmax,scale,beta,nSamples, boundaryPointType=2,epsAdj)
  write.table(pw, paste("uniformRandomStress",eps*100,".csv"))
}

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





