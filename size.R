

sizeAtPowerLaw<-function(n,kmin,kmax,scale, beta,  nSamples, 
                         alpha, bootstrap=FALSE, nSimulation=0){
  kmin=kmin/scale
  kmax=kmax/scale
  
  #calculate density of discrete power law
  p=powerLawDensity(beta,kmin,kmax)
  
  #simulate tests
  i=c(1:nSamples)
  cl=getCluster()
  v=parSapply(cl,i, fullToss,p,n,kmin,kmax,scale=1,alpha)
  stopCluster(cl)
  return(v)
}
