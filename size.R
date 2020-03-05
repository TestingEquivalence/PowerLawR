

sizeAtPowerLaw<-function(n,kmin,kmax,scale, beta,  nSamples, 
                         alpha, bootstrap=FALSE, nSimulation=0, tol=NA){
  kmin=kmin/scale
  kmax=kmax/scale
  
  #calculate density of discrete power law
  p=powerLawDensity(beta,kmin,kmax)
  
  #simulate tests
  i=c(1:nSamples)
  cl=getCluster()
  v=parSapply(cl,i, fullToss,p,n,kmin,kmax,scale=1,alpha, bootstrap, nSimulation,tol)
  #v=sapply(i, fullToss,p,n,kmin,kmax,scale=1,alpha, bootstrap, nSimulation,tol)
  
  stopCluster(cl)
  return(v)
}
MLEatPowerLaw<-function(n,kmin,kmax,scale, beta,  nSamples, 
                         alpha, bootstrap=FALSE, nSimulation=0, tol=NA){
  kmin=kmin/scale
  kmax=kmax/scale
  
  #calculate density of discrete power law
  p=powerLawDensity(beta,kmin,kmax)
  
  #simulate tests
  i=c(1:nSamples)
  cl=getCluster()
  v=parSapply(cl,i, MLEToss,p,n,kmin,kmax)
  stopCluster(cl)
  return(v)
}
