

sizeAtPowerLaw<-function(n,kmin,kmax,scale, beta, eps, nSamples, alpha){
  kmin=kmin/scale
  kmax=kmax/scale
  
  #calculate density of discrete power law
  p=powerLawDensity(beta,kmin,kmax)
  
  #simulate tests
  i=c(1:nSamples)
  cl=getCluster()
  v=parSapply(i, fullToss,p=p,n=n,kmin=kmin,kmax=kmax,scale=1,eps=eps,alpha=alpha)
  stopCluster(cl)
  return(v)
}