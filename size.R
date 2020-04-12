fullToss<-function(i,p, n, kmin,kmax, scale){
  set.seed(i*1000000)
  counting=rmultinom(n=1,size=n,prob=p)
  res=test(counting,kmin,kmax, scale)
  return(res)
}

sizeAtPowerLaw<-function(n,kmin,kmax, scale, beta,nSamples){
  #calculate density of discrete power law
  p=powerLawDensity(beta,kmin=kmin/scale,kmax=kmax/scale)
  
  #simulate tests
  #v=sapply(i, fullToss,p,n,kmin,kmax,scale)
  i=c(1:nSamples)
  cl=getCluster()
  clusterExport(cl,c("test"))
  v=parSapply(cl,i, fullToss,p,n,kmin,kmax,scale)
  stopCluster(cl)
  
  return(v)
}

