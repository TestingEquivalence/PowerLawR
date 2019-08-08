sizeAtPowerLaw<-function(counting,kmin,kmax,scale, beta, eps, nSamples){
  #calcualte counting frequnce
  df=list2freq(counting,kmin,kmax,scale)
  n=sum(df)
  df=df/n
  cdf=cumsum(df)
  kmin=kmin/scale
  kmax=kmax/scale
  
  #calculate density of discrete power law
  p=powerLawDensity(beta,kmin,kmax)
  powerAtPoint(p,n,nSamples,kmin,kmax,scale=1,eps)
}