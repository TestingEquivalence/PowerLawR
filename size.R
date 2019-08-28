testPowerAtPowerLaw<-function(counting,kmin,kmax,scale, beta, eps, nSamples, alpha){
  #compute counting frequnces
  df=list2freq(counting,kmin,kmax,scale)
  
  #calculate sample size
  n=sum(df)
  skmin=kmin/scale
  skmax=kmax/scale
  
  #calculate density of discrete power law
  p=powerLawDensity(beta,kmin,kmax)
  powerAtPoint(p,n,nSamples,skmin,skmax,scale=1,eps,alpha)
}