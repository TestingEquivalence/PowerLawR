
derivative<-function(F1,F2)
{
  Z=2*(F1-F2)
  Z=rev(Z)
  Z=cumsum(Z)
  Z=rev(Z)
  return(Z)
}

asympt_stdev<-function(p,derivative){
  vec = derivative
  vnsq_1  = sum(p*vec*vec)
  
  k=length(p)
  vnsq_2=0
  for (j1 in 1:k)
    for (j2 in 1:k)
      vnsq_2 = vnsq_2 + vec[j1] * vec[j2] * p[j1] * p[j2]
  
  
  vnsq  = vnsq_1 - vnsq_2
  return (sqrt(vnsq))
}

asymptotic_test<-function(alpha, counting, kmin, kmax, scale)
{
  #calcualte counting frequnce
  df=list2freq(counting,kmin,kmax,scale)
  n=sum(df)
  df=df/n
  cdf=cumsum(df)
  kmin=kmin/scale
  kmax=kmax/scale
  
  res = nearestPowerLaw(cdf,kmin,kmax,1,3)
  beta=res$minimum
  min_distance=res$objective
  pLawCDF=powerLawCDF(beta,kmin,kmax)
  
  drv=derivative(cdf,pLawCDF)
  vol=asympt_stdev(df,drv)
  vol=vol/ sqrt(n);
  
  qt=qnorm(1-alpha,0,1)
  min_eps = min_distance*min_distance + qt*vol
  min_eps=sqrt(min_eps)
  
  vec=c(min_eps,min_distance,beta)
  names(vec)=c("min_eps","min_distance","beta")
  return(vec)
}