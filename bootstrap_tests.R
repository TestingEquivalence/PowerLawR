bootstrap_stdev1<-function(p,n,nSimulation,kmin,kmax, tol){
  
  i=c(1:nSimulation)
  f<-function(k){
    v=rmultinom(n=1,size=n,prob=p)
    v=v/sum(v)
    cdf=cumsum(v)
    res = nearestPowerLaw(cdf,kmin,kmax,1,3,tol=tol)

    distance=res$objective
    return(distance*distance)
  }
  
  sample=sapply(i,f)
  return(sqrt(var(sample)))
}

bootstrap_test1<-function(alpha, frequency, kmin, kmax,
                          scale,nSimulation, tol=NA)
{
  #calcualte cdf
  n=sum(frequency)
  p=frequency/n
  cdf=cumsum(p)
  kmin=kmin/scale
  kmax=kmax/scale
  
  res = nearestPowerLaw(cdf,kmin,kmax,1,3,tol)
  beta=res$minimum
  distance=res$objective
  
  vol=bootstrap_stdev2(p,n,nSimulation,kmin,kmax,tol)
  qt=qnorm(1-alpha,0,1)
  min_eps = distance*distance + qt*vol
  min_eps=sqrt(min_eps)
  
  vec=c(min_eps,distance,beta,n)
  names(vec)=c("min_eps","distance","beta","sample_size")
  return(vec)
}

bootstrap_stdev2<-function(p,n,nSimulation,kmin,kmax, tol){
  
  i=c(1:nSimulation)
  f<-function(k){
    v=rmultinom(n=1,size=n,prob=p)
    v=v/sum(v)
    cdf=cumsum(v)
    res = nearestPowerLaw(cdf,kmin,kmax,1,3,tol=tol)
    
    #variance denominator
    #variance denominator
    pLawCDF=powerLawCDF(beta,kmin,kmax)
    drv=derivative(cdf,pLawCDF) 
    vol=asympt_stdev(frequency,drv)/sqrt(n)
    
    return(distance*distance/vol)
  }
  
  sample=sapply(i,f)
  return(sqrt(var(sample)))
}

bootstrap_test2<-function(alpha, frequency, kmin, kmax,
                         scale,nSimulation, tol=NA)
{
  #calcualte cdf
  n=sum(frequency)
  p=frequency/n
  cdf=cumsum(p)
  kmin=kmin/scale
  kmax=kmax/scale
  
  res = nearestPowerLaw(cdf,kmin,kmax,1,3,tol)
  beta=res$minimum
  
  #variance denominator
  pLawCDF=powerLawCDF(beta,kmin,kmax)
  drv=derivative(cdf,pLawCDF) 
  vol=asympt_stdev(frequency,drv)/sqrt(n)
  
  distance=res$objective/vol
  
  vol=bootstrap_stdev2(p,n,nSimulation,kmin,kmax,tol)
  qt=qnorm(1-alpha,0,1)
  min_eps = distance*distance + qt*vol
  min_eps=sqrt(min_eps)
  
  vec=c(min_eps,distance,beta,n)
  names(vec)=c("min_eps","distance","beta","sample_size")
  return(vec)
}
