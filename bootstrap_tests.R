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
  
  vol=bootstrap_stdev1(p,n,nSimulation,kmin,kmax,tol)
  qt=qnorm(1-alpha,0,1)
  min_eps = distance*distance + qt*vol
  min_eps=sqrt(min_eps)
  
  vec=c(min_eps,distance,beta,n)
  names(vec)=c("min_eps","distance","beta","sample_size")
  return(vec)
}

closeRandomPoint<-function(p,n, eps,beta,kmin,kmax){
  repeat{
    v=rmultinom(n=1,size=n,prob=p)
    v=v/n
    cv=cumsum(v)
    res=nearestPowerLaw(cv,kmin,kmax,1,3)
    if (res$objective>eps)
      if (res$minimum>beta*0.8)
        if (res$minimum<beta*1.2)
          return(v)
  }
}

linComb<-function(x,y,a){
  return((1-a)*x+a*y) 
}

linearBoundaryPoint<-function(p,q,eps,kmin,kmax){
  P=cumsum(p)
  Q=cumsum(q)
  
  aim<-function(a){
    lc=linComb(P,Q,a)
    res = nearestPowerLaw(lc,kmin,kmax,1,3)
    beta=res$minimum
    distance=res$objective
    return(distance-eps)
  }
  
  aMin=uniroot(aim, c(0,1))
  return(linComb(p,q,aMin$root))
}


bootstrap_test2<-function(alpha, frequency, kmin, kmax,
                         scale,nSimulation, tol=NA, nDirections)
{
  #calcualte cdf
  n=sum(frequency)
  p=frequency/n
  cdf=cumsum(p)
  kmin=kmin/scale
  kmax=kmax/scale
  
  res = nearestPowerLaw(cdf,kmin,kmax,1,3,tol)
  beta=res$minimum
  
  #generate exterior points
  set.seed(03042020)
  extPoints=list()
  for (i in c(1:nDirections)){
    ep=closeRandomPoint(p,n, eps,beta,kmin,kmax)
    extPoints=c(extPoints,ep)
  }
  
  #generate linear boundary points
  bndPoints=lapply(extPoints, linearBoundaryPoint,q=p,eps,kmin,kmax)
  
  #find closes bnd point
  cbndPoints=lapply(bndPoints, cumsum)
  cp=cumsum(p)
  dst=lapply(bndPoints, l2, F2=cp)
  pos=which.min(dst)
  bndPoint=bndPoints[[pos]]
  
  vec=c(min_eps,distance,beta,n)
  names(vec)=c("min_eps","distance","beta","sample_size")
  return(vec)
}
