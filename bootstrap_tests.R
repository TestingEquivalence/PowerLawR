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


bootstrap_test_base<-function(alpha, p, kmin, kmax,
                              nSimulation, tol, eps, exteriorPoints)
{
  cdf=cumsum(p)
  res = nearestPowerLaw(cdf,kmin,kmax,1,3,tol)
  beta=res$minimum
  distance=res$objective
  
  if (distance>=eps) {
    p_value=1
    return(p_value)
  }
  
  #generate boundary points
  bndPoints=list()
  nDir=length(exteriorPoints)
  for (i in c(1:nDir)){
    bndPoints[[i]]=linearBoundaryPoint(p,q=exteriorPoints[[i]],eps,kmin,kmax)
  }
  
  #find closes bnd point
  cbndPoints=lapply(bndPoints, cumsum)
  cp=cumsum(p)
  dst=lapply(cbndPoints, l2, F2=cp)
  pos=which.min(dst)
  bndPoint=bndPoints[[pos]]
  
  res = nearestPowerLaw(cumsum(bndPoint),kmin,kmax,1,3,tol)
  
  
  #simulate bootstrap distribution
  i=c(1:nSimulation)
  
  f<-function(k){ 
    v=rmultinom(n=1,size=n,prob=bndPoint)
    v=v/sum(v)
    cdf=cumsum(v)
    res = nearestPowerLaw(cdf,kmin,kmax,1,3,tol=tol)
    
    return(res$objective)
  }
  
  sample=sapply(i,f)
  p_value=sum(distance>=sample)/nSimulation
  return(p_value)
}

bootstrap_test2<-function(alpha, frequency, kmin, kmax, scale,
                          nSimulation, nDirections, eps, tol){
  n=sum(frequency)
  p=frequency/n
  kmin=kmin/scale
  kmax=kmax/scale
  
  #generate H0 points
  exteriorPoints=list()
  for (i in c(1:nDirections)){
    exteriorPoints[[i]]=closeRandomPoint(p,n, eps,beta,kmin,kmax)
  }
  
  res=bootstrap_test_base(alpha,p,kmin,kmax,nSimulation,tol,eps,
                          exteriorPoints)
  return(res)
}
  

bootstrap_test3<-function(alpha, frequency, kmin, kmax,
                          scale,nSimulation, tol=NA, nDirections,minEps, maxEps){
  n=sum(frequency)
  p=frequency/n
  kmin=kmin/scale
  kmax=kmax/scale
  
  #generate H0 points
  exteriorPoints=list()
  for (i in c(1:nDirections)){
    exteriorPoints[[i]]=closeRandomPoint(p,n, maxEps,beta,kmin,kmax)
  }
  
  f<-function(eps){
    p_val=bootstrap_test_base(alpha,p,kmin,kmax,nSimulation,tol,eps,
                              exteriorPoints)
    return(p_val-alpha)
  }
  
  res=uniroot(f,lower=minEps,upper=maxEps)
  return(res$root)
}
