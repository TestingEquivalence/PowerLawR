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


bootstrap_test2<-function(frequency, kmin, kmax, scale,
                          nSimulation, nDirections, eps, tol)
{
  #rescale
  n=sum(frequency)
  p=frequency/n
  kmin=kmin/scale
  kmax=kmax/scale
  
  #calculate distance
  cdf=cumsum(p)
  res = nearestPowerLaw(cdf,kmin,kmax,1,3,tol)
  beta=res$minimum
  distance=res$objective
  
  #check if it could work
  if (distance>=eps) {
    p_value=1
    return(p_value)
  }
  
  # generate H0 points
  exteriorPoints=list()
  for (i in c(1:nDirections)){
    exteriorPoints[[i]]=closeRandomPoint(p=p,n=n,eps = eps,beta=beta,kmin = kmin,kmax = kmax)
  }
  
  #generate boundary points
  bndPoints=list()
  nDir=length(exteriorPoints)
  for (i in c(1:nDir)){
    bndPoints[[i]]=linearBoundaryPoint(p=p,q=exteriorPoints[[i]],eps=eps,kmin=kmin,kmax=kmax)
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


bootstrap_test2_1<-function(alpha,frequency, kmin, kmax,
                          scale,nSimulation, tol=NA, 
                          nDirections,minEps=NA, maxEps=NA){
  n=sum(frequency)
  p=frequency/n
  kmin=kmin/scale
  kmax=kmax/scale
  
  #default min eps
  cdf=cumsum(p)
  res = nearestPowerLaw(cdf,kmin,kmax,1,3, tol)
  beta=res$minimum
  distance=res$objective
    
  if (is.na(minEps)){
    minEps=distance
  }
    
  if (is.na(maxEps)){
    maxEps=distance*1.3
  }
  
  f<-function(eps){
    p_val=bootstrap_test2(frequency,kmin,kmax,scale=1,nSimulation,nDirections,eps,tol)
    return(p_val-alpha)
  }
  
  tryCatch({
    res=uniroot(f,lower=minEps,upper=maxEps,tol = alpha/100)
    min_eps=res$root
  },
  error=function(cond){
    vec=c(NA,distance,beta,n)
    names(vec)=c("min_eps","distance","beta","sample_size")
  },
  warning=function(cond){}
  )
  
  vec=c(min_eps,distance,beta,n)
  names(vec)=c("min_eps","distance","beta","sample_size")
  return(vec)
}
