closeRandomPoint<-function(p,n, eps){
  repeat{
    v=rmultinom(n=1,size=n,prob=p)
    v=v/n
    v=cumsum(v)
    res=nearestPowerLaw(v,kmin,kmax,1,3)
    if (res$objective>eps)
      if (res$minimum>1)
        if (res$minimum<3)
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

randomBoundaryPoint<-function(n,eps,kmin,kmax,beta){
  p=powerLawDensity(beta,kmin,kmax)
  q=closeRandomPoint(p,n,eps)
  lp=linearBoundaryPoint(p,q,eps,kmin,kmax)
  return(lp)
}

boundaryPower<-function(n,eps,kmin,kmax,scale,beta,alpha){
  kmin=kmin/scale
  kmax=kmax/scale
  
  i=c(1:100)
  
  f<-function(i){
    randomBoundaryPoint(n,eps,kmin,kmax,beta)
  }
  
  set.seed(01102019)
  rndBndPoints=lapply(i, f)
  
  powerAtPoints(rndBndPoints,n,nSamples = 1000,kmin,kmax,1,eps,alpha)
}

