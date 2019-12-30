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

uniformRandomStress<-function(kmin,kmax,beta,eps){
  repeat{
    v=runif(kmax-kmin+1)
    v=v/sum(v)
    p=powerLawDensity(beta,kmin,kmax)
    v=linComb(p,v,0.2)
    
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

randomBoundaryPoint1<-function(n,eps,kmin,kmax,beta){
  p=powerLawDensity(beta,kmin,kmax)
  q=closeRandomPoint(p,n, eps,beta,kmin,kmax)
  lp=linearBoundaryPoint(p,q,eps,kmin,kmax)
  res=nearestPowerLaw(cumsum(lp),kmin,kmax,1,3)
  return(lp)
}

randomBoundaryPoint2<-function(n,eps,kmin,kmax,beta,p){
  q=closeRandomPoint(p,n,eps,beta)
  lp=linearBoundaryPoint(p,q,eps,kmin,kmax)
  return(lp)
}

randomBoundaryPoint3 <-function(n,eps,kmin,kmax,beta){
  p=powerLawDensity(beta,kmin,kmax)
  q=uniformRandomStress(kmin,kmax,beta,eps)
  lp=linearBoundaryPoint(p,q,eps,kmin,kmax)
  return(lp)
} 

boundaryPower<-function(n,eps,kmin,kmax,scale,beta,alpha, boundaryPointType,p=NULL){
  kmin=kmin/scale
  kmax=kmax/scale
  
  i=c(1:100)
  
  if (boundaryPointType==1){
    f<-function(i){
      randomBoundaryPoint1(n,eps,kmin,kmax,beta)
    }  
  }
  
  if (boundaryPointType==2){
    f<-function(i){
      randomBoundaryPoint2(n,eps,kmin,kmax,beta,p)
    }
  }
  
  if (boundaryPointType==3){
    f<-function(i){
      randomBoundaryPoint3(n,eps,kmin,kmax,beta)
    }
  }
  
  set.seed(01012020)
  rndBndPoints=lapply(i, f)
  
  res=powerAtPoints(rndBndPoints,n,nSamples = 1000,kmin,kmax,1,eps,alpha)
  return(res)
}
