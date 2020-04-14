
uniformRandomStress<-function(kmin,kmax,beta,eps){
  p=powerLawDensity(beta,kmin,kmax)
  repeat{
    v=runif(length(p))
    v=v*p
    v=v/sum(v)
    cv=cumsum(v) 
    res=nearestPowerLaw(cv,kmin,kmax,1,3)
    if (res$objective>eps){
      lp=linearBoundaryPoint(p,v,eps,kmin,kmax)
      res=nearestPowerLaw(cumsum(lp),kmin,kmax,1,3)
      if (res$minimum>beta*0.8)
        if (res$minimum<beta*1.2)
          return(lp)
    }
      
  }
}


powerLawStress<-function(n,eps,kmin,kmax,beta){
  p=powerLawDensity(beta,kmin,kmax)
  q=closeRandomPoint(p,n, eps,beta,kmin,kmax)
  lp=linearBoundaryPoint(p,q,eps,kmin,kmax)
  res=nearestPowerLaw(cumsum(lp),kmin,kmax,1,3)
  return(lp)
}


boundaryPower<-function(n,eps,kmin,kmax,scale,beta,nSamples,
                        boundaryPointType,epsAdj){
  
  i=c(1:100)
  
  if (boundaryPointType==1){
    f<-function(i){
      powerLawStress(n,eps,kmin=kmin/scale,kmax=kmax/scale,beta)
    }  
  }
  
  if (boundaryPointType==2){
    f<-function(i){
      uniformRandomStress(kmin=kmin/scale,kmax=kmax/scale,beta,eps)
    }
  }
  
  set.seed(01012020)
  points=lapply(i, f)
  
  eps=eps*epsAdj
  res=powerAtPoints(points, n,  nSamples,  kmin, kmax,scale,test,eps)
  return(res)
}

