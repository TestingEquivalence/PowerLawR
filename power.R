
uniformRandomStress<-function(kmin,kmax,beta,eps){
  p=powerLawDensity(beta,kmin,kmax)
  repeat{
    v=runif(kmax-kmin+1,0.9,1.1)
    v=v*p
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


boundaryPower<-function(n,eps,kmin,kmax,scale,beta,alpha, boundaryPointType,
                        bootstrap, nSimulation, tol, adjEps=1){
  
  
  kmin=kmin/scale
  kmax=kmax/scale
  
  i=c(1:100)
  
  if (boundaryPointType==1){
    f<-function(i){
      powerLawStress(n,eps,kmin,kmax,beta)
    }  
  }
  
  if (boundaryPointType==2){
    f<-function(i){
      uniformRandomStress(kmin,kmax,beta,eps)
    }
  }
  
  set.seed(01012020)
  rndBndPoints=lapply(i, f)
  
  res=powerAtPoints(rndBndPoints,n,nSamples = 1000,
                    kmin,kmax,1,eps*adjEps,alpha, bootstrap, nSimulation,tol)
  return(res)
}

