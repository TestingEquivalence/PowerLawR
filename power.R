
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


boundaryPower<-function(parameter){
  
  i=c(1:100)
  
  if (parameter$boundaryPointType==1){
    f<-function(i){
      powerLawStress(n=parameter$n,eps=parameter$eps,
                     kmin=parameter$kmin/parameter$scale,
                     kmax=parameter$kmax/parameter$scale,parameter$beta)
    }  
  }
  
  if (parameter$boundaryPointType==2){
    f<-function(i){
      uniformRandomStress(kmin=parameter$kmin/parameter$scale,
                          kmax=parameter$kmax/parameter$scale,
                          parameter$beta,
                          parameter$eps)
    }
  }
  
  set.seed(01012020)
  points=lapply(i, f)
  
  parameter$eps=parameter$eps*parameter$epsAdj
  res=powerAtPoints(points, parameter)
  return(res)
}

