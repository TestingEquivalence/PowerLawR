
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

boundaryPower<-function(parameter){
  
  
  f<-function(i){
    uniformRandomStress(kmin=parameter$kmin/parameter$scale,
                        kmax=parameter$kmax/parameter$scale,
                        parameter$beta,
                        parameter$eps)
  }
  
  # initial seed
  set.seed(01012020)
 
  res=c()
  
  i=c(1:100)
  points=lapply(i, f)
  parameterAdj=parameter
  parameterAdj$eps=parameter$eps*parameter$epsAdj
  res=powerAtPoints(points, parameter=parameterAdj)
  
  return(res)
}

