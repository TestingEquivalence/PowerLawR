powerLawDensity<-function(beta, kmin, kmax){
  v=c(kmin:kmax)
  v=v ^ (-beta)
  v=v /sum(v)
  return (v)
}

powerLawCDF<-function(beta, kmin, kmax){
  v=powerLawDensity(beta,kmin,kmax)
  v=cumsum(v)
  return(v)
}

l2<-function(F1,F2){
  v=F1-F2
  v=v*v
  return(sqrt(sum(v)))
}

nearestPowerLaw<-function(CDF, kmin,kmax,betaLower, betaUpper, tol=NA){
  f<-function(beta){
    F1=powerLawCDF(beta,kmin,kmax)
    return(l2(F1,CDF))
  }
  
  if (is.na(tol)){
    res=optimize(f,c(betaLower,betaUpper),lower=betaLower,upper = betaUpper)
  }
  else{
    res=optimize(f,c(betaLower,betaUpper),lower=betaLower,upper = betaUpper, tol=tol)
  }
    
  return(res)
}

powerLawLikelihood<-function(beta, kmin, kmax,frequency){
  dens=powerLawDensity(beta,kmin,kmax)
  dens=log(dens)
  return(sum(dens*frequency))
}

powerLawMLE<-function(frequency, kmin,kmax,betaLower, betaUpper, tol=NA){
  f<-function(beta){
    -powerLawLikelihood(beta, kmin, kmax,frequency)
  }
  
  if (is.na(tol)){
    res=optimize(f,c(betaLower,betaUpper),lower=betaLower,upper = betaUpper)
  }
  else{
    res=optimize(f,c(betaLower,betaUpper),lower=betaLower,upper = betaUpper, tol=tol)
  }
  return(res)
}