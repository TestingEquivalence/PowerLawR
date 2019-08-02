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

distanceCDF<-function(F1,F2){
  v=F1-F2
  v=v*v
  return(sqrt(sum(v)))
}

nearestPowerLaw<-function(CDF, kmin,kmax,betaLower, betaUpper){
  f<-function(beta){
    F1=powerLawCDF(beta,kmin,kmax)
    return(distanceCDF(F1,CDF))
  }
  
  res=optimize(f,c(betaLower,betaUpper),lower=betaLower,upper = betaUpper)
  return(res)
}