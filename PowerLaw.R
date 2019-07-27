powerLawDensity<-function(beta, kmin, kmax){
  v=c(kmin:kmax)
  Math.Pow(k, -b)
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

distanceDensity<-function(f1,f2){
  F1=cumsum(f1)
  F2=cumsum(f2)
  return(distanceCDF(F1,F2))
}