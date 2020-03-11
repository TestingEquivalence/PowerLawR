library(parallel)

# Calculate the number of cores
getCluster<-function(){
  no_cores <- detectCores() - 1
  
  # Initiate cluster
  cl <- makeCluster(no_cores,'SOCK')
  clusterExport(cl,c("powerLawDensity","powerLawCDF","l2","nearestPowerLaw","derivative",
                     "asympt_stdev","asymptotic_test","bootstrap_stdev","bootstrap_test",
                     "list2freq","fullToss","toss","powerLawMLE","MLEToss","powerLawLikelihood"))
  
  return(cl)
}

fullToss<-function(i,p, n, kmin, kmax,scale, alpha,bootstrap, nSimulation,tol){
  set.seed(i)
  counting=rmultinom(n=1,size=n,prob=p)
  if (bootstrap){
    res=bootstrap_test(alpha,counting,kmin,kmax,scale,nSimulation,tol)
  }
  else {
    res=asymptotic_test(alpha,counting,kmin,kmax,scale,tol)
  }
  return(res)
}

toss<-function(i,p, n, kmin, kmax,scale, eps,alpha, bootstrap, nSimulation,tol){
  counting=rmultinom(n=1,size=n,prob=p)
  if (bootstrap){
    res=bootstrap_test(alpha,counting,kmin,kmax,scale,nSimulation,tol)
  }
  else {
    res=asymptotic_test(alpha,counting,kmin,kmax,scale,tol)
  }
  return(res[1]<=eps)
}

powerAtPoint<-function(p, n,  nSamples,  kmin, kmax,scale, eps,alpha,
                       bootstrap, nSimulation, tol){
  set.seed(01082019)
  i=c(1:nSamples)
  v=sapply(i, toss,p,n,kmin,kmax,scale,eps,alpha, bootstrap, nSimulation, tol)
  return(sum(v==TRUE)/nSamples)
}

powerAtPoints<-function(points, n,  nSamples,  kmin, kmax,scale, eps,alpha,
                        bootstrap, nSimulation,tol){
  cl=getCluster()
  v=parSapply(cl,points,powerAtPoint,n,nSamples,kmin,kmax,
              scale,eps,alpha, bootstrap, nSimulation,tol)
  # v=sapply(points,powerAtPoint,n,nSamples,kmin,kmax,
  #             scale,eps,alpha, bootstrap, nSimulation,tol)
  stopCluster(cl)
  return(v)
}
MLEToss<-function(i,p, n, kmin, kmax){
  set.seed(i)
  frq=rmultinom(n=1,size=n,prob=p)
  res=powerLawMLE(frq,kmin,kmax,1,3)
  return(res$minimum)
} 