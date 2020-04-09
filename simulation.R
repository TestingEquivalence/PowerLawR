library(parallel)

# Calculate the number of cores
getCluster<-function(){
  no_cores <- detectCores() - 1
  
  # Initiate cluster
  cl <- makeCluster(no_cores,'SOCK')
  clusterExport(cl,c("powerLawDensity","powerLawCDF","l2","nearestPowerLaw","derivative",
                     "asympt_stdev","asymptotic_test","bootstrap_test",
                     "list2freq","fullToss","toss","powerLawLikelihood",
                     "bootstrap_stdev1","bootstrap_test1","test","powerLawMLE"))
  
  return(cl)
}



toss<-function(i,p, n, kmin, kmax,scale, eps,alpha, bootstrap, nSimulation,tol){
  counting=rmultinom(n=1,size=n,prob=p)
  if (bootstrap){
    res=bootstrap_test(alpha,counting,kmin,kmax,scale,nSimulation,tol)
    return(res[1]<eps)
  }
  else {
    res=asymptotic_test(alpha,counting,kmin,kmax,scale,tol)
    return(res[1]<eps)
  }
}

powerAtPoint<-function(p, n,  nSamples,  kmin, kmax,scale, eps,alpha,
                       bootstrap, nSimulation, tol){
  set.seed(01082019)
  i=c(1:nSamples)
  v=sapply(i, toss,p,n,kmin,kmax,scale,eps,alpha, bootstrap, nSimulation, tol)
  res=sum(v==TRUE)/nSamples
  return(res)
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