library(parallel)

# Calculate the number of cores
getCluster<-function(){
  no_cores <- detectCores() - 1
  
  # Initiate cluster
  cl <- makeCluster(no_cores,'SOCK')
  clusterExport(cl,c("powerLawDensity","powerLawCDF","l2","nearestPowerLaw","derivative",
                     "asympt_stdev","asymptotic_test","bootstrap_stdev","bootstrap_test",
                     "list2freq","fullToss","toss"))
  
  return(cl)
}

fullToss<-function(i,p, n, kmin, kmax,scale, alpha,bootstrap, nSimulation){
  set.seed(i)
  counting=rmultinom(n=1,size=n,prob=p)
  if (bootstrap){
    res=bootstrap_test(alpha,counting,kmin,kmax,scale,nSimulation)
  }
  else {
    res=asymptotic_test(alpha,counting,kmin,kmax,scale)
  }
  return(res)
}

toss<-function(i,p, n, kmin, kmax,scale, eps,alpha, bootstrap, nSimulation){
  counting=rmultinom(n=1,size=n,prob=p)
  if (bootstrap){
    res=bootstrap_test(alpha,counting,kmin,kmax,scale,nSimulation)
  }
  else {
    res=asymptotic_test(alpha,counting,kmin,kmax,scale)
  }
  return(res[1]<=eps)
}

powerAtPoint<-function(p, n,  nSamples,  kmin, kmax,scale, eps,alpha,
                       bootstrap, nSimulation, bType){
  set.seed(01082019)
  i=c(1:nSamples)
  v=sapply(i, toss,p,n,kmin,kmax,scale,eps,alpha, bootstrap, nSimulation)
  return(sum(v==TRUE)/nSamples)
}

powerAtPoints<-function(points, n,  nSamples,  kmin, kmax,scale, eps,alpha,
                        bootstrap, nSimulation, bType){
  cl=getCluster()
  v=parSapply(cl,points,powerAtPoint,n,nSamples,kmin,kmax,
              scale,eps,alpha, bootstrap, nSimulation)
  stopCluster(cl)
  return(v)
}
 