library(parallel)

# Calculate the number of cores
getCluster<-function(){
  no_cores <- detectCores() - 1
  
  # Initiate cluster
  cl <- makeCluster(no_cores,'SOCK')
  clusterExport(cl,c("powerLawDensity","powerLawCDF","l2","nearestPowerLaw","derivative",
                     "asympt_stdev","asymptotic_test","bootstrap_test",
                     "list2freq","fullToss","powerLawLikelihood",
                     "bootstrap_stdev1","bootstrap_test1","test","powerLawMLE"))
  
  return(cl)
}

powerAtPoint<-function(p, n,  nSamples,  kmin, kmax,scale,test,eps){
  set.seed(01082019)
  points=list()
  for (i in c(1:nSamples)){
    points[[i]]=rmultinom(n=1,size=n,prob=p)  
  }
  
  v=sapply(points, test,kmin,kmax,scale,eps)
  res=sum(v==TRUE)/nSamples
  return(res)
}

powerAtPoints<-function(points, n,  nSamples,  kmin, kmax,scale,test,eps){
  cl=getCluster()
  v=parSapply(cl,points,powerAtPoint,n,nSamples,kmin,kmax,scale,test,eps)
  stopCluster(cl)
  
  #v=sapply(points,powerAtPoint,n,nSamples,kmin,kmax,scale,test,eps)
  return(v)
}