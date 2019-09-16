library(parallel)
source("PowerLaw.R")
source("asymptotic_test.R")

# Calculate the number of cores
getCluster<-function(){
  no_cores <- detectCores() - 1
  
  # Initiate cluster
  cl <- makeCluster(no_cores,'SOCK')
  clusterExport(cl,c("powerLawDensity","powerLawCDF","l2","nearestPowerLaw","derivative","asympt_stdev","asymptotic_test",
                     "list2freq","fullToss","toss"))
  
  return(cl)
}

fullToss<-function(i,p, n, kmin, kmax,scale, alpha){
  set.seed(i)
  counting=rmultinom(n=1,size=n,prob=p)
  res=asymptotic_test(alpha,counting,kmin,kmax,scale)
  return(res)
}

toss<-function(i,p, n, kmin, kmax,scale, eps,alpha){
  counting=rmultinom(n=1,size=n,prob=p)
  res=asymptotic_test(alpha,counting,kmin,kmax,scale)
  return(res[1]<=eps)
}

powerAtPoint<-function(p, n,  nSamples,  kmin, kmax,scale, eps,alpha){
  set.seed(01082019)
  i=c(1:nSamples)
  v=sapply(i, toss,p,n,kmin,kmax,scale,eps,alpha)
  return(sum(v==TRUE)/nSamples)
}

powerAtPoints<-function(points, n,  nSamples,  kmin, kmax,scale, eps,alpha){
  cl=getCluster()
  v=parSapply(cl,points,powerAtPoint,n,nSamples,kmin,kmax,scale,eps,alpha)
  stopCluster(cl)
  return(v)
}
 