library(parallel)

# Calculate the number of cores
getCluster<-function(){
  no_cores <- detectCores() - 1
  
  # Initiate cluster
  cl <- makeCluster(no_cores,'SOCK')
  clusterExport(cl,c("powerLawDensity","powerLawCDF","l2","nearestPowerLaw","derivative",
                     "asympt_stdev","asymptotic_test","bootstrap_test",
                     "list2freq","fullToss","powerLawLikelihood",
                     "bootstrap_stdev1","bootstrap_test1","powerLawMLE",
                     "bootstrap_test2","bootstrap_test3","test",
                     "closeRandomPoint","linComb","linearBoundaryPoint"))
  
  return(cl)
}

test<-function(counting,parameter){
  if (parameter$test=="asymptotic"){
    res=asymptotic_test(alpha = parameter$alpha,frequency = counting,
                        kmin = parameter$kmin,kmax =  parameter$kmax,
                        scale = parameter$scale,tol = parameter$tol)
    return(res[1]<parameter$eps)
  }
  
  if (parameter$test=="bootstrap1"){
    res=bootstrap_test1(alpha = parameter$alpha, frequency = counting,
                        kmin = parameter$kmin,kmax = parameter$kmax,
                        scale = parameter$scale, nSimulation = parameter$nSimulation,
                        tol=parameter$tol)
    return(res[1]<parameter$eps)
  }
  
  if (parameter$test=="bootstrap2"){
    pval=bootstrap_test2(frequency = counting, kmin=parameter$kmin, kmax=parameter$kmax,
                        scale = parameter$scale,nSimulation = parameter$nSimulation,
                        nDirections = parameter$nDirections,eps=parameter$eps,
                        tol=parameter$tol)
    return(pval<parameter$alpha)
  }
  
  return(NA)
}

powerAtPoint<-function(p,parameter){
  set.seed(01082019)
  points=list()
  for (i in c(1:parameter$nSamples)){
    points[[i]]=rmultinom(n=1,size=parameter$n,prob=p)  
  }
  
  v=sapply(points, test,parameter)
  res=sum(v==TRUE)/parameter$nSamples
  return(res)
}

powerAtPoints<-function(points, parameter){
  cl=getCluster()
  v=parSapply(cl,points,powerAtPoint,parameter)
  stopCluster(cl)
  
  #v=sapply(points,powerAtPoint,parameter)
  return(v)
}