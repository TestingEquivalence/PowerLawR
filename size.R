fullToss<-function(i,parameter){
  set.seed(i*1000000)
  counting=rmultinom(n=1,size=parameter$n,prob=parameter$p)
  
  if (parameter$test=="asymptotic"){
    res=asymptotic_test(alpha = parameter$alpha,frequency = counting,
                        kmin = parameter$kmin,kmax =  parameter$kmax,
                        scale = parameter$scale,tol = parameter$tol)
    return(res)
  }
  
  if (parameter$test=="bootstrap1"){
    res=bootstrap_test1(alpha = parameter$alpha, frequency = counting,
                        kmin = parameter$kmin,kmax = parameter$kmax,
                        scale = parameter$scale, nSimulation = parameter$nSimulation,
                        tol=parameter$tol)
    return(res)
  }
  
  if (parameter$test=="bootstrap2"){
    res=bootstrap_test2(frequency = counting, kmin=parameter$kmin, kmax=parameter$kmax,
                        scale = parameter$scale,nSimulation = parameter$nSimulation,
                        nDirections = parameter$nDirections,eps=parameter$eps,
                        tol=parameter$tol)
    return(res)
  }
  
  return(NA)
}

sizeAtPowerLaw<-function(parameter){
  #calculate density of discrete power law
  kmin=parameter$kmin/parameter$scale
  kmax=parameter$kmax/parameter$scale
  beta=parameter$beta
  parameter$p=powerLawDensity(beta,kmin,kmax)
  
  #simulate tests
  #v=sapply(i, fullToss,p,n,kmin,kmax,scale)
  i=c(1:nSamples)
  cl=getCluster()
  v=parSapply(cl,i, fullToss,parameter)
  stopCluster(cl)
  
  return(v)
}

