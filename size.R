fullToss<-function(i,parameter){
  set.seed(i*1000000)
  counting=rmultinom(n=1,size=parameter$n,prob=parameter$p)
  
  if (parameter$test=="asymptotic"){
    res=asymptotic_test(alpha = parameter$alpha,frequency = counting,
                        kmin = parameter$kmin,tol = parameter$tol)
    return(res)
  }
  
  if (parameter$test=="bootstrap"){
    res=bootstrap_test(alpha = parameter$alpha, frequency = counting,
                       kmin = parameter$kmin,
                       nSimulation = parameter$nSimulation, tol=parameter$tol)
    return(res)
  }
  
  if (parameter$test=="MLE"){
    res=powerLawMLE(counting,kmin=parameter$kmin,kmax=parameter$kmax,1,3)
    return(res)
  }

  return(NA)
}

sizeAtPowerLaw<-function(parameter){
  #calculate density of discrete power law
  parameter$kmin=parameter$kmin/parameter$scale
  parameter$kmax=parameter$kmax/parameter$scale
  parameter$scale=1

  parameter$p=powerLawDensity(beta = parameter$beta,
                              kmin = parameter$kmin,
                              kmax = parameter$kmax)
  i=c(1:parameter$nSamples)
  
  # simulate tests
  # v=sapply(i, fullToss,parameter)
  
  cl=getCluster()
  v=parSapply(cl,i, fullToss,parameter)
  stopCluster(cl)
  
  return(v)
}

