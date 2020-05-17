
derivative<-function(F1,F2)
{
  Z=2*(F1-F2)
  Z=rev(Z)
  Z=cumsum(Z)
  Z=rev(Z)
  return(Z)
}

asympt_stdev<-function(p,derivative){
  vec = derivative
  vnsq_1  = sum(p*vec*vec)
  
  k=length(p)
  vnsq_2=0
  for (j1 in 1:k)
    for (j2 in 1:k)
      vnsq_2 = vnsq_2 + vec[j1] * vec[j2] * p[j1] * p[j2]
  
  
  vnsq  = vnsq_1 - vnsq_2
  return (sqrt(vnsq))
}

asymptotic_test<-function(alpha, frequency, kmin, kmax, scale, tol=NA)
{
  #calcualte cdf
  n=sum(frequency)
  frequency=frequency/n
  cdf=cumsum(frequency)
  kmin=kmin/scale
  kmax=kmax/scale
  
  res = nearestPowerLaw(cdf,kmin,kmax,1,3, tol)
  beta=res$minimum
  distance=res$objective
  pLawCDF=powerLawCDF(beta,kmin,kmax)
  
  drv=derivative(cdf,pLawCDF)
  vol=asympt_stdev(frequency,drv)
  vol=vol/ sqrt(n)
  
  qt=qnorm(1-alpha,0,1)
  min_eps = distance*distance + qt*vol
  min_eps=sqrt(min_eps)
  
  vec=c(min_eps,distance,beta,n)
  names(vec)=c("min_eps","distance","beta","sample_size")
  return(vec)
}

bootstrap_test<-function(alpha, frequency, kmin, kmax,
                          scale,nSimulation, tol=NA){
  bootstrap_test1(alpha, frequency, kmin, kmax, scale, nSimulation, tol)  
}
  

fmultiple<-function(row,parameter){
  kmin=parameter$kmins[row[1]]
  kmax=parameter$kmaxs[row[2]]
  frequency=list2freq(parameter$counting,kmin,kmax,parameter$scale)
  
  if (parameter$test=="asymptotic"){
    res=asymptotic_test(alpha = parameter$alpha,frequency,
                        kmin,kmax,
                        scale = parameter$scale,tol = parameter$tol)
  }
  
  if (parameter$test=="bootstrap1"){
    set.seed(30062020)
    res=bootstrap_test1(alpha = parameter$alpha, frequency,
                        kmin,kmax,
                        scale = parameter$scale, nSimulation = parameter$nSimulation,
                        tol=parameter$tol)
  }
  
  if (parameter$test=="bootstrap2"){
    set.seed(30062020)
    res=bootstrap_test2_1(alpha=parameter$alpha,frequency=frequency, 
                          kmin=kmin, kmax=kmax,
                          scale=parameter$scale,
                          nSimulation=parameter$tol, tol=parameter$tol, 
                          nDirections=parameter$nDirections)
  }
  
  return(c(row[1],row[2],res))
}

multiple_test <- function(parameter) {
  nrow=length(parameter$kmins)
  ncol = length(parameter$kmaxs)
  min_eps=matrix(data=NA,nrow,ncol)
  beta=matrix(data=NA,nrow,ncol)
  distance=matrix(data=NA,nrow,ncol)
  sample_size=matrix(data=NA,nrow,ncol)
  
  rownames(min_eps)=parameter$kmins
  rownames(beta)=parameter$kmins
  rownames(distance)=parameter$kmins
  rownames(sample_size)=parameter$kmins
  
  colnames(min_eps)=parameter$kmaxs
  colnames(beta)=parameter$kmaxs
  colnames(distance)=parameter$kmaxs
  colnames(sample_size)=parameter$kmaxs
  
  i=c(1:nrow)
  j=c(1:ncol)
  grd=expand.grid(i,j)
  colnames(grd)=c("i","j")
  
  # cl=getCluster()
  # clusterExport(cl,c("fmultiple"))
  # ls=parApply(cl,grd, 1, fmultiple,parameter)
  # stopCluster(cl)
  
  ls=apply(grd, 1, fmultiple, parameter)
 
  for (rn in c(1:ncol(ls))){
    r=ls[,rn]
    i=r[1]
    j=r[2]
    min_eps[i,j]=r[3]
    distance[i,j]=r[4]
    beta[i,j]=r[5]
    sample_size[i,j]=r[6]
  }

  ls=list(min_eps=min_eps,distance=distance,beta=beta, sample_size=sample_size)
  return(ls)
}
multiple_MLE <- function(parameter) {
  nrow=length(parameter$kmins)
  ncol = length(parameter$kmaxs)
  beta=matrix(data=NA,nrow,ncol)
  sample_size=matrix(data=NA,nrow,ncol)
  rownames(beta)=parameter$kmins
  rownames(sample_size)=parameter$kmins
  colnames(beta)=parameter$kmaxs
  colnames(sample_size)=parameter$kmaxs
  
  
  i=c(1:nrow)
  j=c(1:ncol)
  grd=expand.grid(i,j)
  colnames(grd)=c("i","j")
  
  for (i in c(1:nrow)){
    for (j in c(1:ncol)){
      kmin=parameter$kmins[i]
      kmax=parameter$kmaxs[j]
      frq=list2freq(parameter$counting,kmin,kmax,parameter$scale)
      sample_size[i,j]=sum(frq)
      res=powerLawMLE(frq,kmin/parameter$scale,kmax/parameter$scale,1,3)
      beta[i,j]=res$minimum
    }
  }
  
  ls=list(beta=beta, sample_size=sample_size)
  return(ls)
}
