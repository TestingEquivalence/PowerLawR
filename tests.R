
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

bootstrap_stdev<-function(p,n,nSimulation,kmin,kmax, tol){
  
  i=c(1:nSimulation)
  f<-function(k){
    v=rmultinom(n=1,size=n,prob=p)
    v=v/sum(v)
    cdf=cumsum(v)
    res = nearestPowerLaw(cdf,kmin,kmax,1,3,tol=tol)
    distance=res$objective
    return(distance*distance)
  }
  
  sample=sapply(i,f)
  return(sqrt(var(sample)))
}

asymptotic_test<-function(alpha, frequency, kmin, kmax, scale)
{
  #calcualte cdf
  n=sum(frequency)
  frequency=frequency/n
  cdf=cumsum(frequency)
  kmin=kmin/scale
  kmax=kmax/scale
  
  res = nearestPowerLaw(cdf,kmin,kmax,1,3)
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

bootstrap_test2<-function(alpha, frequency, kmin, kmax,
                         scale,nSimulation, tol=NA)
{
  #calcualte cdf
  n=sum(frequency)
  p=frequency/n
  cdf=cumsum(p)
  kmin=kmin/scale
  kmax=kmax/scale
  
  res = nearestPowerLaw(cdf,kmin,kmax,1,3)
  beta=res$minimum
  distance=res$objective
  
  vol=bootstrap_stdev(p,n,nSimulation,kmin,kmax,tol)
  qt=qnorm(1-alpha,0,1)
  min_eps = distance*distance + qt*vol
  min_eps=sqrt(min_eps)
  
  vec=c(min_eps,distance,beta,n)
  names(vec)=c("min_eps","distance","beta","sample_size")
  return(vec)
}

bootstrap_test<-function(alpha, frequency, kmin, kmax,
                         scale,nSimulation, tol=NA)
{
  #calcualte cdf
  n=sum(frequency)
  p=frequency/n
  cdf=cumsum(p)
  kmin=kmin/scale
  kmax=kmax/scale
  
  #calculate  distance
  res = nearestPowerLaw(cdf,kmin,kmax,1,3)
  beta=res$minimum
  distance=res$objective
  
  #compute bootstrap distribution
  i=c(1:nSimulation)
  f<-function(k){
    v=rmultinom(n=1,size=n,prob=p)
    v=v/sum(v)
    cdf=cumsum(v)
    res = nearestPowerLaw(cdf,kmin,kmax,1,3,tol=tol)
    return(res$objective)
  }
  
  sample=sapply(i,f)
  qt=quantile(sample, probs = c(alpha))
  mu=mean(sample)
  min_eps=distance+distance-qt
  
  vec=c(min_eps,distance,beta,n)
  names(vec)=c("min_eps","distance","beta","sample_size")
  return(vec)
}

fmultiple<-function(row, kmins, kmaxs, alpha, scale, 
                    counting,bootstrap, nSimulation){
  kmin=kmins[row[1]]
  kmax=kmaxs[row[2]]
  frequency=list2freq(counting,kmin,kmax,scale)
  
  if (bootstrap){
    set.seed(30062020)
    res=bootstrap_test(alpha,frequency,kmin,kmax,scale,nSimulation)
  }
  else {
    res=asymptotic_test(alpha,frequency,kmin,kmax,scale)
  }
  
  print(paste("done: ","kmin=",kmin," kmax=", kmax))
  return(c(row[1],row[2],res))
}

multiple_test <- function(alpha, counting, kmins, kmaxs,
                          scale,bootstrap=FALSE, nSimulation=0) {
  nrow=length(kmins)
  ncol = length(kmaxs)
  min_eps=matrix(data=NA,nrow,ncol)
  beta=matrix(data=NA,nrow,ncol)
  distance=matrix(data=NA,nrow,ncol)
  sample_size=matrix(data=NA,nrow,ncol)
  
  
  rownames(min_eps)=kmins
  rownames(beta)=kmins
  rownames(distance)=kmins
  rownames(sample_size)=kmins
  
  
  colnames(min_eps)=kmaxs
  colnames(beta)=kmaxs
  colnames(distance)=kmaxs
  colnames(sample_size)=kmaxs
  
  
  i=c(1:nrow)
  j=c(1:ncol)
  grd=expand.grid(i,j)
  colnames(grd)=c("i","j")
  
  # cl=getCluster()
  # clusterExport(cl,c("fmultiple"))
  # ls=parApply(cl,grd, 1, fmultiple, kmins,kmaxs,alpha,scale, 
  #             counting, bootstrap,nSimulation)
  # stopCluster(cl)
  
  ls=apply(grd, 1, fmultiple, kmins,kmaxs,alpha,scale, 
              counting, bootstrap,nSimulation)
  
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
