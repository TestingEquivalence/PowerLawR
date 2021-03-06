asymptotic="asymptotic"
bootstrap="bootstrap"

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
  
  f<-function(j){
    v=vec[j]*vec
    v=v*p[j]
    v=v*p
    return(sum(v))
  }
  
  vv=sapply(c(1:k),f)
  vnsq_2=sum(vv)

  vnsq  = vnsq_1 - vnsq_2
  return (sqrt(vnsq))
}

#' The asymptotic test is based on the asymptotic distribution of the test statistic. 
#' The test statistic is the minimum Euclidean distance between the empirical CDF 
#' and the family of power laws with given lower and upper bounds.
#' The test should be used carefully because it is approximate 
#' and may be anti-conservative at some points. 
#' In order to obtain a conservative test reducing of alpha  (usually halving) or
#' slight shrinkage of the tolerance parameter may be appropriate.
#' \code{asymptotic_test} asymptotic equivalence test for power law
#' @param alpha significance level
#' @param frequency vector of observed counting frequencies
#' @param kmin lower bound of the power law
#' @param tol optional tolerance parameter for the numerical optimization 
#' to find the closest power law 
#' @return test returns:
#' min_eps - the minimim tolerance parameter fo which H0 can be rejected
#' distance -Euclidean distance between empirical CDF and CDF of the closest power law
#' beta - minimum distance estimate of power law exponent
#' sample_size - sample size
asymptotic_test<-function(alpha, frequency, kmin, tol=NA)
{
  #calcualte cdf
  n=sum(frequency)
  frequency=frequency/n
  cdf=cumsum(frequency)
 
  kmax=length(frequency)+kmin-1
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

fmultiple<-function(row,parameter){
  kmin=parameter$kmins[row[1]]
  kmax=parameter$kmaxs[row[2]]
  frequency=list2freq(parameter$counting,kmin,kmax,parameter$scale)
  
  if (parameter$test=="asymptotic"){
    res=asymptotic_test(alpha = parameter$alpha,frequency,
                        kmin=kmin/parameter$scale)
  }
  
  if (parameter$test=="bootstrap"){
    set.seed(30062020)
    res= bootstrap_test(alpha = parameter$alpha, frequency,
                        kmin=kmin/parameter$scale,
                        nSimulation = parameter$nSimulation)
  }
  
  return(c(row[1],row[2],res))
}

#' 
#' \code{multiple_test} Convenient function to perform multiple equivalence tests on the same data.
#' It also transforms the usual counting data to frequency data. 
#' Usually we observe the counting data only and also do not know 
#' the upper and lower bound of the power low. 
#' In this case we need to transform the counting date to frequencies . 
#' We also may perform multiple equivalence tests for different values of upper and lower bounds. 
#' The convenient function "multiple_test" performs all these tasks efficiently 
#' using multiple cores for computing.
#' @param parameter The parameter should be a list (s3 object) 
#' containing following fields (see example.R):
#' scale - scaling, which maybe necessary to make computations feasible
#' alpha - significance level 
#' nSimulation - number of bootstrap replications  
#' test - string, should be asymptotic or bootstrap
#' kmins - vector of possible lower bounds of power law
#' kmaxs - vector of possible upper bounds of power law
#' counting - counting data  where the upper and lower bounds are unknown 
#' @return test returns list of four tables:
#' beta - table of the estimated beta's (minimum distance estimate of power law exponent)
#' distance - table of Euclidean distances between empirical CDF and CDF of the closest power law
#' sample_size - table of sample sizes
#' min_eps - table of the minimum tolerance parameters for which H0 can be rejected

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
  
  cl=getCluster()
  clusterExport(cl,c("fmultiple"))
  ls=parApply(cl,grd, 1, fmultiple,parameter)
  stopCluster(cl)

  # ls=apply(grd, 1, fmultiple, parameter)
 
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
