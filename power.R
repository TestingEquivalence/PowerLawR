closeRandomPoint<-function(p,n, eps){
  repeat{
    rtab=as.vector(tab)
    v=rmultinom(n=1,size=n,prob=p)
    v=v/n
    v=cumsum(v)
    res=nearestPowerLaw(v,kmin,kmax,1,3)
    if (res$objective>eps)
      if (res$minimum>1)
        if (res$minimum<3)
          return(res)
  }
}

boundaryPoint<-function(i,l,eps,kmin,kmax){
  p=closeRandomPoint(tab,eps,distance)
  n=sum(tab)
  tab=tab/n
  
  if (identical(distance,cond_l2)){
    q=p2triangle(startValue(tab))
  }
  
  if (identical(distance,min_l2)){
    q=min_l22(tab)$par
    q=p2triangle(q)
  }
  
  res=linearBoundaryPoint(p,q,eps,distance)
  return(res)
}
