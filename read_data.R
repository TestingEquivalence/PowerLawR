readVector<-function(pathOfFile){
  v=read.csv2(pathOfFile,header = FALSE)
  return(v$V1)
}

list2freq<-function(v, kmin, kmax, scale){
  #cut the list and rescale
  v=v[v>=kmin]
  v=v[v<=kmax]
  v=v/scale
  
  kmin=as.integer(kmin/scale)
  v=v-kmin
  
  kmax=as.integer(kmax/scale)
  
  nbins=kmax-kmin+1
  
  res=tabulate(v, nbins)
  return(res)
}