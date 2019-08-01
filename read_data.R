readVector<-function(pathOfFile){
  v=read.csv2(pathOfFile,header = FALSE)
  return(v$V1)
}

list2freq<function(v, kmin, kmax, scale){
  #cut the list and rescale
  v=v[v>=kmin]
  v=v[v<=kmax]
  v=v/scale
  v=round(v,0)
  v=as.integer(v)
  
  
  
  return(v)
}