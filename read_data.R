readVector<-function(pathOfFile){
  v=read.csv2(pathOfFile,header = FALSE)
  return(v$V1)
}

# res=readVector("C:\\Users\\Ostrovski\\Google Drive\\Writing\\PowerLaw\\CitySize\\list_ge.csv")
