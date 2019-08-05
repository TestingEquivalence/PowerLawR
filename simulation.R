library(parallel)
source("PowerLaw.R")
source("asymptotic_test.R")

# Calculate the number of cores
getCluster<-function(){
  no_cores <- detectCores() - 1
  
  # Initiate cluster
  cl <- makeCluster(no_cores,'SOCK')
  clusterExport(cl,c("powerLawDensity"))
  
  return(cl)
}

 