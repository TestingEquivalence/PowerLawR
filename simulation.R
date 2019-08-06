library(parallel)
source("PowerLaw.R")
source("asymptotic_test.R")

# Calculate the number of cores
getCluster<-function(){
  no_cores <- detectCores() - 1
  
  # Initiate cluster
  cl <- makeCluster(no_cores,'SOCK')
  clusterExport(cl,c("powerLawDensity","powerLawCDF","distanceCDF","nearestPowerLaw","derivative","asympt_stdev","asymptotic_test",
                     "list2freq"))
  
  return(cl)
}

 