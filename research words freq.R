source("PowerLaw.R")
source("tests.R")
source("power.R")
source("read_data.R")
source("simulation.R")
source("size.R")
source("bootstrap_tests.R")

#vector of city sizes in Germany
words=readVector("C:\\Users\\Ostrovski\\Google Drive\\Writing\\PowerLaw\\WordCounts\\words.csv")

#different k_min values
#kmins=c(500,600,700,800,1000)

#different k_max values, first value coincides with most frequent word
#kmaxs=c(max(words),2e6,5e6,10e6)

parameter=list(scale=1000,alpha=0.05, tol=0.001,nSimulation=1000, nDirections = 100, 
               test="asymptotic", 
               kmins=c(1000,2000,3000, 5000, 10000), 
               kmaxs=c(max(words),2e6,5e6,10e6), 
               counting=words)

parameter=list(scale=100,alpha=0.05, tol=0.001,nSimulation=1000, nDirections = 100, 
               test="asymptotic", 
               kmins=c(100,200,300, 500, 1000), 
               kmaxs=c(max(words),2e6,5e6,10e6), 
               counting=words)


#carry out multiple tests for the power law
#given test is computed for each combination of k_min and k_max 
#use all but one CPU cores

result=multiple_test(parameter)

#write test results
write.table(result$beta,paste("beta.csv"))
write.table(result$distance,paste("distance.csv"))
write.table(result$min_eps,paste("min_eps.csv"))
write.table(result$sample_size,paste("sample_size.csv"))

#MLE Estimator of beta
result=multiple_MLE(parameter)

write.table(result$beta,paste("MLE_beta.csv"))
write.table(result$sample_size,paste("MLE_sample_size.csv"))
