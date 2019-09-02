source("PowerLaw.R")
source("read_data.R")
source("asymptotic_test.R")
source("size.R")
source("simulation.R")

#vector of city sizes in Germany
citySize=readVector("C:\\Users\\Ostrovski\\Google Drive\\Writing\\PowerLaw\\CitySize\\list_ge.csv")

#different k_min values
kmins=c(10000,20000,30000,40000,50000)

#different k_max values, first value coincides with largest city
kmaxs=c(max(citySize),5000000,10000000,20000000)

#scale for the populatuion measurement
#it is necessary for the computational feasibility.
scale=100
#significance level
alpha=0.05

#carry out multiple asymptotic tests for the power law
#an asymptotic test is computed for each combination of k_min and k_max 
#use all but one CPU cores
result=multiple_asymptotic_test(alpha,citySize,kmins,kmaxs, scale)

#write test results
write.table(result$beta,paste("beta_",scale,".txt"))
write.table(result$distance,paste("distance_",scale,".txt"))
write.table(result$min_eps,paste("min_eps_",scale,".txt"))
write.table(result$sample_size,paste("sample_size_",scale,".txt"))


#compute test power at the power law points
alpha=0.05
kmin=10000
kmax=max(citySize)
scale=10000
beta=2
nSamples=1000
eps=0.2

power=testPowerAtPowerLaw(citySize,kmin,kmax,scale,beta,eps,nSamples,alpha)
