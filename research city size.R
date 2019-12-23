source("PowerLaw.R")
source("read_data.R")
source("asymptotic_test.R")
source("size.R")
source("simulation.R")
source("power.R")

#vector of city sizes in Germany
citySize=readVector("C:\\Users\\Ostrovski\\Google Drive\\Writing\\PowerLaw\\CitySize\\list_ge.csv")

#different k_min values
kmins=c(10e3,20e3,30e3,40e3,50e3)

#different k_max values, first value coincides with largest city
kmaxs=c(max(citySize),5e6,10e6,20e6)

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
kmin=20000
kmax=10000000
scale=10000
beta=2
nSamples=1000
n=662

size=sizeAtPowerLaw(n,kmin,kmax,scale,beta,nSamples,alpha)
write.table(size, "size.txt")

#compute test power at the boundary points,
#which are close to the given power law
alpha=0.05
kmin=20000
kmax=10000000
scale=10000
beta=2
nSamples=1000
n=662
eps=0.08

pw=boundaryPower(n,eps,kmin,kmax,scale,beta,alpha)
write.table(pw, "power.txt")
