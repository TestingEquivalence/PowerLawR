source("PowerLaw.R")
source("tests.R")
source("power.R")
source("read_data.R")
source("simulation.R")
source("size.R")

#vector of city sizes in Germany
words=readVector("C:\\Users\\Ostrovski\\Google Drive\\Writing\\PowerLaw\\WordCounts\\words.csv")

#different k_min values
kmins=c(200,300,500,700,1000)

#different k_max values, first value coincides with largest city
kmaxs=c(max(words),2e6,5e6,10e6)

#scale for the populatuion measurement
#it is necessary for the computational feasibility.
scale=10
#significance level
alpha=0.05

#carry out multiple asymptotic tests for the power law
#an asymptotic test is computed for each combination of k_min and k_max 
#use all but one CPU cores
result=multiple_test(alpha,words,kmins,kmaxs, scale, bootstrap = FALSE, nSimulation = 1000)

#write test results
write.table(result$beta,paste("beta_",scale,".txt"))
write.table(result$distance,paste("distance_",scale,".txt"))
write.table(result$min_eps,paste("min_eps_",scale,".txt"))
write.table(result$sample_size,paste("sample_size_",scale,".txt"))


#compute test power at the power law points
alpha=0.05
kmin=20e3
kmax=10e6
scale=1e3
nSamples=1000
n=662

for (beta in c(2.1, 2.2, 2.3, 2.4, 2.5)) {
  size=sizeAtPowerLaw(n,kmin,kmax,scale,beta,nSamples,alpha,
                      bootstrap = FALSE, nSimulation = 1000)
  write.table(t(size), paste("size",beta,".txt"))
}

beta=2.3
eps=0.12

pw=boundaryPower(n,eps,kmin,kmax,scale,beta,alpha, boundaryPointType = 1,
                 bootstrap = FALSE, nSimulation = 1000)
write.table(pw, "powerLawStress.txt")

pw=boundaryPower(n,eps,kmin,kmax,scale,beta,alpha, boundaryPointType = 2,
                 bootstrap = FALSE, nSimulation = 1000)
write.table(pw, "uniformRandomStress.txt")
 
