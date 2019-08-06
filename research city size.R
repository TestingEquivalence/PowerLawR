source("PowerLaw.R")
source("read_data.R")
source("asymptotic_test.R")

#vector of city sizes in Germany
citySize=readVector("C:\\Users\\Ostrovski\\Google Drive\\Writing\\PowerLaw\\CitySize\\list_ge.csv")

kmins=c(10000,20000,50000,100000)
kmaxs=c(max(citySize),5000000,10000000,20000000)
scale=10000
alpha=0.05

result=multiple_asymptotic_test(alpha,citySize,kmins,kmaxs, scale)


write.table(result$beta,paste("beta_",scale,".txt"))
write.table(result$distance,paste("distance_",scale,".txt"))
write.table(result$min_eps,paste("min_eps_",scale,".txt"))

