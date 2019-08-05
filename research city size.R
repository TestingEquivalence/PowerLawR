source("PowerLaw.R")
source("read_data.R")
source("asymptotic_test.R")

#vector of city sizes in Germany
citySize=readVector("C:\\Users\\Ostrovski\\Google Drive\\Writing\\PowerLaw\\CitySize\\list_ge.csv")

kmins=c(10000,20000,50000,100000)
kmaxs=c(max(citySize),5000000,10000000,20000000)
scales=c(10000,1000,100,10,1)
alpha=0.05

result=multiple_asymptotic_test(alpha,citySize,kmins,kmaxs, 10000)


write.table(beta,"beta_100.txt")
write.table(distance,"distance_100.txt")
write.table(min_eps,"min_eps_100.txt")

