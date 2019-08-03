source("PowerLaw.R")
source("read_data.R")

#vector of city sizes in Germany
citySize=readVector("C:\\Users\\Ostrovski\\Google Drive\\Writing\\PowerLaw\\CitySize\\list_ge.csv")

kmin=50000
kmax=5000000
scale=10000
alpha=0.01

res=asymptotic_test(alpha,citySize,kmin,kmax,scale)

cdf=list2freq(citySize,kmin,kmax,scale)
cdf=cdf/sum(cdf)
cdf=cumsum(cdf)

pcdf=powerLawCDF(2.11180419921875,kmin/scale,kmax/scale)
distanceCDF(cdf,pcdf)
