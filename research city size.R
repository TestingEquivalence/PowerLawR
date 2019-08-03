source("PowerLaw.R")
source("read_data.R")
source("asymptotic_test.R")

#vector of city sizes in Germany
citySize=readVector("C:\\Users\\Ostrovski\\Google Drive\\Writing\\PowerLaw\\CitySize\\list_ge.csv")

kmins=c(10000,20000,50000,100000)
kmaxs=c(max(citySize),5000000,10000000,20000000)
scales=c(10000,1000,100,10,1)
alpha=0.05

# res=asymptotic_test(alpha,citySize,kmin,kmax,scale)
# 
# cdf=list2freq(citySize,kmin,kmax,scale)
# cdf=cdf/sum(cdf)
# cdf=cumsum(cdf)
# 
# pcdf=powerLawCDF(2.11180419921875,kmin/scale,kmax/scale)
# distanceCDF(cdf,pcdf)

nrow=length(kmins)
ncol = length(kmaxs)
min_eps=matrix(data=NA,nrow,ncol)
beta=matrix(data=NA,nrow,ncol)
distance=matrix(data=NA,nrow,ncol)

rownames(min_eps)=kmins
rownames(beta)=kmins
rownames(distance)=kmins

colnames(min_eps)=kmaxs
colnames(beta)=kmaxs
colnames(distance)=kmaxs

scale=100

for (i in c(1:nrow)) {
  for (j in c(1:ncol)) {
    res=asymptotic_test(alpha,citySize,kmins[i],kmax[j],scale)
    min_eps[i,j]=res[1]
    distance[i,j]=res[2]
    beta[i,j]=res[3]
  }
}

write.table(beta,"beta_100.txt")
write.table(distance,"distance_100.txt")
write.table(min_eps,"min_eps_100.txt")

