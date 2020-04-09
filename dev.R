source("PowerLaw.R")
source("tests.R")
source("power.R")
source("read_data.R")
source("simulation.R")
source("size.R")
source("bootstrap_tests.R")

alpha=0.05
scale=10e3
nSamples=1000
n=662
beta=2.3
eps=0.10
adjEps=1
kmin=20e3/scale
kmax=10e6/scale


vbeta=c(1:100)
vpoints=list()

set.seed(01012020)
for (i in c(1:100)){
  point=uniformRandomStress(kmin,kmax,beta,eps)
  res=nearestPowerLaw(cumsum(point),kmin,kmax,1,3)
  vpoints[[i]]=point
  vbeta[[i]]=res$minimum
}

summary(vbeta)

k=46 #worst case
k=80 #best case 
pl=powerLawCDF(vbeta[k],kmin,kmax)
pcdf=cumsum(vpoints[[k]])
diff=pl-pcdf
plot(diff)

p=vpoints[[k]]
v_mineps=c(1:1000)
v_dst=c(1:1000)
v_beta=c(1:1000)

set.seed(01082019)

for (i in c(1:1000)){
  counting=rmultinom(n=1,size=n,prob=p)
  res=bootstrap_test(alpha,counting,kmin,kmax,scale=1,tol=0.001)
  v_mineps[i]=res[1]
  v_dst[i]=res[2]
  v_beta[i]=res[3]
}

summary(v_mineps)
summary(v_dst)
summary(v_beta)

hist(v_mineps)
hist(v_dst)
hist(v_beta)

v=v_mineps<eps
sum(v==TRUE)/1000

citySize=readVector("C:\\Users\\Ostrovski\\Google Drive\\Writing\\PowerLaw\\CitySize\\list_ge.csv")
alpha=0.05
kmin=20e3
kmax=10e6
scale=10e3
nSamples=1000
n=662
eps=0.08

frequency=list2freq(citySize,kmin,kmax,scale)
set.seed(30062020)
res=bootstrap_test2(alpha, frequency, kmin, kmax,
                    scale,nSimulation=1000, tol=0.001, nDirections=200)
